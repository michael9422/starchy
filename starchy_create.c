/*
 * starchy_create.c
 *
 * The "starchy" astronomy library is for coordinatizing pictures of stars
 * by comparing their positions with the positions of stars in a catalog.
 * It is an all-sky search, limited to small fields of view (~10 degrees
 * or less).
 *
 *-----------------------------------------------------------------------
 * Copyright (C) 2014 - 2016 FSF
 *
 * This library is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *-----------------------------------------------------------------------
 *
 * This program creates searchable files of stars and triplets using star
 * positions from a catalog(s). The files are:
 *
 *    file               description
 *  --------    ------------------------------------------
 *    star      general star info (relative format)
 *    utree     indexed by star direction vector (3-way kd-tree)
 *    member    star triplet memberships (singly linked lists)
 *    ttree     triplets (3-way kd-tree)
 *
 * As this program creates the index files, constraints are imposed on the
 * angular separation of stars in triplets. These depend on the field of view
 * of the pictures to be coordinatized. Smaller field-of-view pictures require
 * more stars and thus larger index files.
 *
 * Run this program using:
 *
 *	# starchy_create <faint> <sep_min> <sep_max> <member> <int_min>
 *
 * where the parameters are:
 *
 *    faint     - faint limit (V magnitude)
 *    sep_min   - minimum star spacing for triplets (degrees)
 *    sep_max   - maximum star spacing for triplets (degrees)
 *    member    - maximum number of triplets for any star to be a member of
 *    int_min   - interior angle minimum (degrees)
 *
 * Creating an index file from a star catalog can take lots of computer
 * time and memory. Extra contraints may be needed to limit the file size.
 * Note however, that reducing the number of triplets in the index file may
 * cause more failures to identify star fields.
 *
 * The proper motions, if available, are used to adjust the star positions
 * so that they are current at the time of file creation. Adjustments to
 * positions for parallax, precession, nutation, or anything else are not
 * applied here or for the searching routines.
 *
 * For using the Tycho-2 catalog, the default path names are
 *
 *      /usr/local/data/tycho/CD/data/index.dat
 *      /usr/local/data/tycho/CD/data/catalog.dat
 *
 *-----------------------------------------------------------------------
 *
 * File history:
 *
 *	2014-12-04 Begin writing the program
 *
 * To do:
 *
 */

#define _GNU_SOURCE	1

#include <ctype.h>
#include <errno.h>
#include <fcntl.h>
#include <float.h>
#include <inttypes.h>
#include <limits.h>
#include <math.h>
#include <search.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/errno.h>
#include <sys/file.h>
#include <sys/ioctl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/un.h>
#include <time.h>
#include <unistd.h>
#include <utmp.h>

#include <starchy.h>

STARCHY_FILE *sf;

float v_flimit;			// V magnitude faint limit
float sep_min;			// radians
float sep_max;			// radians
int nmax_member;		// star triplet membership limit
float min_interior_angle;	// for triangles formed by triplets (radians)

int n_ua;
struct starchy_u_struct *ua;
struct starchy_u_struct u_min, u_max;

int n_va;
struct starchy_v_struct *va;

int n_triplet;
struct starchy_t_struct *ta;


/*
 * Conversions for arcsec, arcmin, degree, radian
 */
#define SEC2DEG		2.777777778e-4
#define SEC2RAD		4.848136811e-6
#define DEG2RAD		0.01745329252
#define RAD2DEG		57.29577951
#define MIN2RAD		2.9088820e-4
#define RAD2MIN		3437.746771
#define RAD2SEC		206264.8062
#define DEG2MIN		60.0
#define MIN2HRS		0.00111111111
#define MIN2DEG		0.01666666667


//--------------------------------------------------- mag3_f
float mag3_f(float v[3])
{
	return (sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]));
}


//--------------------------------------------------- inner3_f
float inner3_f(float a[3], float b[3])
{
	return (a[0]*b[0] + a[1]*b[1] + a[2]*b[2]);
}


//--------------------------------------------------- cross3_f
void cross3_f(float a[3], float b[3], float c[3])
{
	a[0] =  (b[1]*c[2] - c[1]*b[2]);
	a[1] = -(b[0]*c[2] - c[0]*b[2]);
	a[2] =  (b[0]*c[1] - c[0]*b[1]);
}


//--------------------------------------------------- cross3
void cross3(double a[3], double b[3], double c[3])
{
	a[0] =  (b[1]*c[2] - c[1]*b[2]);
	a[1] = -(b[0]*c[2] - c[0]*b[2]);
	a[2] =  (b[0]*c[1] - c[0]*b[1]);
}


/*--------------------------------------- add_tycho_stars
 * Select stars to be used from the Tycho-2 catalog, and adds them to the star
 * file using the global variable sf->fd_star.
 *
 * Parameter 'v_flimit' is the V magnitude faint limit.
 * Parameter 'year' is used to adjust the positions for proper motion.
 * The return value is the number of records written.
 */
int add_tycho_stars(float v_flimit, float year)
{
	FILE *fp;
	char s[300];
	int32_t i, n;
	struct tycho_index_struct {
		int rn_cat, rn_uapp;
		double rad_min, rad_max, decd_min, decd_max;
	};
	struct tycho_index_struct tyi[9538];
	int n_tyi;
	struct tycho_cat_struct {
		int tyc1, tyc2, tyc3;
		double mrad, mdecd, pmra, pmdec, Bt, Vt;
	};
	struct tycho_cat_struct ty;
	struct starchy_star_struct ss0;
	int n_write;
	double rra, rdec;
	double u0[3], u1[3], u2[3];
	float f0, f1;

	n_write = 0;

	//---------- index.dat
	snprintf(s, sizeof(s), "%s/%s", starchy_tycho_path, starchy_tycho_index);
	fp = fopen(s, "r");
	if (fp == NULL) {
		fprintf(stderr, "%s: failed opening Tycho index\n", __func__);
		return n_write;
	}
	n_tyi = 0;
	while (fread(s, 44, 1, fp) == 1) {
		if (n_tyi >= 9538) break;

		memset(&tyi[n_tyi], 0, sizeof(struct tycho_index_struct));
		tyi[n_tyi].rn_cat = atoi(s + 0);
		tyi[n_tyi].rn_uapp = atoi(s + 8);
		tyi[n_tyi].rad_min = atof(s + 15);
		tyi[n_tyi].rad_max = atof(s + 22);
		tyi[n_tyi].decd_min = atof(s + 29);
		tyi[n_tyi].decd_max = atof(s + 36);
		n_tyi++;
	}
	fclose (fp);

	//----------- catalog.dat
	snprintf(s, sizeof(s), "%s/%s", starchy_tycho_path, starchy_tycho_catalog);
	fp = fopen(s, "r");
	if (fp == NULL) {
		fprintf(stderr, "%s: failed opening Tycho catalog\n", __func__);
		return n_write;
	}

	for (i = 0; i < n_tyi - 1; i++) {
		fseek(fp, 208*(tyi[i].rn_cat - 1), SEEK_SET);
		while (ftell(fp) < 208 * (tyi[i + 1].rn_cat - 1)) {
			memset(s, 0, sizeof(s));
			if (fread(s, 208, 1, fp) != 1) break;
			if (s[13] != ' ') continue;	// full entries only

			memset(&ty, 0, sizeof(ty));
			ty.Vt = atof(s + 123);
			if (ty.Vt > v_flimit) continue;

			ty.tyc1 = atoi(s + 0);
			ty.tyc2 = atoi(s + 5);
			ty.tyc3 = atoi(s + 11);
			ty.mrad = atof(s + 15);		// J2000 mean RA degrees
			ty.mdecd = atof(s + 28);	// J2000 mean Dec degrees
			ty.pmra = atof(s + 41);		// proper motion in mas/year
			ty.pmdec = atof(s + 49);	// proper motion in mas/year
			ty.Bt = atof(s + 110);

			if (ty.tyc1 != (i + 1)) {
				fprintf(stderr, "%s: error in catalog file", __func__);
				break;
			}

			if ((ty.pmra > 200.0) || (ty.pmdec > 200.0)) continue;
#if 0
			// I think that this calculation is incorrect
			ty.mrad += (ty.pmra / 1000.0) * SEC2DEG * (year - 2000.0);
			ty.mdecd += (ty.pmdec / 1000.0) * SEC2DEG * (year - 2000.0);

			rra = ty.mrad * DEG2RAD;
			rdec = ty.mdecd * DEG2RAD;
#else
			// this is better
			rra = ty.mrad * DEG2RAD;
			rdec = ty.mdecd * DEG2RAD;

			u0[0] = cos(rdec) * cos(rra);
			u0[1] = cos(rdec) * sin(rra);
			u0[2] = sin(rdec);

			u1[0] = sin(rdec) * -cos(rra);
			u1[1] = sin(rdec) * -sin(rra);
			u1[2] = cos(rdec);

			cross3(u2, u1, u0);

			f0 = (ty.pmdec / 1000.0) * SEC2RAD * (year - 2000.0);	// PM North
			f1 = (ty.pmra / 1000.0) * SEC2RAD * (year - 2000.0);	// PM East

			u0[0] += f0 * u1[0] + f1 * u2[0];
			u0[1] += f0 * u1[1] + f1 * u2[1];
			u0[2] += f0 * u1[2] + f1 * u2[2];

			rra = atan2(u0[1], u0[0]);
			rdec = atan2(u0[2], sqrt(u0[0] * u0[0] + u0[1] * u0[1]));

			if (isnan(rra) || isnan(rdec)) continue;

			ty.mrad = rra * RAD2DEG;
			ty.mdecd = rdec * RAD2DEG;

			if (ty.mrad < 0.0) ty.mrad += 360.0;
			if (fabs(ty.mdecd) >= 90.0) continue;
#endif
			memset(&ss0, 0, sizeof(ss0));
			ss0.rad = ty.mrad;
			ss0.decd = ty.mdecd;
			ss0.u[0] = cos(rdec) * cos(rra);
			ss0.u[1] = cos(rdec) * sin(rra);
			ss0.u[2] = sin(rdec);
			ss0.v_mag = ty.Vt;
			ss0.b_mag = ty.Bt;
			ss0.catalog = STARCHY_CATALOG_TYCHO;
			strncpy(ss0.catalog_star_id, s, 12);
			n = write(sf->fd_star, &ss0, sizeof(ss0));
			if (n == sizeof(ss0))
				 n_write++;
		}
	}

	fclose (fp);
	return n_write;
}


//--------------------------------------------------- ucomp0
int ucomp0(const void *a, const void *b)
{
	struct starchy_u_struct *p, *q;
	int r;

	p = (struct starchy_u_struct *) a;
	q = (struct starchy_u_struct *) b;

	if (p->u[0] == q->u[0])
		r = 0;
	else
		r = (p->u[0] > q->u[0]) ? +1 : -1;
	return r;
}

//--------------------------------------------------- ucomp1
int ucomp1(const void *a, const void *b)
{
	struct starchy_u_struct *p, *q;
	int r;

	p = (struct starchy_u_struct *) a;
	q = (struct starchy_u_struct *) b;

	if (p->u[1] == q->u[1])
		r = 0;
	else
		r = (p->u[1] > q->u[1]) ? +1 : -1;
	return r;
}

//--------------------------------------------------- ucomp2
int ucomp2(const void *a, const void *b)
{
	struct starchy_u_struct *p, *q;
	int r;

	p = (struct starchy_u_struct *) a;
	q = (struct starchy_u_struct *) b;

	if (p->u[2] == q->u[2])
		r = 0;
	else
		r = (p->u[2] > q->u[2]) ? +1 : -1;
	return r;
}

//--------------------------------------------------- ucomp_sp
int ucomp_sp(const void *a, const void *b)
{
	struct starchy_u_struct *p, *q;

	p = (struct starchy_u_struct *) a;
	q = (struct starchy_u_struct *) b;

	return (p->sp < q->sp) ? -1 : (p->sp > q->sp);
}

//--------------------------------------------------- vcomp_mag
int vcomp_mag(const void *a, const void *b)
{
	int r;
	struct starchy_v_struct *p, *q;

	p = (struct starchy_v_struct *) a;
	q = (struct starchy_v_struct *) b;

	if (p->v_mag == q->v_mag)
		r = 0;
	else
		r = (p->v_mag > q->v_mag) ? +1 : -1;
	return r;
}


/*--------------------------------------- kdtree_u_build
 * This function uses global variables ua and sf->fd_utree.
 * The array ua must contain the stars to be added. Their order is
 * rearranged as the node records are written to the file sf->fd_utree.
 */
off_t kdtree_u_build(int i_min, int n, int depth)
{
	off_t cp;
	int axis, m;
	struct starchy_unode_struct un0;

	if (n <= 0) return 0;

	cp = lseek(sf->fd_utree, 0, SEEK_END);
	axis = depth % 3;

	switch (axis) {
	case 0:
		qsort(&ua[i_min], n, sizeof(ua[0]), ucomp0);
		break;
	case 1:
		qsort(&ua[i_min], n, sizeof(ua[0]), ucomp1);
		break;
	case 2:
		qsort(&ua[i_min], n, sizeof(ua[0]), ucomp2);
		break;
	}

	m = i_min + n/2;	// median
	un0.u[0] = ua[m].u[0];
	un0.u[1] = ua[m].u[1];
	un0.u[2] = ua[m].u[2];
        un0.sp = ua[m].sp;
	un0.l = 0;
	un0.r = 0;
	write(sf->fd_utree, &un0, sizeof(un0));

	un0.l = kdtree_u_build(i_min, m - i_min, depth + 1);
	un0.r = kdtree_u_build(m + 1, n - (m - i_min) - 1, depth + 1);

	lseek(sf->fd_utree, cp, SEEK_SET);
	write(sf->fd_utree, &un0, sizeof(un0));

	return cp;
}


/*--------------------------------------- kdtree_u_search
 * This function uses global variables ua, n_ua and sf->fd_utree.
 * The array ua is filled with the selected stars.
 */
void kdtree_u_search(off_t np, int depth)
{
	int axis, n;
	struct starchy_unode_struct un0;
	bool bool_l, bool_r;

	if (np == 0) return;

	lseek(sf->fd_utree, np, SEEK_SET);
	n = read(sf->fd_utree, &un0, sizeof(un0));
	if (n != sizeof(un0)) {
		fprintf(stderr, "%s: read utree file failed\n", __func__);
		return;
	}

	if ((un0.u[0] > u_min.u[0]) && (un0.u[0] < u_max.u[0]) &&
	    (un0.u[1] > u_min.u[1]) && (un0.u[1] < u_max.u[1]) &&
	    (un0.u[2] > u_min.u[2]) && (un0.u[2] < u_max.u[2])) {
		ua[n_ua].u[0] = un0.u[0];
		ua[n_ua].u[1] = un0.u[1];
		ua[n_ua].u[2] = un0.u[2];
		ua[n_ua].sp = un0.sp;
		n_ua++;
	}
	axis = depth % 3;

	switch (axis) {
	case 0:
		bool_l = (un0.u[0] > u_min.u[0]);
		bool_r = (un0.u[0] < u_max.u[0]);
		break;
	case 1:
		bool_l = (un0.u[1] > u_min.u[1]);
		bool_r = (un0.u[1] < u_max.u[1]);
		break;
	case 2:
		bool_l = (un0.u[2] > u_min.u[2]);
		bool_r = (un0.u[2] < u_max.u[2]);
		break;
	}
	if (bool_l) kdtree_u_search(un0.l, depth + 1);
	if (bool_r) kdtree_u_search(un0.r, depth + 1);

	return;
}


/*-------------------------------------- add_member
 * Insert a member record into the beginning of the list using the global
 * variables sf->fd_member and sf->fd_star.
 */
void add_member(struct starchy_star_struct *a, off_t asp, off_t bsp, off_t csp)
{
	struct starchy_member_struct m;
	off_t mp;

	m.bsp = bsp;
	m.csp = csp;
	m.next = a->member_head;

	mp = lseek(sf->fd_member, 0, SEEK_END);
	write(sf->fd_member, &m, sizeof(m));

	a->member_head = mp;
	a->n_member++;

	lseek(sf->fd_star, asp, SEEK_SET);
	write(sf->fd_star, a, sizeof(*a));

	return;
}


//--------------------------------------------------- tcomp0
int tcomp0(const void *a, const void *b)
{
	struct starchy_t_struct *p, *q;
	int r;

	p = (struct starchy_t_struct *) a;
	q = (struct starchy_t_struct *) b;

	if (p->ab == q->ab)
		r = 0;
	else
		r = (p->ab > q->ab) ? +1 : -1;
	return r;
}

//--------------------------------------------------- tcomp1
int tcomp1(const void *a, const void *b)
{
	struct starchy_t_struct *p, *q;
	int r;

	p = (struct starchy_t_struct *) a;
	q = (struct starchy_t_struct *) b;

	if (p->a == q->a)
		r = 0;
	else
		r = (p->a > q->a) ? +1 : -1;
	return r;
}

//--------------------------------------------------- tcomp2
int tcomp2(const void *a, const void *b)
{
	struct starchy_t_struct *p, *q;
	int r;

	p = (struct starchy_t_struct *) a;
	q = (struct starchy_t_struct *) b;

	if (p->b == q->b)
		r = 0;
	else
		r = (p->b > q->b) ? +1 : -1;
	return r;
}



/*--------------------------------------- kdtree_t_build
 * This function uses global variables ta and sf->fd_ttree.
 * The array ta must contain the triplets to be added. Their order is
 * rearranged as the node records are written to the file sf->fd_ttree.
 */
off_t kdtree_t_build(int i_min, int n, int depth)
{
	off_t cp;
	int axis, m;
	struct starchy_tnode_struct tn0;

	if (n <= 0) return 0;

	cp = lseek(sf->fd_ttree, 0, SEEK_END);
	axis = depth % 3;

	switch (axis) {
	case 0:
		qsort(&ta[i_min], n, sizeof(ta[0]), tcomp0);
		break;
	case 1:
		qsort(&ta[i_min], n, sizeof(ta[0]), tcomp1);
		break;
	case 2:
		qsort(&ta[i_min], n, sizeof(ta[0]), tcomp2);
		break;
	}

	m = i_min + n/2;	// median
	tn0.asp = ta[m].asp;
	tn0.bsp = ta[m].bsp;
	tn0.csp = ta[m].csp;
	tn0.ab = ta[m].ab;
	tn0.a = ta[m].a;
	tn0.b = ta[m].b;
	tn0.l = 0;
	tn0.r = 0;
	write(sf->fd_ttree, &tn0, sizeof(tn0));

	tn0.l = kdtree_t_build(i_min, m - i_min, depth + 1);
	tn0.r = kdtree_t_build(m + 1, n - (m - i_min) - 1, depth + 1);

	lseek(sf->fd_ttree, cp, SEEK_SET);
	write(sf->fd_ttree, &tn0, sizeof(tn0));

	return cp;
}


//---------------------------------------------------- main
int main(int argc, char *argv[], char *env[])
{
	int i, j, k, n;
	float f0;
	double d0, d1, d2, d3, d4, d5;
	struct starchy_star_hdr_struct sh;
	struct starchy_star_struct ss0, ss1, ss2, ss3;
	struct starchy_utree_hdr_struct uh;
	struct starchy_member_hdr_struct mh;
	struct starchy_t_struct triplet0;
	struct starchy_ttree_hdr_struct th;
	time_t time0;
	struct tm tm0;
	off_t sp0, sp1, sp2;
	float v0[3], v1[3], v2[3];
	FILE *fp;
	char fn_tmp[300];

	/*
	 * Parse command line parameters
	 */
	f0 = 8.0;
	if (argc > 1) sscanf(argv[1], "%f", &f0);
	v_flimit = f0;

	f0 = 1.0;
	if (argc > 2) sscanf(argv[2], "%f", &f0);
	sep_min = f0 * DEG2RAD;

	f0 = 4.0;
	if (argc > 3) sscanf(argv[3], "%f", &f0);
	sep_max = f0 * DEG2RAD;

	nmax_member = 12;
	if (argc > 4) sscanf(argv[4], "%d", &nmax_member);

	f0 = 20.0;
	if (argc > 5) sscanf(argv[5], "%f", &f0);
	min_interior_angle = f0 * DEG2RAD;

	//------ open starchy files
	sf = starchy_fopen(true);
	if (sf == NULL) {
		fprintf(stderr, "%s: starchy_fopen failed\n", __func__);
		return -1;
	}

	//--- truncate the starchy index files to zero length
	ftruncate(sf->fd_star, 0);
	ftruncate(sf->fd_utree, 0);
	ftruncate(sf->fd_member, 0);
	ftruncate(sf->fd_ttree, 0);

	/*
	 * Write the star file header
	 */
	memset(&sh, 0, sizeof(sh));
	memcpy(sh.magic, "STARCHY0", 8);
	time0 = time(NULL);
	gmtime_r(&time0, &tm0);
	snprintf(sh.date_s, sizeof(sh.date_s), "%4d-%02d-%02dT%02d:%02d:%02d",
		tm0.tm_year + 1900,
		tm0.tm_mon + 1,
		tm0.tm_mday,
		tm0.tm_hour,
		tm0.tm_min,
		tm0.tm_sec);
	sh.l_rec = sizeof(struct starchy_star_struct);
	sh.n_rec = 0;

	lseek(sf->fd_star, 0, SEEK_SET);
	write(sf->fd_star, &sh, sizeof(sh));

	/*
	 * Add stars from the catalog(s) to the star file.
	 * If different catalogs are used in combination, the coordinate
	 * systems must agree (i.e., all J2000). The results from matching
	 * are based on the catalog coordinates.
	 */
	n = 0;
	f0 = (tm0.tm_year + 1900) + (tm0.tm_mon / 12.0);
	printf("%.1f year\n", f0);
	n += add_tycho_stars(v_flimit, f0);
//	n += add_gsc_stars(v_flimit, f0);
//	n += add_usno_stars(v_flimit, f0);

	lseek(sf->fd_star, 0, SEEK_SET);
	read(sf->fd_star, &sh, sizeof(sh));
	sh.n_rec = n;
	lseek(sf->fd_star, 0, SEEK_SET);
	write(sf->fd_star, &sh, sizeof(sh));
	printf("%d stars\n", sh.n_rec);

	/*
	 * Allocate an array for sorting the stars. This is needed for building
         * the kd-tree 'utree'
	 */
	n_ua = sh.n_rec;
	ua = malloc(n_ua * sizeof(struct starchy_u_struct));
	if (ua == NULL) {
		fprintf(stderr, "%s: malloc ua failed\n", __func__);
		return -1;
	}

	for (i = 0; i < n_ua; i++) {
		ua[i].sp = lseek(sf->fd_star, 0, SEEK_CUR);
		read(sf->fd_star, &ss0, sizeof(ss0));
		ua[i].u[0] = ss0.u[0];
		ua[i].u[1] = ss0.u[1];
		ua[i].u[2] = ss0.u[2];
	}

	/*
	 * Write the utree file header
	 */
	memset(&uh, 0, sizeof(uh));
	memcpy(uh.magic, "STARCHY1", 8);
	strncpy(uh.date_s, sh.date_s, sizeof(uh.date_s));
	uh.l_rec = sizeof(struct starchy_unode_struct);
	uh.root = 0;

	lseek(sf->fd_utree, 0, SEEK_SET);
	write(sf->fd_utree, &uh, sizeof(uh));

	/*
	 * Build a 3-way kd-tree for nearest neighbor searches, indexed by
	 * the components of u, the unit direction vector for the stars.
	 */
	uh.root = kdtree_u_build(0, n_ua, 0);

	lseek(sf->fd_utree, 0, SEEK_SET);
	write(sf->fd_utree, &uh, sizeof(uh));

	/*
	 * Write the member file header
	 */
	memset(&mh, 0, sizeof(mh));
	memcpy(mh.magic, "STARCHY2", 8);
	strncpy(mh.date_s, sh.date_s, sizeof(mh.date_s));
	mh.l_rec = sizeof(struct starchy_member_struct);

	lseek(sf->fd_member, 0, SEEK_SET);
	write(sf->fd_member, &mh, sizeof(mh));

	/*
	 * Select triplets of stars and write them to a temporary file
	 */
	snprintf(fn_tmp, sizeof(fn_tmp), "%s/%s", starchy_path_index, starchy_fname_tmp);
	fp = fopen(fn_tmp, "w+");
	if (fp == NULL) {
		fprintf(stderr, "open temporary triplet file failed\n");
		return -1;
	}

	lseek(sf->fd_star, 0, SEEK_SET);
	read(sf->fd_star, &sh, sizeof(sh));
	sp0 = lseek(sf->fd_star, 0, SEEK_CUR);

	n_triplet = 0;
	for (i = 0; i < sh.n_rec; i++) {		// vertex A star
		lseek(sf->fd_star, sp0, SEEK_SET);
		n = read(sf->fd_star, &ss0, sizeof(ss0));
		if (n != sizeof(ss0)) break;

		/*
		 * Search the utree for nearby stars
		 */
		n_ua = 0;
		d0 = sin(sep_max);

		u_min.u[0] = ss0.u[0] - d0;
		u_min.u[1] = ss0.u[1] - d0;
		u_min.u[2] = ss0.u[2] - d0;

		u_max.u[0] = ss0.u[0] + d0;
		u_max.u[1] = ss0.u[1] + d0;
		u_max.u[2] = ss0.u[2] + d0;

		kdtree_u_search(uh.root, 0);

		/*
		 * Eliminate stars that are too far away
		 */
		for (j = n_ua - 1; j >= 0; j--) {
			d0 = acos(inner3_f(ss0.u, ua[j].u));
			if (d0 < sep_max) continue;

			// eliminate it
			for (k = j; k < n_ua - 1; k++)
				memcpy(&ua[k], &ua[k + 1], sizeof(ua[0]));
			n_ua--;
		}
//		printf("%d stars selected\n", n_ua);

		/*
		 * Sort them so that the brightest stars are used first in
		 * triplets.
		 */
		n_va = n_ua;
		va = malloc(n_va * sizeof(va[0]));
		if (va == NULL) {
			fprintf(stderr, "%s: malloc va failed\n", __func__);
			return -1;
		}
		for (i = 0; i < n_va; i++) {
			va[i].u[0] = ua[i].u[0];
			va[i].u[1] = ua[i].u[1];
			va[i].u[2] = ua[i].u[2];
			va[i].sp = ua[i].sp;

			lseek(sf->fd_star, va[i].sp, SEEK_SET);
			n = read(sf->fd_star, &ss3, sizeof(ss3));

			va[i].v_mag = ss3.v_mag;
		}
		qsort(va, n_va, sizeof(va[0]), vcomp_mag);
		for (i = 0; i < n_va; i++) {
			ua[i].u[0] = va[i].u[0];
			ua[i].u[1] = va[i].u[1];
			ua[i].u[2] = va[i].u[2];
			ua[i].sp = va[i].sp;
		}
		free(va);

		/*
		 * There might be some duplication of triplets, e.g. triplet
		 * ABC and CAB or BCA are both in the kd-tree file.
		 * An attempt to prevent duplicates by imposing the constraint
		 * on the star order {sp(A) < sp(B) < sp(C)} in the star file
		 * resulted in identification failures, and so was abandoned.
		 */
		for (j = 0; j < n_ua; j++) {		// vertex B star
			// The separation angle of stars A & B
			d0 = acos(inner3_f(ss0.u, ua[j].u));
			if (d0 < sep_min) continue;

			sp1 = ua[j].sp;
			lseek(sf->fd_star, sp1, SEEK_SET);
			read(sf->fd_star, &ss1, sizeof(ss1));
			if (ss1.n_member >= nmax_member) continue;
			if (ss1.n_b >= nmax_member / 3) continue;

			for (k = 0; k < n_ua; k++) {	// vertex C star
				if (k == j) continue;

				cross3_f(v0, ss0.u, ua[j].u);
				cross3_f(v1, ss0.u, ua[k].u);
				d4 = mag3_f(v0);
				d5 = mag3_f(v1);
				v0[0] /= d4;
				v0[1] /= d4;
				v0[2] /= d4;

				v1[0] /= d5;
				v1[1] /= d5;
				v1[2] /= d5;

				// avoid triplets that reorder ABC to ACB
				cross3_f(v2, v0, v1);
				if (inner3_f(v2, ss0.u) < 0.0) continue;

				// interior angle for vertex A
				d1 = acos(inner3_f(v0, v1));
				if (d1 < min_interior_angle) continue;
//				if (d1 > 90.0 * DEG2RAD) continue;

				cross3_f(v0, ua[j].u, ss0.u);
				cross3_f(v1, ua[j].u, ua[k].u);
				d4 = mag3_f(v0);
				d5 = mag3_f(v1);
				v0[0] /= d4;
				v0[1] /= d4;
				v0[2] /= d4;

				v1[0] /= d5;
				v1[1] /= d5;
				v1[2] /= d5;
				// interior angle for vertex B
				d2 = acos(inner3_f(v0, v1));
				if (d2 < min_interior_angle) continue;

				d3 = PI - (d1 + d2);
				if (d3 < min_interior_angle) continue;

				sp2 = ua[k].sp;
				lseek(sf->fd_star, sp2, SEEK_SET);
				read(sf->fd_star, &ss2, sizeof(ss2));
				if (ss2.n_member >= nmax_member) continue;
				if (ss2.n_c >= nmax_member / 3) continue;

				if (isnan(d0) || isnan(d1) || isnan(d2))
					continue;
				if (isinf(d0) || isinf(d1) || isinf(d2))
					continue;

				/*
				 * Include this triplet
				 */
				triplet0.asp = sp0;
				triplet0.bsp = sp1;
				triplet0.csp = sp2;
				triplet0.ab = d0;
				triplet0.a = d1;
				triplet0.b = d2;
				fwrite(&triplet0, sizeof(triplet0), 1, fp);

				ss0.n_a++;
				ss1.n_b++;
				ss2.n_c++;

				add_member(&ss0, sp0, sp1, sp2);
				add_member(&ss1, sp1, sp0, sp2);
				add_member(&ss2, sp2, sp0, sp1);

				n_triplet++;

//				if (ss1.n_member >= nmax_member) break;
				break;
			}
			if (ss0.n_member >= nmax_member) break;
			if (ss0.n_a >= nmax_member / 3) break;
		}
		sp0 += sizeof(ss0);
	}
	printf("%d triplets\n", n_triplet);

	/*
	 * Put the triplets from the temporary file into an array
	 */
	ta = malloc(n_triplet * sizeof(struct starchy_t_struct));
	if (ta == NULL) {
		fprintf(stderr, "%s: malloc ta failed\n", __func__);
		return -1;
	}
	rewind(fp);
	fread(ta, sizeof(triplet0), n_triplet, fp);
	fclose(fp);
	unlink(fn_tmp);

	/*
	 * Write the ttree (triplet kd-tree) file header
	 */
	memset(&th, 0, sizeof(th));
	memcpy(th.magic, "STARCHY3", 8);
	strncpy(th.date_s, sh.date_s, sizeof(th.date_s));
	th.l_rec = sizeof(struct starchy_tnode_struct);
	th.min_interior_angle = min_interior_angle;
	th.root = 0;

	lseek(sf->fd_ttree, 0, SEEK_SET);
	write(sf->fd_ttree, &th, sizeof(th));

	/*
	 * Build a 3-way kd-tree for nearest neighbor searches, indexed by
	 * three angles from the star triplets.
	 */
	th.root = kdtree_t_build(0, n_triplet, 0);

	lseek(sf->fd_ttree, 0, SEEK_SET);
	write(sf->fd_ttree, &th, sizeof(th));


	free(ta);
	free(ua);
	starchy_fclose(sf);

	return 0;
}
