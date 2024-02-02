/*
 * starchy.c
 *
 * The "starchy" astronomy library is for coordinatizing pictures of stars
 * by comparing their positions with the positions of stars in a catalog.
 * It is an all-sky search, limited to small fields of view (~10 degrees or
 * less).
 *
 *-----------------------------------------------------------------------
 * Copyright (C) 2014 - 2015 FSF
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
 * These are the library routines.
 * All exported symbols (for linking) must have a "starchy_" prefix.
 * Otherwise, declare functions (or global variables) as static.
 *
 *-----------------------------------------------------------------------
 *
 * File history:
 *
 *	2015-04-23 Put error messages in a buffer instead of printing them.
 *
 *	2014-12-09 Begin writing the program
 *
 * To do:
 *
 *	- The calculations are not precise. For example, the functions
 *	  starchy_image2radec and starchy_radec2image should be inverse
 *	  functions exactly, but they are not. For my purpose purpose though,
 *	  it is good enough.
 *
 *	- Moving the global variables into the starchy_struct would allow
 *	  multi-threaded usage. Currently, threaded usage requires semaphore
 *	  locking for exclusive access.
 *
 *	- Provide an alternative matching function, such as "starchy_match_xy",
 *	  which takes a list of pixel coordinates instead of an image buffer.
 *	  Or, use the same function, but test the image buffer pointer for NULL,
 *	  and also test if the starchy_struct member n_ais is non-zero.
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

/*
 * File names for the starchy files.
 */
char *starchy_path_index	= "/usr/local/share";	// "." for local
char *starchy_fname_lock	= "starchy.lock";
char *starchy_fname_tmp		= "starchy_index.tmp";
char *starchy_fname_star	= "starchy_index.star";
char *starchy_fname_utree	= "starchy_index.utree";
char *starchy_fname_member	= "starchy_index.member";
char *starchy_fname_ttree	= "starchy_index.ttree";

/*
 * Path to the Tycho-2 catalog
 */
char *starchy_tycho_path	= "/usr/local/data";
char *starchy_tycho_index	= "tycho/CD/data/index.dat";
char *starchy_tycho_catalog	= "tycho/CD/data/catalog.dat";

/*
 * Other catalogs (not implemented)
 */
//char *ua2_path = "/usr/local/data";     // USNO A2 catalog path
//char *gsc_north = "http://archive.eso.org/skycat/servers/gsc-server";

/*
 * Global error message buffer
 */
char starchy_err_msg[80];


#define NMAX_TMATCH	3

static bool relaxed_fit;
static struct starchy_tmatch_struct tmatch[NMAX_TMATCH];
static int n_tmatch;

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


/*--------------------------------------- starchy_fopen
 * This function opens a lock file and the four index data files.
 */
STARCHY_FILE *starchy_fopen(bool create)
{
	STARCHY_FILE *p;
	int oflags, pflags;
	int n, fd;

	memset(starchy_err_msg, 0, sizeof(starchy_err_msg));

	//------ allocate starchy file structure
	p = malloc(sizeof(*p));
	if (p == NULL) {
		snprintf(starchy_err_msg, sizeof(starchy_err_msg), "%s: malloc starchy_fp failed", __func__);
		goto fail_0;
	}
	memset(p, 0, sizeof(*p));

	//----- allocate file names
	n = strlen(starchy_path_index) + strlen(starchy_fname_lock);
	p->fname_lock = malloc(n + 10);
	snprintf(p->fname_lock, n + 5, "%s/%s", starchy_path_index, starchy_fname_lock);

	n = strlen(starchy_path_index) + strlen(starchy_fname_star);
	p->fname_star = malloc(n + 10);
	snprintf(p->fname_star, n + 5, "%s/%s", starchy_path_index, starchy_fname_star);

	n = strlen(starchy_path_index) + strlen(starchy_fname_utree);
	p->fname_utree = malloc(n + 10);
	snprintf(p->fname_utree, n + 5, "%s/%s", starchy_path_index, starchy_fname_utree);

	n = strlen(starchy_path_index) + strlen(starchy_fname_member);
	p->fname_member = malloc(n + 10);
	snprintf(p->fname_member, n + 5, "%s/%s", starchy_path_index, starchy_fname_member);

	n = strlen(starchy_path_index) + strlen(starchy_fname_ttree);
	p->fname_ttree = malloc(n + 10);
	snprintf(p->fname_ttree, n + 5, "%s/%s", starchy_path_index, starchy_fname_ttree);

	if ((p->fname_lock == NULL) ||
	    (p->fname_star == NULL) ||
	    (p->fname_utree == NULL) ||
	    (p->fname_member == NULL) ||
	    (p->fname_ttree == NULL)) {
		snprintf(starchy_err_msg, sizeof(starchy_err_msg), "%s: malloc fnames failed", __func__);
		goto fail_1;
	}

	//------ open data files
	if (create)
		oflags = O_RDWR | O_CREAT | O_TRUNC;
	else
		oflags = O_RDWR;
	pflags = S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH;

	//------ lock for exclusive access
	fd = open(p->fname_lock, O_RDWR | O_CREAT, pflags);
	p->fd_lock = fd;
	if (fd < 0) {
		snprintf(starchy_err_msg, sizeof(starchy_err_msg), "%s: open lock file failed", __func__);
		goto fail_1;
	}
	if (flock(fd, LOCK_EX /* | LOCK_NB */) != 0) {
		snprintf(starchy_err_msg, sizeof(starchy_err_msg), "%s: flock failed", __func__);
		close(fd);
		goto fail_1;
	}

	p->fd_star = open(p->fname_star, oflags, pflags);
	p->fd_utree = open(p->fname_utree, oflags, pflags);
	p->fd_member = open(p->fname_member, oflags, pflags);
	p->fd_ttree = open(p->fname_ttree, oflags, pflags);

	if ((p->fd_star < 0) || (p->fd_utree < 0) || (p->fd_member < 0) ||
							 (p->fd_ttree < 0)) {
		snprintf(starchy_err_msg, sizeof(starchy_err_msg),  "%s: failed opening data files", __func__);
		goto fail_2;
	}

	return p;		// success

fail_2:
	if (p->fd_ttree >= 0) close(p->fd_ttree);
	if (p->fd_member >= 0) close(p->fd_member);
	if (p->fd_utree >= 0) close(p->fd_utree);
	if (p->fd_star >= 0) close(p->fd_star);

	flock(p->fd_lock, LOCK_UN);
	close(p->fd_lock);
fail_1:
	if (p->fname_ttree != NULL) free(p->fname_ttree);
	if (p->fname_member != NULL) free(p->fname_member);
	if (p->fname_utree != NULL) free(p->fname_utree);
	if (p->fname_star != NULL) free(p->fname_star);
	if (p->fname_lock != NULL) free(p->fname_lock);

	free(p);
fail_0:
	return NULL;
}

//--------------------------------------- starchy_fclose
void starchy_fclose(STARCHY_FILE *p)
{
	close(p->fd_ttree);
	close(p->fd_member);
	close(p->fd_utree);
	close(p->fd_star);

	flock(p->fd_lock, LOCK_UN);
	close(p->fd_lock);

	free(p->fname_ttree);
	free(p->fname_member);
	free(p->fname_utree);
	free(p->fname_star);
	free(p->fname_lock);

	free(p);
	return;
}


//--------------------------------------------------- mag2_f
static float mag2_f(float w[2])
{
	return (sqrt(w[0]*w[0] + w[1]*w[1]));
}


//--------------------------------------------------- inner2_f
static float inner2_f(float a[2], float b[2])
{
	return (a[0]*b[0] + a[1]*b[1]);
}


//--------------------------------------------------- inner3
static double inner3(double a[3], double b[3])
{
	return (a[0]*b[0] + a[1]*b[1] + a[2]*b[2]);
}


//--------------------------------------------------- inner3_f
static float inner3_f(float a[3], float b[3])
{
	return (a[0]*b[0] + a[1]*b[1] + a[2]*b[2]);
}


//--------------------------------------------------- mag3
static double mag3(double v[3])
{
	return (sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]));
}


//--------------------------------------------------- cross3
static void cross3(double a[3], double b[3], double c[3])
{
	a[0] =  (b[1]*c[2] - c[1]*b[2]);
	a[1] = -(b[0]*c[2] - c[0]*b[2]);
	a[2] =  (b[0]*c[1] - c[0]*b[1]);
}


/*--------------------------------------------------- rsep_rorn
 * Given one stars coordinates rra1, rdec1, and another stars
 * coordinates rra2, rdec2, find the angular separation rsep, and angular
 * orientation rorn, where for small positive RA, rorn = 0; and for
 * positive DEC, rorn = PI/2. All angles are in radians.
 */
static void rsep_rorn(double rra1, double rdec1, double rra2, double rdec2,
                double *sep, double *orn, double *x, double *y)
{
   double rsep, rorn, d1, d2;
   double v1[3], v2[3], u1[3], u2[3];

   v1[0] = cos(rdec1) * cos(rra1);
   v1[1] = cos(rdec1) * sin(rra1);
   v1[2] = sin(rdec1);

   v2[0] = cos(rdec2) * cos(rra2);
   v2[1] = cos(rdec2) * sin(rra2);
   v2[2] = sin(rdec2);

   d1 = inner3(v1, v2);
   if (d1 > 1.0)
      rsep = 0.0;
   else if (d1 < -1.0)
      rsep = PI;
   else
      rsep = acos(d1);		// radian separation

   u1[0] = -sin(rra1);
   u1[1] = cos(rra1);
   u1[2] = 0.0;

   u2[0] = -sin(rdec1) * cos(rra1);
   u2[1] = -sin(rdec1) * sin(rra1);
   u2[2] = cos(rdec1);

   d1 = inner3(v2, u1);
   d2 = inner3(v2, u2);

   if ((d1 != 0.0) || (d2 != 0.0)) rorn = atan2(d2, d1);
   else rorn = 0.0;

   *sep = rsep;
   *orn = rorn;

   *x = rsep * cos(rorn);
   *y = rsep * sin(rorn);
}


/*--------------------------------------- starchy_kdtree_u_search
 * This fills the array p->ua with the selected stars.
 */
void starchy_kdtree_u_search(struct starchy_struct *p, off_t np, int depth)
{
	int axis, n;
	struct starchy_unode_struct un0;
	bool bool_l, bool_r;

	if (np == 0) return;

	lseek(p->sf->fd_utree, np, SEEK_SET);
	n = read(p->sf->fd_utree, &un0, sizeof(un0));
	if (n != sizeof(un0)) {
		snprintf(starchy_err_msg, sizeof(starchy_err_msg), "%s: read utree file failed", __func__);
		return;
	}

	if ((un0.u[0] > p->u_min.u[0]) && (un0.u[0] < p->u_max.u[0]) &&
	    (un0.u[1] > p->u_min.u[1]) && (un0.u[1] < p->u_max.u[1]) &&
	    (un0.u[2] > p->u_min.u[2]) && (un0.u[2] < p->u_max.u[2])) {
		if (p->n_ua < p->nmax_ua) {
			p->ua[p->n_ua].u[0] = un0.u[0];
			p->ua[p->n_ua].u[1] = un0.u[1];
			p->ua[p->n_ua].u[2] = un0.u[2];
			p->ua[p->n_ua].sp = un0.sp;
			p->n_ua++;
		}
	}
	axis = depth % 3;

	switch (axis) {
	case 0:
		bool_l = (un0.u[0] > p->u_min.u[0]);
		bool_r = (un0.u[0] < p->u_max.u[0]);
		break;
	case 1:
		bool_l = (un0.u[1] > p->u_min.u[1]);
		bool_r = (un0.u[1] < p->u_max.u[1]);
		break;
	case 2:
		bool_l = (un0.u[2] > p->u_min.u[2]);
		bool_r = (un0.u[2] < p->u_max.u[2]);
		break;
	}
	if (bool_l) starchy_kdtree_u_search(p, un0.l, depth + 1);
	if (bool_r) starchy_kdtree_u_search(p, un0.r, depth + 1);

	return;
}


/*--------------------------------------- kdtree_t_search
 * This fills the array p->ta with the selected stars.
 */
static void kdtree_t_search(struct starchy_struct *p, off_t np, int depth)
{
	int axis, n;
	struct starchy_tnode_struct tn0;
	bool bool_l, bool_r;

	if (np == 0) return;

	lseek(p->sf->fd_ttree, np, SEEK_SET);
	n = read(p->sf->fd_ttree, &tn0, sizeof(tn0));
	if (n != sizeof(tn0)) {
		snprintf(starchy_err_msg, sizeof(starchy_err_msg), "%s: read ttree file failed", __func__);
		return;
	}

	if ((tn0.ab > p->t_min.ab) && (tn0.ab < p->t_max.ab) &&
	    (tn0.a > p->t_min.a) && (tn0.a < p->t_max.a) &&
	    (tn0.b > p->t_min.b) && (tn0.b < p->t_max.b)) {
		if (p->n_ta < p->nmax_ta) {
			p->ta[p->n_ta].asp = tn0.asp;
			p->ta[p->n_ta].bsp = tn0.bsp;
			p->ta[p->n_ta].csp = tn0.csp;
			p->ta[p->n_ta].ab = tn0.ab;
			p->ta[p->n_ta].a = tn0.a;
			p->ta[p->n_ta].b = tn0.b;
			p->n_ta++;
		}
	}
	axis = depth % 3;

	switch (axis) {
	case 0:
		bool_l = (tn0.ab > p->t_min.ab);
		bool_r = (tn0.ab < p->t_max.ab);
		break;
	case 1:
		bool_l = (tn0.a > p->t_min.a);
		bool_r = (tn0.a < p->t_max.a);
		break;
	case 2:
		bool_l = (tn0.b > p->t_min.b);
		bool_r = (tn0.b < p->t_max.b);
		break;
	}
	if (bool_l) kdtree_t_search(p, tn0.l, depth + 1);
	if (bool_r) kdtree_t_search(p, tn0.r, depth + 1);

	return;
}


//--------------------------------------------------- rowcomp
static int rowcomp(const void *a, const void *b)
{
	uint16_t *p, *q;

	p = (uint16_t *) a;
	q = (uint16_t *) b;
	return (*q < *p) ? -1 : (*q > *p);
}


/*--------------------------------------------------- aicomp
 * Largest to smallest.
 */
static int aicomp(const void *a, const void *b)
{
	int *p, *q;

	p = (int32_t *) a;
	q = (int32_t *) b;
	return (*q < *p) ? -1 : (*q > *p);
}


//--------------------------------------------------- aiscomp
static int aiscomp(const void *a, const void *b)
{
	struct starchy_xy_struct *p, *q;

	p = (struct starchy_xy_struct *) a;
	q = (struct starchy_xy_struct *) b;
	return (q->pixsum < p->pixsum) ? -1 : (q->pixsum > p->pixsum);
}



/*------------------------------- starchy_find_image_stars
 * Find stars in the image buffer b and store them in array p->ais.
 * The image buffer b is clobbered.
 */
int starchy_find_image_stars(struct starchy_struct *p, int x, int y, uint16_t *b)
{
	uint16_t *row;
	int32_t i, j, k, l, m, n, i0, j0;
	double xcen, ycen, d0;
	int32_t pixsum, star_hw, xicen, yicen;
	int ai[100], ai2[100];


	p->n_ais = 0;

	row = malloc(sizeof(row[0]) * x);
	if (row == NULL) {
		snprintf(starchy_err_msg, sizeof(starchy_err_msg), "%s: malloc row failed", __func__);
		return p->n_ais;
	}

	//------- subract median background (bias level & sky brightness)
	for (i = 0; i < x; i++) row[i] = *(b + (y / 2) * x + i);
	qsort(row, x, sizeof(row[0]), rowcomp);

	m = row[x / 2];
	for (j = 0; j < y; j++)
		for (i = 0; i < x; i++)
			*(b + j * x + i) = (*(b + j * x + i) > m) ? *(b + j * x + i) - m : 0;

	//----------------- find stars loop
	m = 16000;
find_loop:
	for (j = 8; j < y - 8; j++) {
		for (i = 8; i < x - 8; i++) {
			if (*(b + j * x + i) < m + m) continue;

			k = 0;
			if (*(b + (j + 1) * x + (i + 0)) > m) k++;
			if (*(b + (j - 1) * x + (i + 0)) > m) k++;
			if (*(b + (j + 0) * x + (i + 1)) > m) k++;
			if (*(b + (j + 0) * x + (i - 1)) > m) k++;
			if (k < 2) continue;

			k = 0;
			if (*(b + j * x + i) > *(b + (j - 1) * x + (i + 0))) k++;
			if (*(b + j * x + i) > *(b + (j + 1) * x + (i + 0))) k++;
			if (*(b + j * x + i) > *(b + (j + 0) * x + (i - 1))) k++;
			if (*(b + j * x + i) > *(b + (j + 0) * x + (i + 1))) k++;
			if (k < 3) continue;

			k = 0;
			if (*(b + j * x + i) > *(b + (j - 1) * x + (i - 1))) k++;
			if (*(b + j * x + i) > *(b + (j + 1) * x + (i - 1))) k++;
			if (*(b + j * x + i) > *(b + (j - 1) * x + (i + 1))) k++;
			if (*(b + j * x + i) > *(b + (j + 1) * x + (i + 1))) k++;
			if (k < 3) continue;

			k = 0;
			if (*(b + j * x + i) >= *(b + (j - 2) * x + (i - 2))) k++;
			if (*(b + j * x + i) >= *(b + (j + 2) * x + (i - 2))) k++;
			if (*(b + j * x + i) >= *(b + (j - 2) * x + (i + 2))) k++;
			if (*(b + j * x + i) >= *(b + (j + 2) * x + (i + 2))) k++;
			if (k < 4) continue;

			k = *(b + (j - 3) * x + i) + *(b + (j + 3) * x + i) +
				 *(b + j * x + (i - 3)) + *(b + j * x + (i + 3));
			if (*(b + j * x + i) <= k / 3) continue;

			k = *(b + (j - 4) * x + (i - 4)) + *(b + (j + 4) * x + (i + 4)) +
				 *(b + (j + 4) * x + (i - 4)) + *(b + (j - 4) * x + (i + 4));
			if (*(b + j * x + i) <= 2 * k / 3) continue;

			n = 0;
			for (k = -4; k < +4; k++) {
				ai[n++] = *(b + (j + 4) * x + (i + k));
				ai[n++] = *(b + (j - k) * x + (i + 4));
				ai[n++] = *(b + (j - 4) * x + (i - k));
				ai[n++] = *(b + (j + k) * x + (i - 4));
			}
			qsort(ai, n, sizeof(ai[0]), aicomp);

			n = 0;
			for (k = -1; k < +1; k++) {
				ai2[n++] = *(b + (j + 1) * x + (i + k));
				ai2[n++] = *(b + (j - k) * x + (i + 1));
				ai2[n++] = *(b + (j - 1) * x + (i - k));
				ai2[n++] = *(b + (j + k) * x + (i - 1));
			}
			qsort(ai2, n, sizeof(ai2[0]), aicomp);

			if (ai2[0] <= ai[0]) continue;
			if (ai2[3] <= ai[9]) continue;
			if (ai2[4] <= ai[25]) continue;

			//-------------------------- find the centroid
			pixsum = 0;
			star_hw = 4;
			xcen = 0.0;
			ycen = 0.0;
			for (i0 = i - star_hw; i0 <= i + star_hw; i0++) {
				for (j0 = j - star_hw; j0 <= j + star_hw; j0++) {
					if ((i0 >= 0) && (i0 < x) && (j0 >= 0) && (j0 < y)) {
						n = *(b + j0 * x + i0);
						pixsum += n;
						xcen += i0 * n;
						ycen += j0 * n;
					}
				}
			}
			if (pixsum > 0) {
				xcen /= pixsum;
				ycen /= pixsum;
			} else {
				xcen = i;
				ycen = j;
			}
			xicen = floor(xcen + 0.5);
			yicen = floor(ycen + 0.5);

			if (!(p->n_ais < STARCHY_NMAX_AIS)) {
				goto end_find_loop;
			}

			//--------- search for and remove duplicate stars
			for (k = 0; k < p->n_ais; k++) {
				d0 = hypot(p->ais[k].xcen - xcen, p->ais[k].ycen - ycen);
				if ((pixsum > 4 * p->ais[k].pixsum) && (d0 < 8.0)) break;
			}
			if (k != p->n_ais) {
				for (l = k; l < p->n_ais - 1; l++) {
					memcpy(&p->ais[l], &p->ais[l+1], sizeof(p->ais[0]));
				}
				p->n_ais--;
			} else {
				for (k = 0; k < p->n_ais; k++) {
					d0 = hypot(p->ais[k].xcen - xcen, p->ais[k].ycen - ycen);
					if (d0 < 4.0) break;
				}
				if (k != p->n_ais) continue;
			}

			p->ais[p->n_ais].xicen = xicen;
			p->ais[p->n_ais].yicen = yicen;
			p->ais[p->n_ais].xcen = xcen;
			p->ais[p->n_ais].ycen = ycen;
			p->ais[p->n_ais].pixsum = pixsum;
			p->n_ais++;
		}
	}
	m /= 2;
	if (((p->n_ais < 20) && (m > 100)) || ((p->n_ais < 15) && (m > 50))) {
		goto find_loop;
	}
end_find_loop:
	qsort(p->ais, p->n_ais, sizeof(p->ais[0]), aiscomp);	// sort by brightness

	free(row);
	return p->n_ais;
}


/*
 * A hash table for star file pointers (binary keys) with sequential access.
 */
struct htf_struct {
	off_t key;		// fname
	struct htf_struct *next, *next_seq;
};

#define NMAX_HTF	2000

static struct htf_struct *htf_head, *htf[NMAX_HTF];
static int htf_collisions;

static int htf_hash(off_t key)
{
	unsigned hash, k;

	hash = 0;
	k = (unsigned) key;
	hash ^= (k <<  0);
	hash ^= (k >>  4);
	hash ^= (k <<  8);
	hash ^= (k >> 12);
	hash ^= (k << 16);
	return (hash % NMAX_HTF);
}

static struct htf_struct *htf_lookup(off_t key)
{
	struct htf_struct *np;

	for (np = htf[htf_hash(key)]; np != NULL; np = np->next)
		if (key == np->key) return np;
	return NULL;
}

static struct htf_struct *htf_insert(off_t key)
{
	struct htf_struct *np;
	int i;

	np = (struct htf_struct *) malloc(sizeof(*np));
	if (np == NULL) {
		snprintf(starchy_err_msg, sizeof(starchy_err_msg), "%s: malloc failed", __func__);
		return NULL;
	}
	np->key = key;
	i = htf_hash(key);
	if (htf[i] != NULL) htf_collisions++;
	np->next = htf[i];
	htf[i] = np;
	np->next_seq = htf_head;
	htf_head = np;
	return np;
}


/*-------------------------------------- get_members
 * This is to get stars in proximity to a particular star, pointed to by 'sp',
 * by recursively searching through the member list, and adding associated
 * stars to a hash table. The 'depth' parameter limits the recursive depth.
 */
static void get_members(struct starchy_struct *p, off_t sp, int depth)
{
	int n;
	struct starchy_star_struct st0;
	struct starchy_member_struct mb0;
	off_t off0;

	if (htf_lookup(sp) != NULL) return;	// already exists?
	if (htf_insert(sp) == NULL) return;	// add sp to hash table

	if (depth <= 0) return;

	lseek(p->sf->fd_star, sp, SEEK_SET);
	n = read(p->sf->fd_star, &st0, sizeof(st0));
	if (n != sizeof(st0)) {
		snprintf(starchy_err_msg, sizeof(starchy_err_msg), "%s: read star file failed", __func__);
		return;
	}

	for (off0 = st0.member_head; off0 != 0; off0 = mb0.next) {
		lseek(p->sf->fd_member, off0, SEEK_SET);
		n = read(p->sf->fd_member, &mb0, sizeof(mb0));
		if (n != sizeof(mb0)) {
			snprintf(starchy_err_msg, sizeof(starchy_err_msg), "%s: read member file failed", __func__);
			break;
		}
		get_members(p, mb0.bsp, depth - 1);
		get_members(p, mb0.csp, depth - 1);
	}
	return;
}

/*-------------------------------------- starchy_try
 * Count the number of coincident stars for the provided triplet.
 * This function fills the ats (catalog stars) array.
 * Parameters i0 and j0 are the indexes for stars A and B in the array ais.
 * Parameters i1 and j1 are the indexes for stars A and B in the array ats.
 * Returns true for success.
 */
bool starchy_try(struct starchy_struct *p, float fov_x, float fov_y,
		 int i0, int j0, struct starchy_t_struct *t, int x, int y,
		 uint16_t *b, struct starchy_tmatch_struct *tm,
		 int *i1, int *j1)
{
	bool ret;
	int i, j, k, n;
	int k1;
	struct htf_struct *htf0;
	struct starchy_star_struct st0, st2;
	struct starchy_utree_hdr_struct uh;
	double rad, decd;
	double d0, d1, d2, d3, d4;
	struct starchy_xy_struct xy0, xy1, xy2, xy3;	// V0, V1, U0, U1
	int n_cand, n_coin;
	bool bool0;

	ret = false;
	p->n_ats = 0;

	//--- star A is the reference direction for rsep_rorn
	lseek(p->sf->fd_star, t->asp, SEEK_SET);
	read(p->sf->fd_star, &st0, sizeof(st0));

	rad = st0.rad;
	decd = st0.decd;

	/*
	 * There are two different ways to fill the array ats with stars
	 * in proximity to the triplet. The first way is to recursively
	 * follow the triplet member linked lists, adding them to a hash
	 * table first, and then to array ats. The second way is to search the
	 * kd-tree 'u'.
	 */
	if (relaxed_fit) goto fill_ats_using_u;

	//--- initialize the hash table
	htf_head = NULL;
	for (i = 0; i < NMAX_HTF; i++) htf[i] = NULL;
	htf_collisions = 0;

	//--- fill the hash table
	get_members(p, t->asp, 3);
	get_members(p, t->bsp, 3);
	get_members(p, t->csp, 3);

	//--- fill the array ats
	d4 = hypot(fov_x * DEG2RAD, fov_y * DEG2RAD);
	n = 0;
	*i1 = -1;
	*j1 = -1;
	bool0 = true;		// first add A (i1) and B (j1) only
scan_hash_table:
	for (htf0 = htf_head; htf0 != NULL; htf0 = htf0->next_seq) {
		if (!(n < STARCHY_NMAX_ATS)) {
			snprintf(starchy_err_msg, sizeof(starchy_err_msg), "%s: reached catalog star limit (use a smaller search range value)", __func__);
			break;
		}
		if (bool0) {
			if ((*i1 >= 0) && (*j1 >= 0))
				break;
			else if (htf0->key == t->asp)
				*i1 = n;
			else if (htf0->key == t->bsp)
				*j1 = n;
			else
				continue;
		} else {
			if ((htf0->key == t->asp) || (htf0->key == t->bsp)) continue;
		}
		lseek(p->sf->fd_star, htf0->key, SEEK_SET);
		read(p->sf->fd_star, &st0, sizeof(st0));

		memset(&p->ats[n], 0, sizeof(p->ats[0]));

		p->ats[n].catalog = st0.catalog;
		strncpy(p->ats[n].catalog_id, st0.catalog_star_id, sizeof(p->ats[0].catalog_id) - 1);

		rsep_rorn(rad*DEG2RAD, decd*DEG2RAD,
			st0.rad*DEG2RAD, st0.decd*DEG2RAD, &d0, &d1, &d2, &d3);
		if (d0 > 0.75 * d4) continue;

		p->ats[n].rad = st0.rad;
		p->ats[n].decd = st0.decd;
		p->ats[n].r = d0;
		p->ats[n].a = d1;
		p->ats[n].x = d2;
		p->ats[n].y = d3;
		p->ats[n].Vt = st0.v_mag;
		p->ats[n].Bt = st0.b_mag;

		n++;
	}
	if (bool0) {
		bool0 = false;
		goto scan_hash_table;
	}
	p->n_ats = n;

	//--- erase the hash table
	while (htf_head != NULL) {
		htf0 = htf_head;
		htf_head = htf_head->next_seq;
		free(htf0);
	}
	goto end_fill_ats;

fill_ats_using_u:
	p->nmax_ua = 10 * STARCHY_NMAX_ATS;
	p->ua = malloc(p->nmax_ua * sizeof(p->ua[0]));
	if (p->ua == NULL) {
		snprintf(starchy_err_msg, sizeof(starchy_err_msg), "%s: malloc failed", __func__);
		goto out;
	}
	p->n_ua = 0;

	d0 = sin(hypot(fov_x * DEG2RAD, fov_y * DEG2RAD));

	p->u_min.u[0] = st0.u[0] - d0;
	p->u_min.u[1] = st0.u[1] - d0;
	p->u_min.u[2] = st0.u[2] - d0;

	p->u_max.u[0] = st0.u[0] + d0;
	p->u_max.u[1] = st0.u[1] + d0;
	p->u_max.u[2] = st0.u[2] + d0;

	lseek(p->sf->fd_utree, 0, SEEK_SET);
	read(p->sf->fd_utree, &uh, sizeof(uh));

	starchy_kdtree_u_search(p, uh.root, 0);

	/*
	 * Eliminate stars that are too far away
	 */
	for (j = p->n_ua - 1; j >= 0; j--) {
		d1 = acos(inner3_f(st0.u, p->ua[j].u));
		if (d1 < d0) continue;

		// eliminate it
		for (k = j; k < p->n_ua - 1; k++)
			memcpy(&p->ua[k], &p->ua[k + 1], sizeof(p->ua[0]));
		p->n_ua--;
	}
//	printf("%d stars selected\n", p->n_ua);

	//--- fill the array ats
	n = 0;

	//--- add star A
	lseek(p->sf->fd_star, t->asp, SEEK_SET);
	read(p->sf->fd_star, &st2, sizeof(st2));

	memset(&p->ats[n], 0, sizeof(p->ats[0]));

	p->ats[n].catalog = st2.catalog;
	strncpy(p->ats[n].catalog_id, st2.catalog_star_id, sizeof(p->ats[0].catalog_id) - 1);

	rsep_rorn(rad * DEG2RAD, decd * DEG2RAD,
		st2.rad * DEG2RAD, st2.decd * DEG2RAD, &d0, &d1, &d2, &d3);

	p->ats[n].rad = st2.rad;
	p->ats[n].decd = st2.decd;
	p->ats[n].r = d0;
	p->ats[n].a = d1;
	p->ats[n].x = d2;
	p->ats[n].y = d3;
	p->ats[n].Vt = st2.v_mag;
	p->ats[n].Bt = st2.b_mag;
	*i1 = n;
	n++;

	//--- add star B
	lseek(p->sf->fd_star, t->bsp, SEEK_SET);
	read(p->sf->fd_star, &st2, sizeof(st2));

	memset(&p->ats[n], 0, sizeof(p->ats[0]));

	p->ats[n].catalog = st2.catalog;
	strncpy(p->ats[n].catalog_id, st2.catalog_star_id, sizeof(p->ats[0].catalog_id) - 1);

	rsep_rorn(rad * DEG2RAD, decd * DEG2RAD,
		st2.rad * DEG2RAD, st2.decd * DEG2RAD, &d0, &d1, &d2, &d3);

	p->ats[n].rad = st2.rad;
	p->ats[n].decd = st2.decd;
	p->ats[n].r = d0;
	p->ats[n].a = d1;
	p->ats[n].x = d2;
	p->ats[n].y = d3;
	p->ats[n].Vt = st2.v_mag;
	p->ats[n].Bt = st2.b_mag;
	*j1 = n;
	n++;

	//--- add other stars
	for (i = 0; i < p->n_ua; i++) {
		if (!(n < STARCHY_NMAX_ATS)) {
			snprintf(starchy_err_msg, sizeof(starchy_err_msg), "%s: reached catalog star limit (use a smaller search range value)", __func__);
			break;
		}
		if ((p->ua[i].sp == t->asp) || (p->ua[i].sp == t->bsp))
			continue;

		lseek(p->sf->fd_star, p->ua[i].sp, SEEK_SET);
		read(p->sf->fd_star, &st2, sizeof(st2));

		memset(&p->ats[n], 0, sizeof(p->ats[0]));

		p->ats[n].catalog = st2.catalog;
		strncpy(p->ats[n].catalog_id, st2.catalog_star_id, sizeof(p->ats[0].catalog_id) - 1);

		rsep_rorn(rad * DEG2RAD, decd * DEG2RAD,
			st2.rad * DEG2RAD, st2.decd * DEG2RAD, &d0, &d1, &d2, &d3);

		p->ats[n].rad = st2.rad;
		p->ats[n].decd = st2.decd;
		p->ats[n].r = d0;
		p->ats[n].a = d1;
		p->ats[n].x = d2;
		p->ats[n].y = d3;
		p->ats[n].Vt = st2.v_mag;
		p->ats[n].Bt = st2.b_mag;
		n++;
	}
	p->n_ats = n;

	free(p->ua);
end_fill_ats:

	//--- compare image stars with catalog stars
	if ((*i1 < 0) || (*j1 < 0)) {
		snprintf(starchy_err_msg, sizeof(starchy_err_msg), "%s: A or B not in array ats", __func__);
		goto out;
	}

	xy0.x = p->ais[j0].x - p->ais[i0].x;	// dV of two image stars
	xy0.y = p->ais[j0].y - p->ais[i0].y;
	xy0.r = hypot(xy0.x, xy0.y);
	xy0.a = atan2(xy0.y, xy0.x);

	xy2.x = p->ats[*j1].x - p->ats[*i1].x;	// dV of two catalog stars
	xy2.y = p->ats[*j1].y - p->ats[*i1].y;
	xy2.r = hypot(xy2.x, xy2.y);
	xy2.a = atan2(xy2.y, xy2.x);

	n_cand = 0;
	n_coin = 0;
	for (k = 0; k < p->n_ats; k++) {	// count stars that coincide
		if ((k == *i1) || (k == *j1)) continue;

		xy3.x = p->ats[k].x - p->ats[*i1].x;
		xy3.y = p->ats[k].y - p->ats[*i1].y;
		xy3.r = hypot(xy3.x, xy3.y);
		xy3.a = atan2(xy3.y, xy3.x);

		xy1.x = xy3.r * (xy0.r/xy2.r) * cos(xy3.a + (xy0.a - xy2.a)) + p->ais[i0].x;
		xy1.y = xy3.r * (xy0.r/xy2.r) * sin(xy3.a + (xy0.a - xy2.a)) + p->ais[i0].y;

		i = (xy1.x * RAD2DEG / fov_x) * x + 0.5;
		j = (xy1.y * RAD2DEG / fov_y) * y + 0.5;
		if ((i < 4) || (i >= x - 4)) continue;
		if ((j < 4) || (j >= y - 4)) continue;

		n_cand++;
		if (*(b + j * x + i) != 0) {
			k1 = *(b + j * x + i) - 1;
			if ((k1 != i0) && (k1 != j0)) n_coin++;
		}
//		if ((n_cand > 20) && (n_coin == 0)) break;
	}

	if (((n_cand >=  3) && (n_coin >= (4 * n_cand) / 5)) ||
	    ((n_cand >=  8) && (n_coin >= (3 * n_cand) / 5)) ||
	    ((n_cand >= 15) && (n_coin >= (2 * n_cand) / 5)) ||
	    ((n_cand >= 30) && (n_coin >= (1 * n_cand) / 5)) ||
	    ((n_cand >= 48) && (n_coin >= (1 * n_cand) / 6)))
		ret = true;
	else if ((relaxed_fit == true) && (n_coin > (p->n_ais/10 + 1)))
		ret = true;

	memcpy(&tm->t, t, sizeof(tm->t));
	tm->i0 = i0;
	tm->j0 = j0;
	tm->i1 = *i1;
	tm->j1 = *j1;
	tm->n_coin = n_coin;
	tm->scale = (xy0.r / xy2.r);
	tm->rot = (xy0.a - xy2.a);
	tm->n_ats = p->n_ats;
	for (i = 0; i < p->n_ats; i++)
		memcpy(&tm->ats[i], &p->ats[i], sizeof(tm->ats[0]));
#if 0
	if (ret)
		printf("%3d n_cand, %3d n_coin, %6.3f scale, %6.2f rot, %6.1lf %6.1lf\n",
			n_cand, n_coin, tm->scale, tm->rot, rad, decd);
#endif
out:
	return ret;
}


//--------------------------------------------------- tmatchcomp_coin
static int tmatchcomp_coin(const void *a, const void *b)
{
   struct starchy_tmatch_struct *p, *q;

   p = (struct starchy_tmatch_struct *) a;
   q = (struct starchy_tmatch_struct *) b;
   return (q->n_coin - p->n_coin);
}



/*-------------------------------------- starchy_match
 * The way the image buffer parameter b is passed to this function is probably
 * not standard C.
 */
bool starchy_match(struct starchy_struct *p, float fov_x, float fov_y,
			float fov_err, int x, int y, uint16_t *b0)
{
	uint16_t *b;
	bool ret;
        struct starchy_ttree_hdr_struct th;
	int i, j, k, l, n0;
	float w0[2], w1[2], w2[2], w3[2];	// 2-D vectors
	double d0, d1, d2, d3, d4, d5;
	bool bool0;
	double scale, rot, rad_center, decd_center;
	double scale_avg, rot_avg;
	int i0, j0, i1, j1, k1;
	struct starchy_xy_struct xy0, xy1, xy2, xy3;	// V0, V1, U0, U1

	ret = false;

	//--- allocate unmodified image buffer storage
	b = malloc(sizeof(b[0]) * x * y);
	if (b == NULL) {
		snprintf(starchy_err_msg, sizeof(starchy_err_msg), "%s: malloc b failed", __func__);
		return ret;
	}
	memcpy(b, b0, sizeof(b[0]) * x * y);

	p->sf = starchy_fopen(false);
	if (p->sf == NULL) {
		goto out_0;
	}
	lseek(p->sf->fd_ttree, 0, SEEK_SET);
	read(p->sf->fd_ttree, &th, sizeof(th));

	p->nmax_ta = 1000;
	p->ta = malloc(p->nmax_ta * sizeof(struct starchy_t_struct));
	if (p->ta == NULL) {
		snprintf(starchy_err_msg, sizeof(starchy_err_msg), "%s: malloc ta failed", __func__);
		goto out_1;
	}

	//--- find image stars
	starchy_find_image_stars(p, x, y, b);		// b is clobbered
	if (p->n_ais < 5) {
		snprintf(starchy_err_msg, sizeof(starchy_err_msg), "%s: too few image stars", __func__);
		goto out_2;
	}
	for (i = 0; i < p->n_ais; i++) {
		p->ais[i].x = fov_x * DEG2RAD * (p->ais[i].xcen/x);
		p->ais[i].y = fov_y * DEG2RAD * (p->ais[i].ycen/y);
		p->ais[i].r = hypot(p->ais[i].x, p->ais[i].y);
		p->ais[i].a = atan2(p->ais[i].y, p->ais[i].x);
	}

	//--- try matching just the brightest image stars first
	n_tmatch = 0;
	n0 = (p->n_ais > 10) ? 10 : p->n_ais;
	relaxed_fit = false;

make_blobs:
	memset(b, 0, x * y * sizeof(b[0]));
	for (k = 0; k < n0; k++) {
		i = p->ais[k].xicen;
		j = p->ais[k].yicen;

		if ((i < 3) || (i >= x - 3)) continue;
		if ((j < 3) || (j >= y - 3)) continue;

		*(b + (j + 1) * x + (i - 1)) = k + 1;
		*(b + (j + 1) * x + (i + 0)) = k + 1;
		*(b + (j + 1) * x + (i + 1)) = k + 1;
		*(b + (j + 0) * x + (i - 1)) = k + 1;
		*(b + (j + 0) * x + (i + 0)) = k + 1;
		*(b + (j + 0) * x + (i + 1)) = k + 1;
		*(b + (j - 1) * x + (i - 1)) = k + 1;
		*(b + (j - 1) * x + (i + 0)) = k + 1;
		*(b + (j - 1) * x + (i + 1)) = k + 1;

		if (relaxed_fit == false) continue;

		/*
		 * A bigger star-pixel-target-region allows for less precise
		 * matching, but probably also causes more incorrect matches.
		 */
		*(b + (j + 2) * x + (i - 1)) = k + 1;
		*(b + (j + 2) * x + (i + 0)) = k + 1;
		*(b + (j + 2) * x + (i + 1)) = k + 1;
		*(b + (j + 1) * x + (i - 2)) = k + 1;
		*(b + (j + 1) * x + (i + 0)) = k + 1;
		*(b + (j + 1) * x + (i + 2)) = k + 1;
		*(b + (j + 0) * x + (i - 2)) = k + 1;
		*(b + (j + 0) * x + (i + 0)) = k + 1;
		*(b + (j + 0) * x + (i + 2)) = k + 1;
		*(b + (j - 1) * x + (i - 2)) = k + 1;
		*(b + (j - 1) * x + (i + 0)) = k + 1;
		*(b + (j - 1) * x + (i + 2)) = k + 1;
		*(b + (j - 2) * x + (i - 1)) = k + 1;
		*(b + (j - 2) * x + (i + 0)) = k + 1;
		*(b + (j - 2) * x + (i + 1)) = k + 1;
	}


	/*
	 * For combinations of three image stars, search for matching
	 * triplets in the index files. Triplets of image stars must be
	 * searched for in every combination.
	 */

	for (i = 0; i < n0; i++) {			// star A index
		for (j = 0; j < n0; j++) {		// star B index
			if (j == i) continue;

			d4 = (p->ais[j].xcen - p->ais[i].xcen);
			d5 = (p->ais[j].ycen - p->ais[i].ycen);

			// too close together?
			if (hypot(d4, d5) < 0.1 * hypot(x, y)) continue;

			// vector AB
			w0[0] = d4/x * fov_x * DEG2RAD;
			w0[1] = d5/y * fov_y * DEG2RAD;

			d0 = mag2_f(w0);		// angular separation

			p->t_min.ab = (1 - fov_err/100) * d0;
			p->t_max.ab = (1 + fov_err/100) * d0;

			w0[0] /= d0;
			w0[1] /= d0;

			for (k = 0; k < n0; k++) {	// star C index
				if ((k == i) || (k == j)) continue;

				// vector AC
				w1[0] = (p->ais[k].xcen - p->ais[i].xcen)/x * fov_x * DEG2RAD;
				w1[1] = (p->ais[k].ycen - p->ais[i].ycen)/y * fov_y * DEG2RAD;
				d1 = mag2_f(w1);

				w1[0] /= d1;
				w1[1] /= d1;

				d3 = acos(inner2_f(w0, w1));	// vertex A angle

				p->t_min.a = 0.995 * d3;
				p->t_max.a = 1.005 * d3;

				if (p->t_max.a < th.min_interior_angle) continue;

				// vector BC
				w2[0] = (p->ais[k].xcen - p->ais[j].xcen)/x * fov_x * DEG2RAD;
				w2[1] = (p->ais[k].ycen - p->ais[j].ycen)/y * fov_y * DEG2RAD;
				d2 = mag2_f(w2);

				w2[0] /= d2;
				w2[1] /= d2;

				w3[0] = -w0[0];
				w3[1] = -w0[1];

				d3 = acos(inner2_f(w2, w3));	// vertex B angle

				p->t_min.b = 0.995 * d3;
				p->t_max.b = 1.005 * d3;

				//---- kd-tree search
				p->n_ta = 0;
				kdtree_t_search(p, th.root, 0);
				if (p->n_ta == 0) continue;

//				printf("%d triplets selected\n", p->n_ta);

				for (l = 0; l < p->n_ta; l++) {
					if (n_tmatch >= NMAX_TMATCH) break;

					bool0 = starchy_try(p, fov_x, fov_y, i, j, &p->ta[l], x, y, b, &tmatch[n_tmatch], &i1, &j1);
					if (bool0) n_tmatch++;
				}
				if (n_tmatch >= NMAX_TMATCH) goto match_continue;
			}
		}
	}
	if (n_tmatch > 0) goto match_continue;

	/*
	 * match failed, try again?
	 */
	if (relaxed_fit == true) goto out_2;

	relaxed_fit = true;
	n0 = p->n_ais;
	goto make_blobs;

match_continue:
//	printf("%d matches\n", n_tmatch);

	//---- select the match with the most coincident stars
	qsort(tmatch, n_tmatch, sizeof(tmatch[0]), tmatchcomp_coin);

	//--- copy the catalog stars from the best match
	p->n_ats = tmatch[0].n_ats;
	for (i = 0; i < p->n_ats; i++)
		memcpy(&p->ats[i], &tmatch[0].ats[i], sizeof(p->ats[0]));

	scale = tmatch[0].scale;
	rot = tmatch[0].rot;
	p->n_coin = tmatch[0].n_coin;

	i0 = tmatch[0].i0;
	j0 = tmatch[0].j0;
	i = tmatch[0].i1;
	j = tmatch[0].j1;

	//--- store a list of matching stars positions in array ams[].
	scale_avg = 0.0;
	rot_avg = 0.0;

	xy0.x = p->ais[j0].x - p->ais[i0].x;	// dV of two image stars
	xy0.y = p->ais[j0].y - p->ais[i0].y;
	xy0.r = hypot(xy0.x, xy0.y);
	xy0.a = atan2(xy0.y, xy0.x);

	xy2.x = p->ats[j].x - p->ats[i].x;	// dV of two Tycho stars
	xy2.y = p->ats[j].y - p->ats[i].y;
	xy2.r = hypot(xy2.x, xy2.y);
	xy2.a = atan2(xy2.y, xy2.x);

	p->n_ams = 0;
	p->n_ans = 0;
	for (k = 0; k < p->n_ats; k++) {	// catalog star index
		if (k == i) continue;

		xy3.x = p->ats[k].x - p->ats[i].x;
		xy3.y = p->ats[k].y - p->ats[i].y;
		xy3.r = hypot(xy3.x, xy3.y);
		xy3.a = atan2(xy3.y, xy3.x);

		xy1.x = xy3.r * (xy0.r / xy2.r) * cos(xy3.a + (xy0.a - xy2.a)) + p->ais[i0].x;
		xy1.y = xy3.r * (xy0.r / xy2.r) * sin(xy3.a + (xy0.a - xy2.a)) + p->ais[i0].y;

		i1 = (xy1.x * RAD2DEG / fov_x) * x + 0.5;
		j1 = (xy1.y * RAD2DEG / fov_y) * y + 0.5;
		if (!((i1 >= 0) && (i1 < x) && (j1 >= 0) && (j1 < y))) continue;

		if (*(b + j1 * x + i1) == 0) goto tm_ans;

		k1 = *(b + j1 * x + i1) - 1;		// image star index
		if (k1 == i0) continue;

		if (!(p->n_ams < STARCHY_NMAX_AMS)) continue;

		p->ams[p->n_ams].pixsum = p->ais[k1].pixsum;
		p->ams[p->n_ams].xicen = p->ais[k1].xicen;
		p->ams[p->n_ams].yicen = p->ais[k1].yicen;
		p->ams[p->n_ams].xcen = p->ais[k1].xcen;
		p->ams[p->n_ams].ycen = p->ais[k1].ycen;

		p->ams[p->n_ams].rad  = p->ats[k].rad;
		p->ams[p->n_ams].decd = p->ats[k].decd;
		strcpy(p->ams[p->n_ams].catalog_id, p->ats[k].catalog_id);

		xy1.x = p->ais[k1].x - p->ais[i0].x;
		xy1.y = p->ais[k1].y - p->ais[i0].y;
		xy1.r = hypot(xy1.x, xy1.y);
		xy1.a = atan2(xy1.y, xy1.x);

		//--- too close together?
		if (xy1.r < 0.2 * fov_x * DEG2RAD) continue;

		if (isnan(xy1.r / xy3.r) || isinf(xy1.r / xy3.r)) continue;

		scale_avg += (xy1.r / xy3.r);
		d0 = (xy1.a - xy3.a);
		if (p->n_ams == 0) {
			d1 = d0;
		} else {
			while (d0 > d1 + PI) d0 -= 2 * PI;
			while (d0 < d1 - PI) d0 += 2 * PI;
		}
		rot_avg += d0;

		p->n_ams++;
		continue;

tm_ans:
		//----- array of non-matching catalog stars
		if (!(p->n_ans < STARCHY_NMAX_ANS)) continue;

		p->ans[p->n_ans].xcen = (xy1.x * RAD2DEG / fov_x) * x;
		p->ans[p->n_ans].ycen = (xy1.y * RAD2DEG / fov_y) * y;
		p->ans[p->n_ans].rad  = p->ats[k].rad;
		p->ans[p->n_ans].decd = p->ats[k].decd;
		strcpy(p->ans[p->n_ans].catalog_id, p->ats[k].catalog_id);
		p->n_ans++;
	}

	if (p->n_ams > 0) {
		scale_avg /= p->n_ams;
		rot_avg /= p->n_ams;
		scale = scale_avg;
		rot = rot_avg;
	}

	for (i = 0; i < p->n_ams; i++) {
		p->ams[i].u[0] = cos(p->ams[i].decd * DEG2RAD) * cos(p->ams[i].rad * DEG2RAD);
		p->ams[i].u[1] = cos(p->ams[i].decd * DEG2RAD) * sin(p->ams[i].rad * DEG2RAD);
		p->ams[i].u[2] = sin(p->ams[i].decd * DEG2RAD);
	}


	//---------------------------- calculate catalog match results
	p->scale = scale;	// actual fov = (estimated fov)/scale
	p->rot = rot;		// radians

	p->xdim = x;
	p->ydim = y;

	p->fov_x = fov_x/scale;
	p->fov_y = fov_y/scale;

	//------------ FITS WCS
	strcpy(p->radesys, "'ICRS    '");
	strcpy(p->ctype1, "'RA---TAN'");
	p->crval1 = p->ats[tmatch[0].i1].rad;		// degrees
	p->cdelt1 = (fov_x / scale) / x;		// degrees/pixel
	p->crpix1 = p->ais[tmatch[0].i0].xcen;
	p->pc1_1  =  cos(rot);
	p->pc1_2  =  sin(rot);
	p->pc2_1  = -sin(rot);
	p->pc2_2  =  cos(rot);
	strcpy(p->ctype2, "'DEC--TAN'");
	p->crval2 = p->ats[tmatch[0].i1].decd;		// degrees
	p->cdelt2 = (fov_y / scale) / y;		// degrees/pixel
	p->crpix2 = p->ais[tmatch[0].i0].ycen;

	starchy_image2radec(p, x/2, y/2, &rad_center, &decd_center);

	p->rad_center = rad_center;
	p->decd_center = decd_center;

	ret = true;

out_2:
	free(p->ta);
out_1:
	p->n_ta = 0;
	p->n_ua = 0;
	p->ta = NULL;
	p->ua = NULL;

	starchy_fclose(p->sf);
out_0:
	free(b);
	return ret;
}


/*----------------------------------------- starchy_radec2image
 * Convert (RA, Dec) to (x, y) pixel coordinates using the
 * FITS WCS values determined by the matching function.
 */
void starchy_radec2image(p, rad, decd, fx, fy)
struct starchy_struct *p;
double rad, decd;
float *fx, *fy;
{
	double d0, d1, d2, d3, d4, d5, d6, d7;

	rsep_rorn(p->crval1 * DEG2RAD, p->crval2 * DEG2RAD,
		rad * DEG2RAD, decd * DEG2RAD, &d0, &d1, &d2, &d3);

	d2 *= RAD2DEG / p->cdelt1;
	d3 *= RAD2DEG / p->cdelt2;

	d4 = d2 * p->pc1_1 + d3 * p->pc2_1;
	d5 = d2 * p->pc1_2 + d3 * p->pc2_2;

	d6 = d4 + p->crpix1;
	d7 = d5 + p->crpix2;

	*fx = d6;
	*fy = d7;
}



/*----------------------------------------- starchy_image2radec
 * Convert (x, y) pixel coordinates to (RA, Dec) using the
 * FITS WCS values determined by the matching function.
 *
 */
int starchy_image2radec(p, fx, fy, rad, decd)
struct starchy_struct *p;
float fx, fy;
double *rad, *decd;
{
	double d0, d1, d4, d5, d6, d7, d8, d9;
	double u0[3], u1[3], u2[3], u3[3], u4[3];


	//-------------
	if (p->crval2 >= 90.0) {
		snprintf(starchy_err_msg, sizeof(starchy_err_msg), "%s: decd = 90 not handled", __func__);
		return -1;
	}

	d0 = p->crval1 * DEG2RAD;
	d1 = p->crval2 * DEG2RAD;

	u0[0] = cos(d1) * cos(d0);
	u0[1] = cos(d1) * sin(d0);
	u0[2] = sin(d1);

	d0 = (p->crval1 - 90.0) * DEG2RAD;
	u1[0] = cos(d0);
	u1[1] = sin(d0);
	u1[2] = 0.0;

	cross3(u2, u1, u0);

	u1[0] = -u1[0];
	u1[1] = -u1[1];
	u1[2] = -u1[2];

#if 0
	//------- same as below
	d4 = p->pc1_1 * p->pc2_2 - p->pc1_2 * p->pc2_1;
	d6 =  p->pc2_2 / d4;
	d7 = -p->pc2_1 / d4;
	d8 = -p->pc1_2 / d4;
	d9 =  p->pc1_1 / d4;
#else
	d6 =  p->pc1_1;
	d7 =  p->pc1_2;
	d8 =  p->pc2_1;
	d9 =  p->pc2_2;
#endif

	d4 = fx - p->crpix1;
	d5 = fy - p->crpix2;

	d0 = d4 * d6 + d5 * d7;
	d1 = d4 * d8 + d5 * d9;

	d0 *= p->cdelt1;	// degrees
	d1 *= p->cdelt2;	// degrees

	d0 *= DEG2RAD;
	d1 *= DEG2RAD;

	u3[0] = d0*u1[0] + d1*u2[0];
	u3[1] = d0*u1[1] + d1*u2[1];
	u3[2] = d0*u1[2] + d1*u2[2];

	u4[0] = u0[0] + u3[0];
	u4[1] = u0[1] + u3[1];
	u4[2] = u0[2] + u3[2];

	d4 = mag3(u4);
	u4[0] /= d4;
	u4[1] /= d4;
	u4[2] /= d4;

	if ((u4[0] != 0.0) || (u4[1] != 0.0)) {
		*rad = atan2(u4[1], u4[0]) * RAD2DEG;
	} else {
		snprintf(starchy_err_msg, sizeof(starchy_err_msg), "%s: calculation error", __func__);
		return -1;
	}
	if (*rad < 0.0) *rad += 360.0;

	if ((u4[2] <= 1.0) && (u4[2] >= -1.0)) {
		*decd = asin(u4[2]) * RAD2DEG;
	} else {
		snprintf(starchy_err_msg, sizeof(starchy_err_msg), "%s: calculation error", __func__);
		return -1;
	}

	return 0;
}
