/*
 * starchy.h
 *
 * This header file is for the "starchy" astronomy library for coordinatizing
 * pictures of stars.
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
 * Tentatively, single precision floating point is used for the data in the
 * index files, instead of double precision, for smaller file sizes.
 *
 * The structures defined here for this library are somewhat disorganized.
 *
 */
#ifndef STARCHY_H
#define STARCHY_H	1

#define _GNU_SOURCE	1

#include <fcntl.h>
#include <float.h>
#include <inttypes.h>
#include <math.h>
#include <semaphore.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/file.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#ifndef PI
 #define PI		3.14159265
#endif


#define STARCHY_CATALOG_TYCHO		0
#define STARCHY_CATALOG_GSC		1
#define STARCHY_CATALOG_USNO_A2		2


/****
 * These are structures for the 'star' file. The format begins with one header
 * record, followed by the star records (unordered, sequential).
 */

struct starchy_star_hdr_struct {
	char magic[8];
	char date_s[30];	// "2014-12-04T00:00:00";
	int32_t l_rec;		// size of struct starchy_star_struct
	int32_t n_rec;		// count of star records
};

struct starchy_star_struct {
	double rad, decd;	// RA/Dec in degrees
	float u[3];		// unit direction vector
	float v_mag;
	float b_mag;
	char catalog;		// STARCHY_CATALOG_<cat> value
	char catalog_star_id[16];
	off_t member_head;	// triplet membership list head pointer
	int8_t n_member;
	int8_t n_a, n_b, n_c;
};

/****
 * These are structures used for the kd-tree named utree, a tree organized by
 * the components of the unit direction vectors of the stars. The file has one
 * header record, followed by the unode records. This is used for (crudely)
 * selecting stars in proximity to a specific star (direction). The related
 * structure starchy_u_struct is used in the routines, but not in the file.
 */

struct starchy_utree_hdr_struct {
	char magic[8];
	char date_s[30];
	int32_t l_rec;
	off_t root;
};

struct starchy_u_struct {
	float u[3];		// unit direction vector
	off_t sp;		// star file record pointer
};

struct starchy_unode_struct {
	float u[3];
	off_t sp;		// star file record pointer
	off_t l, r;		// left, right tree offset
};

/****
 * This is used for a temporary array.
 */
struct starchy_v_struct {
	float u[3];
	float v_mag;
	off_t sp;		// star file record pointer
};

/****
 * These structures are for the 'member' file. It is singly linked lists of
 * records identifying the other two verticies of triplets that stars are
 * members of. The list head pointers are in the star file. The file format
 * begins with one header record, followed member records. This is used for
 * finding stars in proximity to a specific star.
 */

struct starchy_member_hdr_struct {
	char magic[8];
	char date_s[30];
	int32_t l_rec;		// length of starchy_member_struct
};

struct starchy_member_struct {
	off_t bsp, csp;		// star file record pointers for verticies B, C
	off_t next;		// member file record pointer
};

/****
 * These structures are for a 3-way kd-tree named 'ttree', the triplet tree.
 * It is organized by three radian angles. For a triplet with stars A B C, the
 * first angle is the separation of A and B on the sky. The next one is the
 * interior angle for vertex A of the triangle ABC. The final angle
 * is for vertex B.
 * The structure starchy_t_struct is used in the routines but not in the file.
 */

struct starchy_ttree_hdr_struct {
	char magic[8];
	char date_s[30];
	int32_t l_rec;
	float min_interior_angle;
	off_t root;
};

struct starchy_t_struct {
	off_t asp, bsp, csp;
	float ab;
	float a;
	float b;
};

struct starchy_tnode_struct {
	off_t asp, bsp, csp;		// star record pointers (triplet)
	float ab;			// ab separation angle (kd-tree axis 0)
	float a;			// angle a (kd-tree axis 1)
	float b;			// angle b (kd-tree axis 2)
	off_t l, r;
};

/****
 * A typedef wrapper for access to the starchy index files.
 * The fd_lock file is used to enforce exclusive access to
 * index files in a directory.
 */

struct starchy_file_struct {
	char *fname_lock;
	char *fname_star;
	char *fname_utree;
	char *fname_member;
	char *fname_ttree;

	int fd_lock;
	int fd_star;
	int fd_utree;
	int fd_member;
	int fd_ttree;
};

typedef struct starchy_file_struct STARCHY_FILE;

/****
 * This structure is used for both image stars and catalog stars, but there
 * probably should be two different structures.
 */

struct starchy_xy_struct {
	int32_t pixsum;			// pixel sum (less background)
	int32_t xicen, yicen;		// center (pixels)
	double xcen, ycen;		// star centroid (pixels)
	double x, y;			// radians
	double a, r;			// polar (radians)
	double rad, decd;		// Tycho J2000 (degrees)
	float Vt, Bt;			// Tycho magnitudes
//	float ri, mi, mi_min, nn;	// image radius
	double u[3];
	char catalog;
	char catalog_id[16];
};


#define STARCHY_NMAX_AIS	50	// max image stars
#define STARCHY_NMAX_ATS	500	// max catalog stars
#define STARCHY_NMAX_AMS	50	// max matched stars
#define STARCHY_NMAX_ANS	50	// max not matched stars


struct starchy_struct {
	STARCHY_FILE *sf;

	int n_ua, nmax_ua;
	struct starchy_u_struct *ua;
	struct starchy_u_struct u_min, u_max;

	int n_ta, nmax_ta;
	struct starchy_t_struct *ta;
	struct starchy_t_struct t_min, t_max;

	/*
	 * If these arrays were dynamically allocated instead, then to be
	 * accessible by programs using the library, they would have to be
	 * freed by the calling program.
	 */
	struct starchy_xy_struct ais[STARCHY_NMAX_AIS];	// image stars
	struct starchy_xy_struct ats[STARCHY_NMAX_ATS];	// catalog stars
	struct starchy_xy_struct ams[STARCHY_NMAX_AMS];	// matched stars
	struct starchy_xy_struct ans[STARCHY_NMAX_ANS];	// not matched stars

	int32_t n_ais;
	int32_t n_ats;
	int32_t n_ams;
	int32_t n_ans;
	int32_t n_coin;
//	int32_t n_cand;
//	int32_t n_tmatch;
//	int32_t n_center;

	double scale, rot, rad_center, decd_center;
	int32_t xdim, ydim;
	double fov_x, fov_y;
	/*
	 * FITS file WCS values (by keyword).
	 * Note: the crpix<n> values for these library functions are for arrays
	 * beginning with an index of 0, not 1. This may be inconsistent with
	 * the FITS standard. So, it is probably necessary to add +1 to these
	 * values when writing a FITS file with these keyword fields in the
	 * header.
	 */
	char radesys[40], ctype1[40], ctype2[40];
	double crval1, cdelt1, crpix1;
	double crval2, cdelt2, crpix2;
	double pc1_1, pc1_2, pc2_1, pc2_2;

//	char err_msg[80];
};


/****
 * A structure for one potential match
 */
struct starchy_tmatch_struct {
	struct starchy_t_struct t;
	int i0, j0;
	int i1, j1;
	int n_coin;
	double scale, rot;
	int n_ats;
	struct starchy_xy_struct ats[STARCHY_NMAX_ATS];
};


/****
 * Exported (global) variables.
 * Using global variables is not thread safe, and it might not even be
 * multi-process safe (for example, for a shared library?). With these,
 * programs using the library can override the default index file names and
 * path names assigned in "starchy.c".
 * Also, error messages are returned in a global, fixed length buffer.
 */
extern char *starchy_tycho_path;
extern char *starchy_tycho_index;
extern char *starchy_tycho_catalog;

extern char *starchy_path_index;
extern char *starchy_fname_lock;
extern char *starchy_fname_tmp;
extern char *starchy_fname_star;
extern char *starchy_fname_utree;
extern char *starchy_fname_member;
extern char *starchy_fname_ttree;

extern char starchy_err_msg[];


/****
 * Library function prototypes
 */
extern STARCHY_FILE *starchy_fopen(bool);

extern void starchy_fclose(STARCHY_FILE *);

extern void starchy_kdtree_u_search(struct starchy_struct *p, off_t np,
				int depth);

extern bool starchy_match(struct starchy_struct *p,
				float fov_x, float fov_y, float fov_err,
				int x, int y, uint16_t *b);

extern void starchy_radec2image(struct starchy_struct *p,
				double rad, double decd,
				float *fx, float *fy);

extern int starchy_image2radec(struct starchy_struct *p,
				float fx, float fy,
				double *rad, double *decd);

#endif
