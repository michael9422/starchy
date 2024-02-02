/*
 * starchy_search.c
 *
 * This is a test program for the starchy library.
 * Compile it with:
 *
 *	# gcc starchy_search.c -Wall -lstarchy -lm -o starchy_search
 *
 * Run this program using:
 *
 *	# ./strachy_search <filename>
 *
 *-----------------------------------------------------------------------
 * Copyright (C) 2014 -2015 FSF
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


//---------------------------------------------------- main
int main(int argc, char *argv[], char *env[])
{
	struct starchy_struct sty;
	int i, j;
	FILE *fp;
	char hdr[180][80];
	uint8_t *b8;
	uint16_t *b16;
	bool bool0;
	double d0, d1;
	float f0, f1;
	char f_simple;
	int f_bitpix, f_naxis, f_naxis1, f_naxis2;
	bool f_end;
	float fov;		// degrees

	//------ read a sample FITS file
	if (argc > 1) {
		fp = fopen(argv[1], "r");
	} else {
		fprintf(stderr, "usage: # ./starchy_search <filename>\n");
		return -1;
	}
	if (fp == NULL) {
		fprintf(stderr, "%s: fopen FITS failed\n", __func__);
		return -1;
	}

	//----------- set FOV
	if (strcmp(argv[1], "wide50mm-008.fits") == 0)
		fov = 6.7;
	else
		fov = 7.3;

	printf("FOV = %.1f degrees\n", fov);

	/*
	 * For this example program and the included example FITS files, the
	 * header is not parsed. In general though, a FITS file header size is
	 * not a fixed length, and must be parsed for image format information.
	 */
	fread(hdr, 1, 36*80, fp);

	f_end = false;
	f_bitpix = 0;
	f_naxis = 0;
	f_naxis1 = 0;
	f_naxis2 = 0;

	for (i = 0; i < 36; i++) {
		sscanf(hdr[i], "SIMPLE  = %c", &f_simple);
		sscanf(hdr[i], "BITPIX  = %d", &f_bitpix);
		sscanf(hdr[i], "NAXIS   = %d", &f_naxis);
		sscanf(hdr[i], "NAXIS1  = %d", &f_naxis1);
		sscanf(hdr[i], "NAXIS2  = %d", &f_naxis2);
		if (strncmp(hdr[i], "END     ", 8) == 0) {
			f_end = true;
			break;
		}
	}

	if (f_end != true) {
		fprintf(stderr, "FITS header missing END record\n");
		fclose(fp);
		return -1;
	}

	if (f_bitpix != 8) {
		fprintf(stderr, "FITS header wrong BITPIX values\n");
		fclose(fp);
		return -1;
	}

	if ((f_naxis != 2) || (f_naxis1 == 0) || (f_naxis2 == 0)) {
		fprintf(stderr, "FITS header bad or missing NAXIS records\n");
		fclose(fp);
		return -1;
	}

	b8 = malloc(f_naxis1 * f_naxis2);
	b16 = malloc(2 * f_naxis1 * f_naxis2);

	fread(b8, 1, f_naxis1 * f_naxis2, fp);
	fclose(fp);

	for (i = 0; i < f_naxis1 * f_naxis2; i++)
		b16[i] = b8[i] << 8;

	/*
	 * The default path to the index data files can be overridden by
	 * reassigning global variables used by the starchy library.
	 * This could allow using different data files depending
	 * on the field of view of the pictures. The default data index
	 * files are installed in directory "/usr/local/share".
	 */
//	starchy_path_index = "/custom/data/path/";

	/*
	 * The parameters to the matching routine are
	 *
	 *	sty	- match info (starchy_struct *)
	 *	fov_x	- X-axis field of view degrees (float)
	 *	fov_y	- Y-axis field of view degrees (float)
	 *	fov_err	- field of view size uncertainty in percent (float)
	 *	x	- X-axis buffer size (int)
	 *	y	- Y-axis buffer size (int)
	 *	b	- image buffer pointer (uint16_t *)
	 *
	 * Note that a larger size uncertainty value 'fov_err' results in a
	 * slower search.
	 */
	printf("searching...");
	fflush(stdout);
	bool0 = starchy_match(&sty, fov, fov*480/640, 1.0, f_naxis1, f_naxis2, b16);
	printf("\n");

	/****
	 * If match failed, it may be necessary to reverse one axis for
	 * matching to succeed, so try again.
	 */
	if (bool0 == false) {
		printf("search failed, trying reversed parity...");
		fflush(stdout);
		for (j = 0; j < f_naxis2; j++)
			for (i = 0; i < f_naxis1; i++)
				*(b16 + j*f_naxis1 + i) = *(b8 + (f_naxis2 - 1 - j)*f_naxis1 + i) << 8;

		bool0 = starchy_match(&sty, fov, fov*480/640, 1.0, f_naxis1, f_naxis2, b16);
		printf("\n");
	}

	printf("%d image stars\n", sty.n_ais);
	if (bool0) {
		printf("search successful %s\n", starchy_err_msg);
	} else {
		printf("search failed %s\n", starchy_err_msg);
		goto out_1;
	}
	printf("%d catalog stars\n", sty.n_ats);
	printf("%d coincident stars\n", sty.n_coin);
	printf("\n");
	printf("FITS header keywords:\n");
	printf("\n");
	printf("RA      = %20.6lf / degrees\n", sty.rad_center);
	printf("DEC     = %20.6lf / degrees\n", sty.decd_center);
	printf("RADESYS = %-20s\n", sty.radesys);
	printf("CTYPE1  = %-20s\n", sty.ctype1);
	printf("CRVAL1  = %20.6lf / degrees\n", sty.crval1);
	printf("CDELT1  = %20e / deg/pix\n", sty.cdelt1);
	printf("CRPIX1  = %20.1lf\n", sty.crpix1 + 1.0);
	printf("PC1_1   = %20e\n", sty.pc1_1);
	printf("PC1_2   = %20e\n", sty.pc1_2);
	printf("PC2_1   = %20e\n", sty.pc2_1);
	printf("PC2_2   = %20e\n", sty.pc2_2);
	printf("CTYPE2  = %-20s\n", sty.ctype2);
	printf("CRVAL2  = %20.6lf / degrees\n", sty.crval2);
	printf("CDELT2  = %20e / deg/pix\n", sty.cdelt2);
	printf("CRPIX2  = %20.1lf\n", sty.crpix2 + 1.0);
	printf("\n");


	f0 = f_naxis1/2;
	f1 = f_naxis2/2;
	starchy_image2radec(&sty, f0, f1, &d0, &d1);
	printf("RA(x=%.1f, y=%.1f) \t = %9.3lf (degrees)\n", f0, f1, d0);
	printf("Dec(x=%.1f, y=%.1f)\t = %9.3lf (degrees)\n", f0, f1, d1);
	printf("\n");

	starchy_radec2image(&sty, d0, d1, &f0, &f1);
	printf("x(RA=%.3lf, Dec=%.3lf)\t = %7.2f\n", d0, d1, f0);
	printf("y(RA=%.3lf, Dec=%.3lf)\t = %7.2f\n", d0, d1, f1);
	printf("\n");

out_1:
	free(b16);
	free(b8);
	return 0;
}
