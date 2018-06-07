/*
 * This file is part of libmusic - frequency detection library.
 *
 * Copyright (c) 2018 Data And Signal - IT Solutions
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 *   Redistributions of source code must retain the above copyright
 *   notice, this list of conditions and the following disclaimer.
 *
 *   Redistributions in binary form must reproduce the above
 *   copyright notice, this list of conditions and the following
 *   disclaimer in the documentation and/or other materials provided
 *   with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 * COPYRIGHT HOLDERS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 * OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * Piotr Gregor <piotr@dataandsignal.com>
 * Data And Signal - IT Solutions
 *
 */


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>

#include "lm.h"


/**
 * Test of MUSIC agorithm.
 *
 * Input is DTMF "1" = { 697 Hz, 1209 Hz }. 8000 Hz, 16 samples: [10, 25].
 */

/* Reference data. */
#define SIGNALS_EXPECTED_N 2
#define SIGNALS_N (2 * (SIGNALS_EXPECTED_N))								/* Need to account for complex signals. */
#define CORR_ORD ((2 * (SIGNALS_N)) - 1)
#define XN 16

double ref_array_X[XN] = { -0.2071,   -0.7942,   -1.1108,   -0.6395,    0.5197,    1.6468,    1.9312,    1.1106,   -0.3025,   -1.3984,   -1.5514, -0.8580,    0.0096,    0.3922,    0.1753,   -0.1748 };

/* MATLAB result
 *
 * >> y2(10:25)
 * ans =
 *   Columns 1 through 11
 *      -0.2071   -0.7942   -1.1108   -0.6395    0.5197    1.6468    1.9312    1.1106   -0.3025   -1.3984   -1.5514
 *  Columns 12 through 16
 *      -0.8580    0.0096    0.3922    0.1753   -0.1748
 *
 * r_x =
 *		0.3702		0.6437		0.5489		0.1732	   -0.2132	   -0.3703	   -0.2647	   -0.0690
 *	   -0.1008		0.3702		0.6437		0.5489		0.1732	   -0.2132	   -0.3703	   -0.2647
 *	   -0.4661	   -0.1008		0.3702		0.6437		0.5489		0.1732	   -0.2132	   -0.3703
 *	   -0.5171	   -0.4661	   -0.1008		0.3702		0.6437		0.5489		0.1732	   -0.2132
 *	   -0.2860	   -0.5171	   -0.4661	   -0.1008		0.3702		0.6437		0.5489		0.1732
 *		0.0032	   -0.2860	   -0.5171	   -0.4661	   -0.1008		0.3702		0.6437		0.5489
 *		0.1307		0.0032	   -0.2860	   -0.5171	   -0.4661	   -0.1008		0.3702		0.6437
 *		0.0584		0.1307		0.0032	   -0.2860	   -0.5171	   -0.4661	   -0.1008		0.3702
 *	   -0.0583		0.0584		0.1307		0.0032	   -0.2860	   -0.5171	   -0.4661	   -0.1008
 *
 * eigenvects (U) (not needed in MUSIC) =
 *
 *		-0.0846    0.3488   -0.4748    0.5541    0.1778    0.1959    0.5104    0.0903
 *		-0.3811    0.2414   -0.5224   -0.0641   -0.4035   -0.3970   -0.2920   -0.3374
 *      -0.5024   -0.0882   -0.2875   -0.3743    0.5031    0.0787   -0.1739    0.4770
 *      -0.3285   -0.4445   -0.1048   -0.2422   -0.4236    0.5067    0.3984   -0.1731
 *       0.0505   -0.5803   -0.1801    0.0589    0.4140   -0.4871    0.2815   -0.3700
 *       0.3934   -0.3885   -0.4023    0.1241   -0.3768   -0.1516   -0.0926    0.5861
 *       0.4873    0.0101   -0.4507   -0.1995    0.2398    0.4709   -0.3288   -0.3635
 *       0.3041    0.3564   -0.1024   -0.6568   -0.0480   -0.2430    0.5190    0.0891
 *
 * eigenvals =  (squared)
 *		
 *		5.3295
 *      4.3849
 *      0.4334
 *      0.2227
 *		0.0000
 *      0.0000
 *      0.0000
 *      0.0000
 *
 * eigenvects V =
 *
 *		-0.0846    0.3488   -0.4748    0.5541    0.1778    0.1959    0.5104    0.0903
 *		-0.3811    0.2414   -0.5224   -0.0641   -0.4035   -0.3970   -0.2920   -0.3374
 *		-0.5024   -0.0882   -0.2875   -0.3743    0.5031    0.0787   -0.1739    0.4770
 *		-0.3285   -0.4445   -0.1048   -0.2422   -0.4236    0.5067    0.3984   -0.1731
 *		 0.0505   -0.5803   -0.1801    0.0589    0.4140   -0.4871    0.2815   -0.3700
 *		 0.3934   -0.3885   -0.4023    0.1241   -0.3768   -0.1516   -0.0926    0.5861
 *		 0.4873    0.0101   -0.4507   -0.1995    0.2398    0.4709   -0.3288   -0.3635
 *		 0.3041    0.3564   -0.1024   -0.6568   -0.0480   -0.2430    0.5190    0.0891
 */

void cb_dbl(double *x)
{
	fprintf(stderr, "\t%f", *x);
}

void cb_dbl_row_end(void)
{
	fprintf(stderr, "\n");
}


int main (void)
{
	int ym = XN - CORR_ORD;	/* correlation result matrix row dimension */
	int yn = CORR_ORD + 1;	/* correlation result matrix column dimension */

	uint16_t x[XN];			/* samples */

	int		M = ym;
	int		N = yn;
	double	*VT = NULL;
	int		LDVT = yn ;
	double	*WORK = NULL;
	int		LWORK = 4 * M * N * M *N + 6 * M * N + lm_max(M, N);
	
	int i = 0, j = 0;

	/* Peak detection */
	double Fs = 8000;
	double dt = 1/Fs;
	double t = 0;
	double omega = 0;
	double a_real[yn];	/* delays, real(exp(j*omega*t)) */
	double a_imag[yn];	/* delays, imag(exp(j*omega*t)) */
	double av_real = 0;	/* real result of a * V multiplication */
	double av_imag = 0;	/* imaginary result of a * V multiplication */
	uint16_t f = 0;
	double peak = 0;
	double peak_min, peak_max, peak_697, peak_1209;	/* test */

	int res;

	lm_detector_t d;


	memset(&d, 0, sizeof(d));

	d.fs = 8000;
	d.dt = 1.0 / (double) d.fs;
	d.x_n = XN;

	d.signals_n = 2;				/* we are searching for 2 frequencies */
	d.svd.JOBU = 'N';
	d.svd.JOBVT = 'A';				/* we want eigenvalues */

	res = lm_detect(&d, &ref_array_X[0], XN);
	if (res != 0) {
		printf("lm_detect failed with error %d\n", res);
		goto lmt7_fail_sys;
	}
	lm_test(res == 0, "lm_detect failed");

	lm_test(d.svd.M == M, "lm_detect failed: M");
	lm_test(d.svd.N == N, "lm_detect failed: N");
	VT = d.svd.VT;

	fprintf(stderr, "Result VT (rowwise):\n");
	lm_walk_dbl_arr_rowwise(VT, LDVT, N, cb_dbl, cb_dbl_row_end);
	fprintf(stderr, "Result VT (colwise):\n");
	lm_walk_dbl_arr_colwise(VT, LDVT, N, cb_dbl, cb_dbl_row_end);

	/* Peak detection */

	f = 0;
	peak_min = DBL_MAX;
	peak_max = DBL_MIN;
	peak_697 = peak_1209 = DBL_MIN;

	while (f < 4000) {

		i = 0;
		t = 0;
		while (i < N) {
			omega = 2 * M_PI * f;
			a_real[i] = cos(omega * t);
			a_imag[i] = sin(omega * t);
			t += dt;
			++i;
		}

		i = 0;
		av_real = av_imag = 0;
		while (i < N) {
			av_real += a_real[i] * (VT[4 + i * yn] + VT[5 + i * yn] + VT[6 + i * yn] + VT[7 + i * yn]);
			av_imag += a_imag[i] * (VT[4 + i * yn] + VT[5 + i * yn] + VT[6 + i * yn] + VT[7 + i * yn]);
			++i;
		}

		peak = 1 / (av_real * av_real + av_imag * av_imag);

		printf("peak [%04u] = %f\n", f, peak);
		printf("peak got from lm_detect [%04u] = %f\n", f, d.detections[f].peak);

		if (peak > peak_max) {
			peak_max = peak;
		}
		if (peak < peak_min) {
			peak_min = peak;
		}
		if (f == 697) {
			peak_697 = peak;
		}
		if (f == 1209) {
			peak_1209 = peak;
		}

		f++;
	}

	printf("\npeak_min = %f\n", peak_min);
	printf("peak_max = %f\n", peak_max);
	printf("peak_697 = %f\n", peak_697);
	printf("peak_1209 = %f\n", peak_1209);
	printf("fabs((peak_697 - peak_1209) / (peak_max - peak_min)) = %f, %f, %f\n", fabs(peak_697 - peak_1209), fabs(peak_max - peak_min), fabs(peak_697 - peak_1209) / (peak_max - peak_min));
	printf("fabs(peak_697 - peak_min) / (peak_max - peak_min) = %f, %f, %f\n", fabs(peak_697 - peak_min), peak_max - peak_min, fabs(peak_697 - peak_min) / (peak_max - peak_min));
	printf("fabs(peak_1209 - peak_min) / (peak_max - peak_min) = %f, %f, %f\n", fabs(peak_1209 - peak_min), peak_max - peak_min, fabs(peak_1209 - peak_min) / (peak_max - peak_min));

	printf(".info.peak_min = %f\n", d.info.peak_min);
	printf(".info.peak_max = %f\n", d.info.peak_max);
	printf(".info.peak_min_freq = %u\n", d.info.peak_min_freq);
	printf(".info.peak_max_freq = %u\n", d.info.peak_max_freq);
	printf(".info.peak_697 = %f\n", d.info.peak_697);
	printf(".info.peak_1209 = %f\n", d.info.peak_1209);
	
	lm_test(peak_697 > peak_min || peak_1209 == peak_min, "Minimum peak not less than peak 697 Hz and not less than peak 1209 Hz");
	lm_test(peak_697 == peak_max || peak_1209 == peak_max, "Maximum peak not at 697 Hz and not at 1209 Hz");
	lm_test(fabs(peak_697 - peak_1209) < 0.25 * fabs(peak_max - peak_min), "Peaks at 697 Hz and 1209 Hz differ more than 0.25 * fabs(max - min)");
	lm_test((fabs(peak_697 - peak_min) > 0.75 * fabs(peak_max - peak_min)), "Peak at 697 Hz not more than 0.75 * fabs(max - min)");
	lm_test((fabs(peak_1209 - peak_min) > 0.75 * fabs(peak_max - peak_min)), "Peak at 1209 Hz not more than 0.75 * fabs(max - min)");

	lm_deinit(&d);

	return 0;

lmt7_fail_sys:
	fprintf(stderr, "Couldn't complete the test due to system error\n");
	return -1;
}
