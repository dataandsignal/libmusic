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
 * Test Correlation.
 *
 * Input is DTMF "1" = { 697 Hz, 1209 Hz }. 8000 Hz, 16 samples: [10, 25].
 */

/* Reference data. */
#define SIGNALS_EXPECTED_N 2
#define SIGNALS_N (2 * (SIGNALS_EXPECTED_N))								/* Need to account for complex signals. */
#define CORR_ORD ((2 * (SIGNALS_N)) - 1)
#define N 16
double ref_array_X[N] = { -0.2071,   -0.7942,   -1.1108,   -0.6395,    0.5197,    1.6468,    1.9312,    1.1106,   -0.3025,   -1.3984,   -1.5514, -0.8580,    0.0096,    0.3922,    0.1753,   -0.1748 };
double ref_array_Corr[N - CORR_ORD][CORR_ORD + 1] = {
	{ 0.3702,	0.6437,		0.5489,		0.1732,		-0.2132,	-0.3703,	-0.2647,	-0.0690	},
	{ -0.1008,  0.3702,		0.6437,		0.5489,		0.1732,		-0.2132,	-0.3703,	-0.2647	},
	{ -0.4661,  -0.1008,	0.3702,		0.6437,		0.5489,		0.1732,		-0.2132,	-0.3703	},
	{ -0.5171,  -0.4661,	-0.1008,	0.3702,		0.6437,		0.5489,		0.1732,		-0.2132	},
	{ -0.2860,  -0.5171,	-0.4661,	-0.1008,    0.3702,		0.6437,		0.5489,		0.1732	},
	{ 0.0032,	-0.2860,	-0.5171,	-0.4661,	-0.1008,    0.3702,		0.6437,		0.5489	},
	{ 0.1307,	0.0032,		-0.2860,	-0.5171,	-0.4661,	-0.1008,    0.3702,		0.6437	},
	{ 0.0584,   0.1307,		0.0032,		-0.2860,	-0.5171,	-0.4661,	-0.1008,    0.3702	},
	{ -0.0583,  0.0584,		0.1307,		0.0032,		-0.2860,	-0.5171,	-0.4661,	-0.1008 }
};

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
 *		0.3702		0.6437		0.5489		0.1732		-0.2132		-0.3703		-0.2647		-0.0690
 *		-0.1008		0.3702		0.6437		0.5489		0.1732		-0.2132		-0.3703		-0.2647
 *		-0.4661		-0.1008		0.3702		0.6437		0.5489		0.1732		-0.2132		-0.3703
 *		-0.5171		-0.4661		-0.1008		0.3702		0.6437		0.5489		0.1732		-0.2132
 *		-0.2860		-0.5171		-0.4661		-0.1008		0.3702		0.6437		0.5489		0.1732
 *		0.0032		-0.2860		-0.5171		-0.4661		-0.1008		0.3702		0.6437		0.5489
 *		0.1307		0.0032		-0.2860		-0.5171		-0.4661		-0.1008		0.3702		0.6437
 *		0.0584		0.1307		0.0032		-0.2860		-0.5171		-0.4661		-0.1008		0.3702
 *		-0.0583		0.0584		0.1307		0.0032		-0.2860		-0.5171		-0.4661		-0.1008
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
	double *X;				/* samples */
	double *Y;				/* result correlation matrix of X */
	int ym = N - CORR_ORD;	/* result matrix row dimension */
	int yn = CORR_ORD + 1;	/* result matrix column dimension */

	/* helper variables */
	int i = 0, j = 0;

	fprintf(stderr, "Reference array Corr:\n");
	lm_walk_dbl_arr_rowwise(&ref_array_Corr[0][0], ym, yn, cb_dbl, cb_dbl_row_end);
	
	X = calloc(1 * N, sizeof(double));
	if (!X) {
		goto lmt5_fail_sys;
	}
	for (i = 0; i < N; ++i) {
		X[i] = ref_array_X[i];
	}

	Y = lm_corr(X, CORR_ORD, N, NULL);
	lm_test(Y != NULL, "lm_corr failed");
	
	fprintf(stderr, "Result Y (rowwise):\n");
	lm_walk_dbl_arr_rowwise(Y, ym, yn, cb_dbl, cb_dbl_row_end);
	fprintf(stderr, "Result Y (colwise):\n");
	lm_walk_dbl_arr_colwise(Y, ym, yn, cb_dbl, cb_dbl_row_end);

	free(X);
	free(Y);

	return 0;

lmt5_fail_sys:
	fprintf(stderr, "Couldn't complete the test due to system error\n");
	return -1;
}
