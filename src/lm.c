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


#include "config.h"
#include "lm.h"


double* lm_corr(double *s, uint16_t corr_ord, uint16_t n, double *inline_where)
{
	double *y = inline_where;
	uint16_t i = 0, j = 0;
	uint16_t ym = 0, yn = 0;
	double sq = 0;


	if (n < corr_ord) {
		return NULL;
	}
	
	ym = n - corr_ord;									/* for 16 sample vector: corr_ord=7, ym=9, yn=8 */
	yn = corr_ord + 1;
	sq = sqrt(n - corr_ord);

	if (y == NULL) {
		y = calloc(ym * yn, sizeof(double));
		if (!y) {
			return NULL;
		}
	}

	/* Produce delayed frames, and swap columns, rightmost is delay == 0, leftmost is delay == corr_ord.
	 * Scale them in place. */

	i = 0;
	j = 0;
	while (i < ym) {
		j = 0;
		while (j < yn) {

			y[i * yn + yn - 1 - j] = s[(i + j) % n] / sq;
			++j;
		}
		++i;
	}

	return y;
}

int lm_lapack_query(lm_detector_t *d, uint16_t x_n)
{
	/*	Arguments to LAPACK's dgeslm_.
	 *
	 *	DGESDD computes the singular value decomposition (SVD) of a real
	 *	M-by-N matrix A, optionally computing the left and right singular
	 *  vectors.  If singular vectors are desired, it uses a
	 *  divide-and-conquer algorithm.
	 *
	 *  The SVD is written
	 *
	 *          A = U * SIGMA * transpose(V)
	 *
	 *  where SIGMA is an M-by-N matrix which is zero except for its
	 *  min(m,n) diagonal elements, U is an M-by-M orthogonal matrix, and
	 *  V is an N-by-N orthogonal matrix.  The diagonal elements of SIGMA
	 *  are the singular values of A; they are real and non-negative, and
	 *  are returned in descending order.  The first min(m,n) columns of
	 *  U and V are the left and right singular vectors of A.
	 *
	 *  Note that the routine returns VT = V**T, not V.
	 *
	 *	The divide and conquer algorithm makes very mild assumptions about
	 *  floating point arithmetic. It will work on machines with a guard
	 *  digit in alm/subtract, or on those binary machines without guard
	 *  digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
	 *  Cray-2. It could conceivably fail on hexadecimal or decimal machines
	 *  without guard digits.
	 */
	char	JOBU;				/* JOBU is CHARACTER*1
								   Specifies options for computing all or part of the matrix U:
								   = 'A':  all M columns of U are returned in array U:
								   = 'S':  the first min(m,n) columns of U (the left singular
								   vectors) are returned in the array U;
								   = 'O':  the first min(m,n) columns of U (the left singular
								   vectors) are overwritten on the array A;
								   = 'N':  no columns of U (no left singular vectors) are
								   computed. */
	char	JOBVT;				/* JOBVT is CHARACTER*1
								   Specifies options for computing all or part of the matrix
								   V**T:
								   = 'A':  all N rows of V**T are returned in the array VT;
								   = 'S':  the first min(m,n) rows of V**T (the right singular
								   vectors) are returned in the array VT;
								   = 'O':  the first min(m,n) rows of V**T (the right singular
								   vectors) are overwritten on the array A;
								   = 'N':  no rows of V**T (no right singular vectors) are
								   computed. JOBVT and JOBU cannot both be 'O'. */		  
	int		M;					/* M is INTEGER. The number of rows of the input matrix A.  M >= 0. */
	int		N;					/* N is INTEGER. The number of columns of the input matrix A.  N >= 0. */
	double	*A, *A2;			/* A is DOUBLE PRECISION array, dimension (LDA,N)
								   On entry, the M-by-N matrix A. On exit, if JOBZ = 'O',
								   A is overwritten with the first N columns of U (the left singular vectors, stored columnwise) if M >= N;
								   A is overwritten with the first M rows of V**T (the right singular vectors, stored rowwise) otherwise.
								   if JOBZ .ne. 'O', the contents of A are destroyed. */
	int		LDA;				/* LDA is INTEGER. The leading dimension of the array A.  LDA >= max(1,M). */
	double	*S;					/* S is DOUBLE PRECISION array, dimension (min(M,N)). The singular values of A, sorted so that S(i) >= S(i+1). */
	double	*U;					/* U is DOUBLE PRECISION array, dimension (LDU,UCOL)
								   UCOL = M if JOBZ = 'A' or JOBZ = 'O' and M < N;
								   UCOL = min(M,N) if JOBZ = 'S'.
								   If JOBZ = 'A' or JOBZ = 'O' and M < N, U contains the M-by-M
								   orthogonal matrix U;
								   if JOBZ = 'S', U contains the first min(M,N) columns of U
								   (the left singular vectors, stored columnwise);
								   if JOBZ = 'O' and M >= N, or JOBZ = 'N', U is not referenced. */
	int		LDU;				/* LDU is INTEGER
								   The leading dimension of the array U.  LDU >= 1; if
								   JOBZ = 'S' or 'A' or JOBZ = 'O' and M < N, LDU >= M. */
	double	*VT;				/* VT is DOUBLE PRECISION array, dimension (LDVT,N)
								   If JOBZ = 'A' or JOBZ = 'O' and M >= N, VT contains the
								   N-by-N orthogonal matrix V**T;
								   if JOBZ = 'S', VT contains the first min(M,N) rows of
								   V**T (the right singular vectors, stored rowwise);
								   if JOBZ = 'O' and M < N, or JOBZ = 'N', VT is not referenced. */
	int		LDVT;				/* LDVT is INTEGER
								   The leading dimension of the array VT.  LDVT >= 1;
								   if JOBZ = 'A' or JOBZ = 'O' and M >= N, LDVT >= N;
								   if JOBZ = 'S', LDVT >= min(M,N). */
	double	*WORK;				/* WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
								   On exit, if INFO = 0, WORK(1) returns the optimal LWORK; */
	int		LWORK;				/* LWORK is INTEGER
								   The dimension of the array WORK. LWORK >= 1.
								   If LWORK = -1, a workspace query is assumed.  The optimal
								   size for the WORK array is calculated and stored in WORK(1),
								   and no other work except argument checking is performed.

								   Let mx = max(M,N) and mn = min(M,N).
								   If JOBZ = 'N', LWORK >= 3*mn + max( mx, 7*mn ).
								   If JOBZ = 'O', LWORK >= 3*mn + max( mx, 5*mn*mn + 4*mn ).
								   If JOBZ = 'S', LWORK >= 4*mn*mn + 7*mn.
								   If JOBZ = 'A', LWORK >= 4*mn*mn + 6*mn + mx.
								   These are not tight minimums in all cases; see comments inside code.
								   For good performance, LWORK should generally be larger;
								   a query is recommended. */
	int		INFO;				/* INFO is INTEGER
								   = 0:  successful exit.
								   < 0:  if INFO = -i, the i-th argument had an illegal value.
								   > 0:  DBDSDC did not converge, updating process failed. */
	double WORK_QUERY = 0;

	uint8_t corr_ord = 0;			/* correlation order */
	int ym = 0;						/* correlation result matrix row dimension */
	int yn = 0;						/* correlation result matrix column dimension */

	int i = 0, j = 0;


	if (d == NULL) {
		return -1;
	}

	corr_ord = 4 * d->signals_n - 1;
	ym = x_n - corr_ord;	/* correlation result matrix row dimension */
	yn = corr_ord + 1;		/* correlation result matrix column dimension */
	d->svd.M = ym;
	d->svd.N = yn;

	/* Alloc all memory needed in SVD computations. */

	if (d->svd.A == NULL) {	
		A = calloc(ym * yn, sizeof(double));
		if (!A) {
			goto lm_m_l_q_fail_sys;
		}
		d->svd.A = A;
	}

	if (d->svd.A2 == NULL) {	
		A2 = calloc(ym * yn, sizeof(double));
		if (!A2) {
			goto lm_m_l_q_fail_sys;
		}
		d->svd.A2 = A2;
	}

	if (d->svd.S == NULL) {	
		S = calloc(lm_min(ym, yn), sizeof(double));
		if (!S) {
			goto lm_m_l_q_fail_sys;
		}
		d->svd.S = S;
	}

	if (d->svd.U == NULL) {	
		U = calloc(ym * ym, sizeof(double));
		if (!S) {
			goto lm_m_l_q_fail_sys;
		}
		d->svd.U = U;
	}

	if (d->svd.VT == NULL) {	
		VT = calloc(yn * yn, sizeof(double));
		if (!VT) {
			goto lm_m_l_q_fail_sys;
		}
		d->svd.VT = VT;
	}

	/* Call dgesvd_ with lwork = -1 to query optimal workspace size. */

	WORK_QUERY = 0;
	d->svd.LDA = ym;
	d->svd.LDU = ym;
	d->svd.LDVT = yn;
	d->svd.LWORK = 4 * ym * yn * ym * yn + 6 * ym * yn + lm_max(ym, yn);

	LWORK = -1;
	dgesvd_(&d->svd.JOBU, &d->svd.JOBVT, &d->svd.M, &d->svd.N, d->svd.A, &d->svd.LDA, d->svd.S, d->svd.U, &d->svd.LDU, d->svd.VT, &d->svd.LDVT, &WORK_QUERY, &LWORK, &INFO);
	if (INFO != 0) {
		if (INFO < 0) {
			/* Error on LAPACK's dgesvd_ query: the INFO-th argument had illegal value. */
			return -4;
		} else {
			/* Error on LAPACK's dgesvd_ query: "DBDSDC didn't converge, updating process failed. */
			return -5;
		}
	}

	LWORK = (int) WORK_QUERY;

	if (d->svd.WORK == NULL) {	
		WORK = calloc(LWORK, sizeof(double));
		if (!WORK) {
			goto lm_m_l_q_fail_sys;
		}
		d->svd.WORK = WORK;
	}

	d->svd.LWORK = LWORK;
	
	if (d->detections == NULL) {
		
		/* By default span all frequencies up to half sampling rate. */

		d->detections = calloc(d->fs / 2, sizeof(lm_detection_t));
		if (!d->detections) {
			goto lm_m_l_q_fail_sys;
		}

		i = 0;
		while (i < d->fs / 2) {
			d->detections[i].freq = i;
			++i;
		}
		
		d->detections_n = d->fs / 2;
	}

	if (d->a_real == NULL) {
		d->a_real = calloc(d->svd.N, sizeof(double));
		if (!d->a_real) {
			goto lm_m_l_q_fail_sys;
		}
	}

	if (d->a_imag == NULL) {
		d->a_imag = calloc(d->svd.N, sizeof(double));
		if (!d->a_imag) {
			goto lm_m_l_q_fail_sys;
		}
	}

	d->svd.prealloced = 1;

	return 0;

lm_m_l_q_fail_sys:

	if (A != NULL)
		free(A);
	if (A2 != NULL)
		free(A2);
	if (S != NULL)
		free(S);
	if (U != NULL)
		free(U);
	if (VT != NULL)
		free(VT);
	if (WORK != NULL)
		free(WORK);

	return -2;
}

int lm_detect(lm_detector_t *d, double *x, uint16_t x_n)
{
	int			i = 0, j = 0;
	int			INFO = 0;
	int			err = 0;
	uint16_t	corr_ord = 0;
	
	/* If want info. Peak detection */
	double t = 0;
	double omega = 0;
	double av_real = 0;	/* real result of a * V multiplication */
	double av_imag = 0;	/* imaginary result of a * V multiplication */
	uint16_t f = 0;
	double peak = 0;
	double peak_min = DBL_MAX, peak_max = DBL_MIN;
	double peak_min_freq = DBL_MAX, peak_max_freq = DBL_MIN;
	uint16_t ym = 0;
	uint16_t yn = 0;


	if (d == NULL) {
		return -1;
	}
	
	corr_ord = 4 * d->signals_n - 1;

	/* Need to query LAPACK if it hasn't been done yet. */

	if (d->svd.prealloced == 0) {
		
		d->svd.M = x_n - corr_ord;
		d->svd.N = corr_ord + 1;
		err = lm_lapack_query(d, x_n);
		if (err	!= 0) {
			return err;
		}
	}

	/* Compute correlation matrix. */	
	d->svd.A = lm_corr(x, corr_ord, x_n, d->svd.A);
	if (d->svd.A == NULL) {
		return -2;
	}
	
	ym = d->svd.M;
	yn = d->svd.N;
	
	/* LAPACKaize correlation matrix, LAPACK wants it colwise.
	 * In A2 matric store correlation values collwise. */
	for (i = 0; i < ym; ++i) {
		for (j = 0; j < d->svd.N; ++j) {
			d->svd.A2[j * ym + i] = d->svd.A[i * yn + j];
		}
	}
	
	if (!d->detections) {

		/* By default span all frequencies up to half sampling rate. */

		d->detections = calloc(d->fs / 2, sizeof(lm_detection_t));
		if (!d->detections) {
			return -3;
		}

		i = 0;
		while (f < d->fs / 2) {
			d->detections[f].freq = i;
			++i;
		}

		d->detections_n = d->fs / 2;
	}

	/* Compute SVD. */
	dgesvd_(&d->svd.JOBU, &d->svd.JOBVT, &d->svd.M, &d->svd.N, d->svd.A2, &d->svd.LDA, d->svd.S, d->svd.U, &d->svd.LDU, d->svd.VT, &d->svd.LDVT, d->svd.WORK, &d->svd.LWORK, &INFO);
	if (INFO != 0) {
		if (INFO < 0) {
			/* Error on LAPACK's dgesvd_ query: the INFO-th argument had illegal value. */
			return -4;
		} else {
			/* Error on LAPACK's dgesvd_ query: "DBDSDC didn't converge, updating process failed. */
			return -5;
		}
	}

	/* Run peak detection through all detection frequencies. */

	j = 0;
	peak_min = DBL_MAX;
	peak_max = DBL_MIN;

	while (j < d->detections_n) {

		f = d->detections[j].freq;

		i = 0;
		t = 0;
		while (i < yn) {
			omega = 2 * M_PI * f;
			d->a_real[i] = cos(omega * t);
			d->a_imag[i] = sin(omega * t);
			t += d->dt;
			++i;
		}

		i = 0;
		av_real = av_imag = 0;
		while (i < yn) {
			av_real += d->a_real[i] * (d->svd.VT[4 + i * yn] + d->svd.VT[5 + i * yn] + d->svd.VT[6 + i * yn] + d->svd.VT[7 + i * yn]);
			av_imag += d->a_imag[i] * (d->svd.VT[4 + i * yn] + d->svd.VT[5 + i * yn] + d->svd.VT[6 + i * yn] + d->svd.VT[7 + i * yn]);
			++i;
		}

		peak = 1 / (av_real * av_real + av_imag * av_imag);
		d->detections[j].peak = peak;

		if (peak > peak_max) {
			peak_max = peak;
			peak_max_freq = f;
		}
		if (peak < peak_min) {
			peak_min = peak;
			peak_min_freq = f;
		}

		if (f == 697) {
			d->info.peak_697 = peak;
		} else if (f == 770) {
			d->info.peak_770 = peak;
		} else if (f == 852) {
			d->info.peak_852 = peak;
		} else if (f == 941) {
			d->info.peak_941 = peak;
		}

		if (f == 1209) {
			d->info.peak_1209 = peak;
		} else if (f == 1336) {
			d->info.peak_1336 = peak;
		} else if (f == 1477) {
			d->info.peak_1477 = peak;
		} else if (f == 1633) {
			d->info.peak_1633 = peak;
		}

		j++;
	}

	d->info.peak_max = peak_max;
	d->info.peak_max_freq = peak_max_freq;
	d->info.peak_min = peak_min;
	d->info.peak_min_freq = peak_min_freq;


	return 0;

lm_d_m_fail_sys:
	
	return -2;
}

uint8_t lm_dtmf_decision(lm_detector_t *d)
{
	uint8_t dtmf_idx = 0;
	double peak_diff = 0;
	double peak = DBL_MIN;

	uint8_t i = 0;
	uint8_t peak_row_idx = 0;
	uint8_t peak_col_idx = 4;
	
	uint8_t peak_row_idx2 = 0;
	uint8_t peak_col_idx2 = 4;

	double row_val = 0;
	double col_val = 0;
	
	double row_val2 = 0;
	double col_val2 = 0;

	double row_2_peak = 0; 
	double col_2_peak = 0; 


	if (!d || !d->detections || d->detections_n < 1) {
		return 0;
	}

	/* Check if maximum peak is above peak threshold 2.
	peak_diff = d->info.peak_max - d->info.peak_min;
	if (peak_diff < d->peak_max_threshold2) {
		return 0;
	}*/
	
	/**
	 * Check if there are 2 DTMF frequencies with peak above peak threshold:
	 * one frequency must belong to the row frequencies, other to the column
	 * frequencies.
	 */

	/* Search for maximum peak in row frequencies. */
	i = 0;
	peak = d->info.peak_min;
	while (i < 4) {
		if (d->detections[i].peak > peak) {
			peak = d->detections[i].peak;
			peak_row_idx = i;
		} else {
			if (d->detections[i].peak > row_2_peak) {
				row_2_peak = d->detections[i].peak;
				peak_row_idx2 = i;
			}
		}
		++i;
	}

	/* Search for maximum peak in column frequencies. */
	peak = d->info.peak_min;
	while (i < 8) {
		if (d->detections[i].peak > peak) {
			peak = d->detections[i].peak;
			peak_col_idx = i - 4;
		} else {
			if (d->detections[i].peak > col_2_peak) {
				col_2_peak = d->detections[i].peak;
				peak_col_idx2 = i - 4;
			}
		}
		++i;
	}

	row_val = d->detections[peak_row_idx].peak / d->info.peak_min;
	col_val = d->detections[4 + peak_col_idx].peak / d->info.peak_min;
	row_val2 = d->detections[peak_row_idx2].peak / d->info.peak_min;
	col_val2 = d->detections[4 + peak_col_idx2].peak / d->info.peak_min;

#if DEBUG_LM	
	printf("lmusic: peak max=%f min=%f, max row/col: %f/%f max2 row/col: %f/%f, val row/col: %f/%f val2 row/col: %f/%f, digit candidate = '%c'\n",
			d->info.peak_max, d->info.peak_min,
			d->detections[peak_row_idx].peak, d->detections[4 + peak_col_idx].peak,
			d->detections[peak_row_idx2].peak, d->detections[4 + peak_col_idx2].peak,
			row_val, col_val, row_val2, col_val2,
			lm_dtmf_idx_2_char[1 + peak_row_idx * 4 + peak_col_idx]);
#endif

	/* Check if maximum peak is above peak threshold 3.
	peak_diff = d->info.peak_max - d->info.peak_min;
	if (peak_diff < d->peak_max_threshold3) {
		return 0;
	}*/

	/* Both peaks must be above threshold 2 and min peak must be below threshold 1.
	if ((d->detections[peak_row_idx].peak < d->peak_max_threshold2) 
			|| (d->detections[4 + peak_col_idx].peak < d->peak_max_threshold2)
			|| (d->info.peak_min > d->peak_max_threshold1)) {
		return 0;
	}*/

	/* Both peaks must be above threshold 2 and min peak must be below threshold 1. */
	if ((row_val < d->peak_max_threshold1) 
			|| (col_val < d->peak_max_threshold1)) {
		return 0;
	}
	
	/* Both peaks must be more than CROSSTALK times second peaks. */
	if ((row_val / row_val2 < LM_DTMF_CROSSTALK_THRESHOLD) 
			|| (col_val / col_val2 < LM_DTMF_CROSSTALK_THRESHOLD)) {
		return 0;
	}

	/**
	 * DTMF association to indices
	 *
	 * Tones:
	 *			1209 Hz	|	1336 Hz		| 1477 Hz	|	1633 Hz
	 *	697 Hz	'1'			'2'				'3'			'A'
	 *	770 Hz	'4'			'5'				'6'			'B'
	 *	852 Hz	'7'			'8'				'9'			'C'
	 *	941 Hz	'*'			'0'				'#'			'D'
	 *
	 *	Indices
	 *			1|2|3|4
	 *			5|6|7|8
	 *			9|10|11|12
	 *			13|14|15|16
	 */
	dtmf_idx = 1 + peak_row_idx * 4 + peak_col_idx;

	return dtmf_idx;
}

int lm_dtmf_detect(lm_detector_t *d, double *x, uint16_t x_n)
{
	uint8_t idx = 0;

	if (!d) {
		return -1;
	}

	if (lm_detect(d, x, x_n) != 0) {
		return -2;
	}

	return (int) lm_dtmf_decision(d);
}

/*
int lm_init(lm_detector_t *d, uint16_t fs, uint16_t x_n)
{
	int corr_ord = 4 * d->signals_n - 1;	/ Detector must have number of expected signals set. /
	int ym = x_n - corr_ord;				/ correlation result matrix row dimension, for 16 sample vector: 9 /
	int yn = corr_ord + 1;					/ correlation result matrix column dimension, for 16 sample vector: 8 /

	double	*VT = NULL;
	double	*WORK = NULL;
	double WORK_QUERY = 0;
	double *A = NULL;						/ Input to SVD: correlation matrix columnwise for LAPACK. /
	double *S = NULL;
	double *U = NULL;
	int INFO = 0;
		
	int i = 0, j = 0;


	if (!d) {
		return -1;
	}

	d->fs = fs;
	d->x_n = x_n;
	d->dt = dt;
	
	d->prealloced = 1;
	d->svd.JOBU = 'N';
	d->svd.JOBVT = 'A';
	d->svd.LWORK = 4 * M * N * M *N + 6 * M * N + lm_max(M, N);
	d->svd.M = ym;
	d->svd.N = yn;
	d->svd.LDVT = yn;
	d->svd.LDA = ym;
	d->svd.LDU = ym;

	/ Call dgesvd_ with lwork = -1 to query optimal workspace size. /
	
	/ prealloc space for the correlation matrix /
	A = calloc(ym * yn, sizeof(double));
	if (!A) {
		goto lm_d_m_d_i_fail_sys;
	}

	S = calloc(lm_min(ym, yn), sizeof(double));
	if (!S) {
		goto lm_d_m_d_i_fail_sys;
	}

	VT = calloc(yn * yn, sizeof(double));
	if (!VT) {
		goto lm_d_m_d_i_fail_sys;
	}

	LWORK = -1;
	dgesvd_(&d->svd.JOBU, &d->svd.JOBVT, &d->svd.M, &d->svd.N, A, &d->svd.LDA, S, U, &d->svd.LDU, VT, &d->svd.LDVT, &WORK_QUERY, &d->svd.LWORK, &INFO);
	if (INFO != 0) {
		if (INFO < 0) {
			/Error on LAPACK's dgesvd_ query: the INFO-th argument had illegal value. /
			return -3;
		} else {
			/Error on LAPACK's dgesvd_ query: "DBDSDC didn't converge, updating process failed. /
			return -4;
		}
	}

	d->svd.LWORK = (int) WORK_QUERY;
	WORK = calloc(d->svd.LWORK, sizeof(double));
	if (!WORK) {
		goto lm_d_m_d_i_fail_sys;
	}

	d->svd.WORK = WORK;
	d->svd.A = A;
	d->svd.U = NULL;
	d->svd.S = S;
	d->svd.VT = VT;

	if (d.want_info) {
		
		/ if want to know max and min peak /
		detections = calloc(fs/2, sizeof(lm_detection_t));
		if (!detections) {
			return -2;
		}

		i = 0;
		while (i < fs/2) {
			detections[i].freq = i;
			++i;
		}
	}
	
	d->detections = &detections[0];

	return 0;

lm_d_m_d_i_fail_sys:

	if (A != NULL)
		free(A);
	if (Y != NULL)
		free(Y);
	if (S != NULL)
		free(S);
	if (U != NULL)
		free(U);
	if (VT != NULL && !(d->want_evec))
		free(VT);
	return -2;
}
*/

void lm_deinit(lm_detector_t *d)
{
	if (!d) {
		return;
	}

	if (d->svd.A != NULL) {
		free(d->svd.A);
		d->svd.A = NULL;
	}
	if (d->svd.A2 != NULL) {
		free(d->svd.A2);
		d->svd.A2 = NULL;
	}
	if (d->svd.S != NULL) {
		free(d->svd.S);
		d->svd.S = NULL;
	}
	if (d->svd.U != NULL) {
		free(d->svd.U);
		d->svd.U = NULL;
	}
	if (d->svd.VT != NULL) {
		free(d->svd.VT);
		d->svd.VT = NULL;
	}
	if (d->svd.WORK != NULL) {
		free(d->svd.WORK);
		d->svd.WORK = NULL;
	}
	if (d->detections != NULL) {
		free(d->detections);
		d->detections = NULL;
	}
	if (d->a_real != NULL) {
		free(d->a_real);
		d->a_real = NULL;
	}
	if (d->a_imag != NULL) {
		free(d->a_imag);
		d->a_imag = NULL;
	}
}

int lm_standard_init(lm_detector_t *d, uint16_t fs, uint16_t x_n, uint8_t signals_n, lm_detection_t *detections, uint16_t detections_n)
{
	int f = 0;


	if (!d) {
		return -1;
	}
	
	memset(d, 0, sizeof(lm_detector_t));
	
	d->svd.prealloced = 0;
	d->svd.JOBU = 'N';
	d->svd.JOBVT = 'A';

	d->fs = fs;
	d->dt = 1.0 / (double) fs;
	d->x_n = x_n;
	
	d->signals_n = signals_n;		/* we are searching for @signals_n number of frequencies */

	if (!detections) {

		/* By default span all frequencies up to half sampling rate. */

		d->detections = calloc(fs / 2, sizeof(lm_detection_t));
		if (!d->detections) {
			return -2;
		}

		f = 0;
		while (f < fs / 2) {
			d->detections[f].freq = f;
			++f;
		}

		d->detections_n = fs / 2;
	} else {
		d->detections = detections;
		d->detections_n = detections_n;
	}

	/* Complete the setup. */
	return lm_lapack_query(d, x_n);
}

void lm_standard_deinit(lm_detector_t *d)
{
	lm_deinit(d);
}

int lm_dtmf_init(lm_detector_t *d, uint16_t fs, uint16_t x_n)
{
	int f = 0;


	if (!d) {
		return -1;
	}
	
	memset(d, 0, sizeof(lm_detector_t));

	d->svd.prealloced = 0;
	d->svd.JOBU = 'N';
	d->svd.JOBVT = 'A';
	
	d->fs = fs;
	d->dt = 1.0 / (double) fs;
	d->x_n = x_n;
	
	d->signals_n = 2;				/* we are searching for 2 frequencies */

	d->peak_max_threshold1 = LM_DTMF_PEAK_THRESHOLD1;
	d->peak_max_threshold2 = LM_DTMF_PEAK_THRESHOLD2;
	d->peak_max_threshold3 = LM_DTMF_PEAK_THRESHOLD3;

	/* Configure detection at DTMF frequencies. */

	d->detections = calloc(8, sizeof(lm_detection_t));
	if (!d->detections) {
		return -2;
	}

	f = 0;
	while (f < 8) {
		d->detections[f].freq = lm_dtmf_freqs[f];
		++f;
	}

	d->detections_n = 8;

	/* DTMF detection searches for any combination of the two DTMF frequencies. */
	return lm_lapack_query(d, x_n);
}

void lm_dtmf_deinit(lm_detector_t *d)
{
	lm_deinit(d);
}

const char *lm_get_version_string(void)
{
	/*
	 * Simply return the autotools generated string.
	 */
	return LM_VER_STRING;
}

unsigned int lm_get_version(void)
{
	uint32_t major = 0, minor = 0, micro = 0;
	uint32_t rv = 0;
	int parse_rv;
	
	/*
	 * Parse the autotools generated version.
	 */
	parse_rv = sscanf(LM_VERSION, "%u.%u.%u", &major, &minor, &micro);
	if (parse_rv != 3) {
		/*
		 * We're expected to parse all 3 version levels.
		 * If not, then this must not be an official release.
		 * Return all zeros on the version
		 */
		return 0;
	}
	
	/*
	 * We allow 8 bits for the major and minor, while
	 * allowing 16 bits for the micro.  16 bits for the micro
	 * may be beneficial for a continuous delivery model
	 * in the future.
	 */
	rv |= (major & 0xFF) << 24;
	rv |= (minor & 0xFF) << 16;
	rv |= micro & 0xFFFF;
	return rv;
}

int lm_transpose(double *a, int n)
{
	int i = 0, j = 0;
	double tmp = 0;

	while (i < n) {
		j = i;
		while (j < n) {
			tmp = a[i * n + j];
			a[i * n + j] = a[j * n + i];
			a[j * n + i] = tmp;
			++j;
		}
		++i;
	}
}

/* Testing. */

void lm_die(const char *reason, const char *file, int line)
{
	fprintf(stderr, "Test failed: %s. %s:%d\n", reason, file, line);
	exit(EXIT_FAILURE);
}

void lm_walk_32bit_arr_rowwise(uint32_t *arr, int m, int n, void (*cb)(uint32_t *))
{
	int i = 0, j = 0;

	while (i < m) {
		j = 0;
		while (j < n) {
			cb(&arr[i * n + j]);
			++j;
		}
		++i;
	}
}

void lm_walk_32bit_arr_colwise(uint32_t *arr, int m, int n, void (*cb)(uint32_t *))
{
	int i = 0, j = 0;

	while (j < n) {
		i = 0;
		while (i < m) {
			cb(&arr[i * n + j]);
			++i;
		}
		++j;
	}
}

void lm_walk_dbl_arr_rowwise(double *arr, int m, int n, void (*cb)(double *), void (*cb_row_end)(void))
{
	int i = 0, j = 0;

	while (i < m) {
		j = 0;
		while (j < n) {
			cb(&arr[i * n + j]);
			++j;
		}
		cb_row_end();
		++i;
	}
}

void lm_walk_dbl_arr_colwise(double *arr, int m, int n, void (*cb)(double *), void (*cb_col_end)(void))
{
	int i = 0, j = 0;

	while (j < n) {
		i = 0;
		while (i < m) {
			cb(&arr[i * n + j]);
			++i;
		}
		cb_col_end();
		++j;
	}
}

void lm_walk_dbl_arr2_rowwise(double **arr, int m, int n, void (*cb)(double *), void (*cb_row_end)(void))
{
	int i = 0, j = 0;

	while (i < m) {
		j = 0;
		while (j < n) {
			cb(&arr[i][j]);
			++j;
		}
		cb_row_end();
		++i;
	}
}

double lm_epsilon(void) {
	typedef union {
		long long i64;
		double d64;
	} dbl_64;

	dbl_64 s;

	s.d64 = 1.;
	s.i64++;
	return (s.d64 - 1.);
}

int lm_cmp_dbl_arr_rowwise(double *arr1, double *arr2, int m, int n)
{
	int i = 0, j = 0;
	double epsilon = lm_epsilon();

	while (i < m) {
		j = 0;
		while (j < n) {
			if (fabs(arr1[i * n + j] - arr2[i * n + j]) > epsilon) {
				return -1;
			}
			++j;
		}
		++i;
	}

	return 0;
}

int lm_cmp_dbl_arr_row_col(double *arr1, double *arr2, int m, int n)
{
	int i = 0, j = 0;
	double epsilon = 0.0000005;
	double a, b;

	fprintf(stderr, "epsilon: %.20f\n", epsilon);
	while (i < m) {
		j = 0;
		while (j < n) {
			a = arr1[i * n + j];
			b = arr2[i + j * n];
			fprintf(stderr, "[%d][%f]=?[%d][%f]\t[%.20f][%.20f]\n", i*n+j, arr1[i*n+j], i + j * n, arr2[i + j * n], a - b, fabs(a-b));
			if (fabs(arr1[i * n + j] - arr2[i + j * n]) > epsilon) {
				return -1;
			}
			++j;
		}
		++i;
	}

	return 0;
}

int lm_cmp_dbl_arr_diag_row(double *arr1, double *arr2, int m, int n)
{
	int i = 0;
	double epsilon = 0.0000001;
	double a, b;

	fprintf(stderr, "epsilon: %.20f\n", epsilon);
	while (i < m) {
			a = arr1[i * n + i];
			b = arr2[i];
			fprintf(stderr, "[%d][%f]=?[%d][%f]\t[%.20f][%.20f]\n", i*n+i, arr1[i*n+i], i, arr2[i], a - b, fabs(a-b));
			if (fabs(a - b) > epsilon)  {
				return -1;
			}
		++i;
	}

	return 0;
}
