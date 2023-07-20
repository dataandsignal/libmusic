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


#ifndef LM_LM_H
#define LM_LM_H


#ifdef __cplusplus
extern "C" {
#endif


#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <inttypes.h>

#include <math.h>
#include <float.h>

#define LM_VER_STRING PACKAGE_STRING
#define LM_VERSION PACKAGE_VERSION

#define LM_DTMF_SIGNALS_EXPECTED_N 2
#define LM_DTMF_SIGNALS_N (2 * (LM_DTMF_SIGNALS_EXPECTED_N))								/* Need to account for complex signals. */
#define LM_CORR_ORD ((2 * (LM_DTMF_SIGNALS_N)) - 1)

#define LM_DTMF_PEAK_THRESHOLD 1000
#define LM_DTMF_PEAK_THRESHOLD1 500
#define LM_DTMF_PEAK_THRESHOLD2 20
#define LM_DTMF_PEAK_THRESHOLD3 90

#define LM_DTMF_CROSSTALK_THRESHOLD 4

/* 10dB */
#define LM_DTMF_FREQ_POWER_DIFF_THRESHOLD 3.1623

#define LM_DTMF_ENERGY_RATIO 0.95

typedef enum lm_status {
    lm_status_false = 0,
    lm_status_true = 1
} lm_status_t;

typedef struct lm_detection_s {
	uint16_t	freq;
	double		peak;
} lm_detection_t;

static const uint16_t lm_dtmf_freqs[] = {
	697, 770, 852, 941, 1209, 1336, 1477, 1633
};

static const uint16_t lm_dtmf_row_freqs[] = {
	697, 770, 852, 941
};

static const uint16_t lm_dtmf_col_freqs[] = {
	1209, 1336, 1477, 1633
};
 
/**
 * DTMF association to indices, used consistenly throughout lmusic.
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
static const char lm_dtmf_idx_2_char[17] = {
	'X', '1', '2', '3', 'A', '4', '5', '6', 'B', '7', '8', '9', 'C', '*', '0', '#', 'D'
};

typedef struct lm_detection_info_s {
	double peak_min;
	double peak_max;
	uint16_t peak_min_freq;
	uint16_t peak_max_freq;
	double peak_697;
	double peak_770;
	double peak_852;
	double peak_941;
	double peak_1209;
	double peak_1336;
	double peak_1477;
	double peak_1633;
} lm_detection_info_t;

typedef struct lm_detector_s {
	uint16_t		fs;				/* sampling frequency in Hz */
	uint16_t		x_n;			/* configured length of the detection block in samples */
	double			dt;				/* time between two consecutive samples in seconds = 1/@fs */

	struct lm_svd {
		uint8_t	prealloced;			/* 1 if memory and parameters needed by LAPACK are set.
									   lm_detect doesn't need to query LAPACK for the optimal values of parameters.
									   Otherise LAPACK is run twice: once to query for the params, and second to compute SVD.
									   If set to 0 all values in @svd array are ignored. */
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
		int		M;					/* M is INTEGER. The number of rows of the input matrix A.  M >= 0. Correlation result matrix row dimension. */
		int		N;					/* N is INTEGER. The number of columns of the input matrix A.  N >= 0. Correlation result matrix column dimension */
		double	*A, *A2;					/* A is DOUBLE PRECISION array, dimension (LDA,N)
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
	} svd;

	lm_detection_t	*detections;	/* Points to array of detections to be done. */
	uint8_t			signals_n;		/* Number of expected signals in the sample vector. */
	uint16_t		detections_n;	/* Number of detections of interest in @detections array. */

	double *a_real;
	double *a_imag;

	lm_detection_info_t	info;		/* If on entry @want_info is set to 1, then on exit @info contains details about result detection. */
	
	double	peak_max_threshold1;	/* Upper threshold on MIN, min peak must be below it. */
	double	peak_max_threshold2;	/* Lower threshold on MAX, both row and col max freqs must be above it. */
	double	peak_max_threshold3;	/* Upper threshold on MAX, only max peak must be above it. */
} lm_detector_t;

/**
 * Returns a covariance matrix unbiased estimator.
 *
 * Partition a length-N signal vector @s into nonoverlapping data segments (frames) of length
 * @n - @corr_ord.
 * Each data frame occupies one column of output matrix.
 * No windowing of the data is done. This improves the accuracy.
 * The estimate is unbiased. Result matrix has dimension @n - @corr_ord by @corr_ord + 1.
 *
 * @s - input signal samples vector
 * @corr_ord - correlation order
 * @n - number of samples in @s
 * @inline_where - if not NULL then result array is created at this memory, otherwise
 *					memory for result array is dynamically allocated
 *
 * TODO:	Aalm option for computing autocorrelation matrix instead of covariance.
 *			Convenient when a positive definite estimate is required. However, it is
 *			biased estimator and it windows the data.
 */
double* lm_corr(double *s, uint16_t corr_ord, uint16_t n, double *inline_where);

/**
 * Make a query to the LAPACK for optimal matrices dimensions.
 */
int lm_lapack_query(lm_detector_t *d, uint16_t x_n);

/**
 * Detect number of frequencies using MUSIC algorithm.
 *
 * Use this function for detection of any number of signals.
 *
 * Return 0 on success, negative values on error:
 *	-1 - bad input
 *	-2 - couldn't get memory
 *	-3 - generation of correlation matrix failed
 *	-4 - call to LAPACKS dgesvd_ quering optimal values failed: the INFO-th argument had illegal value.
 *	-5 - call to LAPACKS dgesvd_ quering optimal values failed: DBDSDC didn't converge, updating process failed.
 *	-6 - @want_info is 1, but @detections contains less than Fs/2 entries
 *
 * @x - (in) array of the samples
 * @x_n - (in) number of samples in the array @x
 * @d - (in/out) on entry describes what is to be detected, on exit contains result
 */
int lm_detect(lm_detector_t *d, double *x, uint16_t x_n);

/**
 * Make a decision if DTMF is detected.
 *
 * Return the index of the detected DTMFs. 0 if nothing has been detected.
 *
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
uint8_t lm_dtmf_decision(lm_detector_t *d);

/**
 * Detect DTMF frequencies using MUSIC algorithm.
 *
 * Use this function for detection of two frequencies.
 * Return the index of the detected DTMFs (1-16).
 * 0 if nothing has been detected. Negative value on error.
 *
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
 *	
 *
 * @x - array of the samples
 * @x_n - number of samples in the array @x
 */
int lm_dtmf_detect(lm_detector_t *d, double *x, uint16_t x_n);

void lm_deinit(lm_detector_t *d);

/**
 * Detect frequencies using MUSIC algorithm.
 *
 * Use this function for detection of any number of signals.
 */
int lm_standard_init(lm_detector_t *d, uint16_t fs, uint16_t x_n, uint8_t signals_n, lm_detection_t *detections, uint16_t detections_n);

void lm_standard_deinit(lm_detector_t *d);

/**
 * Init DTMF descriptor for detection with MUSIC algorithm.
 * Allocate memory for SVD computation, configure info structure,
 * setup signals number and dimensions of U, S, V matrices.
 *
 * Return 0 on success, negative values on error.
 *	-1 - bad input
 *	-2 - couldn't get memory
 *	-3 - call to LAPACKS dgesvd_ quering optimal values failed: the INFO-th argument had illegal value.
 *	-4 - call to LAPACKS dgesvd_ quering optimal values failed: DBDSDC didn't converge, updating process failed.
 *
 * @d - detector to be initialised
 * @fs - sampling frequency
 * @x_n - detection block in samples (number of samples in a sample vector passed to lm_dtmf_detect)
 */
int lm_dtmf_init(lm_detector_t *d, uint16_t fs, uint16_t x_n);

/**
 * Revert what lm_dtmf_init did: free memory.
 */
void lm_dtmf_deinit(lm_detector_t *d);

/* LAPACK */

/* dgesvd_ prototype */
extern void dgesvd_(char* jobu, char* jobvt, int* m, int* n, double* a, int* lda, double* s, double* u, int* ldu, double* vt, int* ldvt, double* work, int* lwork, int* info);

const char *lm_get_version_string(void);
unsigned int lm_get_version(void);
void lm_die(const char *s, const char *file, int line);

#define lm_circ_inc(p) p += 1; \
						p %= LM_LVL_GROUP_SIZE_MAX;

#define lm_test(x, m) if (!(x)) lm_die(m, __FILE__, __LINE__);

#define lm_min(x, y) ((x) < (y) ? (x) : (y))
#define lm_max(x, y) ((x) < (y) ? (y) : (x))

#define LM_FL_RCV 1u
#define LM_FL_SND 2u

#define lm_set_flag(p, f) (p)->flags |= (f)
#define lm_clear_flag(p, f) (p)->flags &= ~(f)
#define lm_test_flag(p, f) (p)->flags & (f)

/* Transpose square matrix @a in place. */
void lm_transpose(double *a, int n);

/* For testing. */

void lm_walk_32bit_arr_rowwise(uint32_t *arr, int m, int n, void (*cb)(uint32_t *));
void lm_walk_32bit_arr_colwise(uint32_t *arr, int m, int n, void (*cb)(uint32_t *));
void lm_walk_dbl_arr_rowwise(double *arr, int m, int n, void (*cb)(double *), void (*cb_row_end)(void));
void lm_walk_dbl_arr_colwise(double *arr, int m, int n, void (*cb)(double *), void (*cb_row_end)(void));

void lm_walk_dbl_arr2_rowwise(double **arr, int m, int n, void (*cb)(double *), void (*cb_row_end)(void));

/* Compute double EPSILON value on this machine. */
double lm_epsilon(void);

/* Cmp both arrays rowwise. */
int lm_cmp_dbl_arr_rowwise(double *arr1, double *arr2, int m, int n);

/* Cmp first array rowwise, second colwise. */
int lm_cmp_dbl_arr_row_col(double *arr1, double *arr2, int m, int n);

/* Cmp first array on diagonal, second rowwise. */
int lm_cmp_dbl_arr_diag_row(double *arr1, double *arr2, int m, int n);

#ifdef __cplusplus
}
#endif

#endif /* LM_LM_H */
