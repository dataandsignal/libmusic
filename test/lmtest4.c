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
 * Test LAPACK's dgesvd_.
 *
 * Example 2.
 *
 * DGESDD computes the singular value decomposition (SVD) of a real
 * M-by-N matrix A, optionally computing the left and right singular
 * vectors.  If singular vectors are desired, it uses a
 * divide-and-conquer algorithm.
 *
 * The SVD is written
 *
 *          A = U * SIGMA * transpose(V)
 *
 * where SIGMA is an M-by-N matrix which is zero except for its
 * min(m,n) diagonal elements, U is an M-by-M orthogonal matrix, and
 * V is an N-by-N orthogonal matrix.  The diagonal elements of SIGMA
 * are the singular values of A; they are real and non-negative, and
 * are returned in descending order.  The first min(m,n) columns of
 * U and V are the left and right singular vectors of A.
 *
 * Note that the routine returns VT = V**T, not V.
 *
 * The divide and conquer algorithm makes very mild assumptions about
 * floating point arithmetic. It will work on machines with a guard
 * digit in alm/subtract, or on those binary machines without guard
 * digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
 * Cray-2. It could conceivably fail on hexadecimal or decimal machines
 * without guard digits.
 */

/* Reference data. */
double ref_array_A[3][3] = {
	{ 1, 2, 3},
	{ 2, 4, 5 },
	{ 3, 5, 6 }
};

double ref_array_U[3][3] = {
	{ -0.327985, -0.736976, -0.591009 },
	{ -0.591009, -0.327985, 0.736976 },
	{ -0.736976, 0.591009, -0.327985 }
};

double ref_array_Sigma[3][1] = {
	{ 11.344814 },
	{ 0.515729 },
	{ 0.170915 }
};

double ref_array_VT[3][3] = {
	{ -0.327985, -0.591009, -0.736976 },
	{ 0.736976, 0.327985, -0.591009 },
	{ -0.591009, 0.736976, -0.327985 }
};

/* MATLAB result
 *
 *	>> A = [ 1, 2, 3; 2, 4, 5; 3, 5, 6]
 *
 *	A = 
 *		1     2     3
 *		2     4     5
 *      3     5     6
 *
 *	>> [U, S, V] = svd(A)
 *
 *	U =
 *      -0.3280   -0.7370   -0.5910
 *      -0.5910   -0.3280    0.7370
 *      -0.7370    0.5910   -0.3280
 *
 *	S =
 *		11.3448		0			0
 *      0			0.5157      0
 *      0			0			0.1709
 *
 *	V =
 *		-0.3280    0.7370   -0.5910
 *      -0.5910    0.3280    0.7370
 *      -0.7370   -0.5910   -0.3280
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
	/*	Arguments to LAPACK's dgesvd_.
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
	double	*A;					/* A is DOUBLE PRECISION array, dimension (LDA,N)
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

	/* helper variables */
	int i = 0, j = 0;
	double WORK_QUERY = 0;
	
	
	/* Call dgesvd_ with lwork = -1 to query optimal workspace size. */

	JOBU = 'A';
	JOBVT = 'A';
	M = 3;
	N = 3;
	LDA = 3;			/* (out) */
	LDU = 3;			/* (out) */
	S = NULL;			/* (don't care) */
	U = NULL;			/* (don't care) */
	VT = NULL;			/* (don't care) */
	LDVT = 3;			/* (out) */
	WORK = NULL;		/* (out) , because LWORK is 0 do not care */
	LWORK = 4 * M * N * M *N + 6 * M * N + lm_max(M, N);
	
	A = calloc(M * N, sizeof(double));
	if (!A) {
		goto lmt4_fail_sys;
	}
	for (i = 0; i < M; ++i) {
		for (j = 0; j < N; ++j) {
			A[j * M + i] = ref_array_A[i][j];
		}
	}

	S = calloc(lm_min(M, N), sizeof(double));
	if (!S) {
		goto lmt4_fail_sys;
	}
	
	U = calloc(LDU * M, sizeof(double));
	if (!U) {
		goto lmt4_fail_sys;
	}

	VT = calloc(LDVT * N, sizeof(double));
	if (!VT) {
		goto lmt4_fail_sys;
	}
	
	fprintf(stderr, "Reference array A:\n");
	lm_walk_dbl_arr_rowwise(&ref_array_A[0][0], M, N, cb_dbl, cb_dbl_row_end);
	
	fprintf(stderr, "Reference array U (rowwise):\n");
	lm_walk_dbl_arr_rowwise(&ref_array_U[0][0], M, M, cb_dbl, cb_dbl_row_end);
	
	fprintf(stderr, "Reference array Sigma (rowwise):\n");
	lm_walk_dbl_arr_rowwise(&ref_array_Sigma[0][0], lm_min(M, N), 1, cb_dbl, cb_dbl_row_end);

	fprintf(stderr, "Reference array V (rowwise):\n");
	lm_walk_dbl_arr_rowwise(&ref_array_VT[0][0], N, N, cb_dbl, cb_dbl_row_end);
	
	LWORK = -1;
	dgesvd_("A", "A", &M, &N, A, &LDA, S, U, &LDU, VT, &LDVT, &WORK_QUERY, &LWORK, &INFO);
	if (INFO != 0) {
		if (INFO < 0) {
			fprintf(stderr, "Error on LAPACK's dgesvd_ query: \"the %d-th argument had illegal value\"\n", INFO);
		} else {
			fprintf(stderr, "Error on LAPACK's dgesvd_ query: \"DBDSDC didn't converge, updating process failed\"\n");
		}
		return -1;
	}

	LWORK = (int) WORK_QUERY;
	WORK = calloc(LWORK, sizeof(double));
	if (!WORK) {
		goto lmt4_fail_sys;
	}

	fprintf(stderr, "LAPACK's dgesvd_ query optimal results: LDA %d, LDU %d, LDVT %d, LWORK %d, WORK_QUERY %f\n", LDA, LDU, LDVT, LWORK, WORK_QUERY);
	fprintf(stderr, "Rest of params: M %d, N %d\n", M, N);

	fprintf(stderr, "Matrix A submitted to LAPACK:\n");
	lm_walk_dbl_arr_rowwise(A, M, N, cb_dbl, cb_dbl_row_end);

	/* Compute SVD. */
	dgesvd_(&JOBU, &JOBVT, &M, &N, A, &LDA, S, U, &LDU, VT, &LDVT, WORK, &LWORK, &INFO);
	if (INFO != 0) {
		if (INFO < 0) {
			fprintf(stderr, "Error on LAPACK's dgesvd_ query: \"the %d-th argument had illegal value\"\n", INFO);
		} else {
			fprintf(stderr, "Error on LAPACK's dgesvd_ query: \"DBDSDC didn't converge, updating process failed\"\n");
		}
		return -1;
	}

	fprintf(stderr, "LAPACK's dgesvd_ SVD completed\n");

	fprintf(stderr, "Result A:\n");
	lm_walk_dbl_arr_rowwise(A, M, N, cb_dbl, cb_dbl_row_end);
	
	fprintf(stderr, "Result U (rowwise):\n");
	lm_walk_dbl_arr_rowwise(U, LDU, M, cb_dbl, cb_dbl_row_end);
	fprintf(stderr, "Result U (colwise):\n");
	lm_walk_dbl_arr_colwise(U, LDU, M, cb_dbl, cb_dbl_row_end);


	fprintf(stderr, "Result S (rowwise):\n");
	lm_walk_dbl_arr_rowwise(S, lm_min(M, N), 1, cb_dbl, cb_dbl_row_end);

	fprintf(stderr, "Result VT (rowwise):\n");
	lm_walk_dbl_arr_rowwise(VT, LDVT, N, cb_dbl, cb_dbl_row_end);
	fprintf(stderr, "Result VT (colwise):\n");
	lm_walk_dbl_arr_colwise(VT, LDVT, N, cb_dbl, cb_dbl_row_end);

	free(WORK);
	free(A);
	free(S);
	free(U);
	free(VT);

	return 0;

lmt4_fail_sys:
	fprintf(stderr, "Couldn't complete the test due to system error\n");
	return -1;
}
