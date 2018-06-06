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
 * Test of DTMF decision making with MUSIC agorithm, lm_dtmf_detect.
 *
 * Input is DTMF "1" = { 697 Hz, 1209 Hz }. 8000 Hz, 16 samples: [10, 25].
 */

#define Fs 8000
#define XN 16

double ref_array_X[XN] = { -0.2071,   -0.7942,   -1.1108,   -0.6395,    0.5197,    1.6468,    1.9312,    1.1106,   -0.3025,   -1.3984,   -1.5514, -0.8580,    0.0096,    0.3922,    0.1753,   -0.1748 };


int main (int argc, char **argv)
{
	int dtmf = 0;


	lm_detector_t d;

	/**
	 * Init standard MUSIC DTMF detector.
	 * Detects only 8 DTMF frequencies.
	 */
	lm_test(lm_dtmf_init(&d, Fs, XN) == 0, "lm_dtmf_init failed"); 

	dtmf = lm_dtmf_detect(&d, ref_array_X, XN);
	if (dtmf != 1) {
		printf("lm_dtmf_detect failed got=%d but expected=1\n", dtmf);
		goto lmt10_fail_sys;
	}
	lm_test(dtmf == 1, "lm_dtmf_detect failed");

	/* Peak detection */

	printf(".info.peak_min = %f\n", d.info.peak_min);
	printf(".info.peak_max = %f\n", d.info.peak_max);
	printf(".info.peak_min_freq = %u\n", d.info.peak_min_freq);
	printf(".info.peak_max_freq = %u\n", d.info.peak_max_freq);
	printf(".info.peak_697 = %f\n", d.info.peak_697);
	printf(".info.peak_770 = %f\n", d.info.peak_770);
	printf(".info.peak_852 = %f\n", d.info.peak_852);
	printf(".info.peak_941 = %f\n", d.info.peak_941);
	printf(".info.peak_1209 = %f\n", d.info.peak_1209);
	printf(".info.peak_1336 = %f\n", d.info.peak_1336);
	printf(".info.peak_1477 = %f\n", d.info.peak_1477);
	printf(".info.peak_1633 = %f\n", d.info.peak_1633);

	lm_deinit(&d);

	return 0;

lmt10_fail_sys:
	fprintf(stderr, "Couldn't complete the test due to system error\n");
	return -1;
}
