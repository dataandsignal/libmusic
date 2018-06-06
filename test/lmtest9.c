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
 * Test of DTMF detection with MUSIC agorithm.
 * Test all results for DTMF frequencies are same:
 * 
 * 1. using lm_standard_init for standard detector but with frequencies limited to DTMF frequencies.
 * 2. using lm_standard_init for standard detector with default settings.
 * 3. using lm_dtmf_init for default DTMF detection
 *
 * Input is external file.
 */

#define Fs 8000
#define XN 16
#define SIGNALS_N 2


// convert 16 bit ints to doubles
void intToFloat(int16_t *input, double *output, int length)
{
    int i;
 
    for (i = 0; i < length; i++) {
        output[i] = (double)input[i] / ((double) 1.0);
    }
}

int main (int argc, char **argv)
{
	uint16_t x[XN];			/* samples */

	int i = 0, j = 0;

	int res;

	/* Testing data */
    int16_t intData[XN];
    double inputData[XN];
    int numWords = 0;
    int sampleCount = 0;
    const char *inFileName = "fraction-131-140_16bit_s.raw";
    FILE *inFile = NULL;
	uint32_t frame_n = 0;

	lm_detector_t d;
	lm_detection_t *detections = NULL;
	
	lm_detection_info_t info1;

	detections = calloc(8, sizeof(lm_detection_t));
	if (!detections) {
		printf("Can't get memory for testing\n");
		goto lmt9_fail_sys;
	}

	i = 0;
	while (i < 8) {
		detections[i].freq = lm_dtmf_freqs[i];
		++i;
	}

    if (argc == 2) {
        inFileName = argv[1];
    }
	
	printf("input file name = %s\n",inFileName);
 
    inFile = fopen(inFileName,"rb");
    if (!inFile) {
		printf("Cannot open input file %s\n",inFileName);
		goto lmt9_fail_sys;
    }
 
	/**
	 * Init standard MUSIC detector but with frequency range limited to 8 frequencies
	 * (DTMF frequencies) specified in @detections array.
	 * Detects only those 8 frequencies.
	 */
	lm_test(lm_standard_init(&d, Fs, XN, SIGNALS_N, detections, 8) == 0, "lm_standard_init failed"); 

	frame_n = 0;
    sampleCount = 0;

    numWords = fread(intData, sizeof(int16_t), XN, inFile);

	// until end of file
	while(numWords == XN)
	{
		frame_n++;
		printf("\nframe [%u]:\n", frame_n);

		intToFloat(intData, inputData, numWords);

		sampleCount += XN;
		printf("samples_n = [%d, %d)\n", sampleCount - XN, sampleCount);

		res = lm_detect(&d, inputData, XN);
		if (res != 0) {
			printf("lm_detect failed with error %d\n", res);
			goto lmt9_fail_sys;
		}
		lm_test(res == 0, "lm_detect failed");

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
		
		numWords = fread(intData, sizeof(int16_t), XN, inFile );
	}

    printf("\nFinished. sampleCount = %d\n",sampleCount);
    fclose( inFile );
	
	/* Copy detection info for next test. */
	memcpy(&info1, &d.info, sizeof(lm_detection_info_t));
	
	lm_deinit(&d);
	
	/* Run test again, now with detections == NULL, i.e. want all frequencies in range [0, Fs / 2). */
	
	printf("input file name = %s\n",inFileName);
 
    inFile = fopen(inFileName,"rb");
    if (!inFile) {
		printf("Cannot open input file %s\n",inFileName);
		goto lmt9_fail_sys;
    }

	memset(&d, 0, sizeof(d));

	/**
	 * Init standard MUSIC detector.
	 * Detects all frequencies in range [0, Fs / 2).
	 */
	lm_test(lm_standard_init(&d, Fs, XN, SIGNALS_N, NULL, 0) == 0, "lm_standard_init failed"); 

	frame_n = 0;
    sampleCount = 0;

    numWords = fread(intData, sizeof(int16_t), XN, inFile);

	// until end of file
	while(numWords == XN)
	{
		frame_n++;
		printf("\nframe [%u]:\n", frame_n);

		intToFloat(intData, inputData, numWords);

		sampleCount += XN;
		printf("samples_n = [%d, %d)\n", sampleCount - XN, sampleCount);

		res = lm_detect(&d, inputData, XN);
		if (res != 0) {
			printf("lm_detect failed with error %d\n", res);
			goto lmt9_fail_sys;
		}
		lm_test(res == 0, "lm_detect failed");

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
		
		numWords = fread(intData, sizeof(int16_t), XN, inFile );
	}

    printf("\nFinished. sampleCount = %d\n",sampleCount);
    fclose(inFile);

	/* Copy detection info for next test. */
	memcpy(&info1, &d.info, sizeof(lm_detection_info_t));

	/* Test compare. */
	lm_test(info1.peak_697 == d.info.peak_697, "peak error: 697 Hz (standard detector: default)");
	lm_test(info1.peak_770 == d.info.peak_770, "peak error: 770 Hz (standard detector: default)");
	lm_test(info1.peak_852 == d.info.peak_852, "peak error: 852 Hz (standard detector: default)");
	lm_test(info1.peak_941 == d.info.peak_941, "peak error: 941 Hz (standard detector: default)");
	lm_test(info1.peak_1209 == d.info.peak_1209, "peak error: 1209  Hz (standard detector: default)");
	lm_test(info1.peak_1336 == d.info.peak_1336, "peak error: 1336 Hz (standard detector: default)");
	lm_test(info1.peak_1477 == d.info.peak_1477, "peak error: 1477 Hz (standard detector: default)");
	lm_test(info1.peak_1633 == d.info.peak_1633, "peak error: 1633 Hz (standard detector: default)");
	
	lm_deinit(&d);
	
	/* Run test again, now using standard DTMF detector. */

	printf("input file name = %s\n",inFileName);
 
    inFile = fopen(inFileName,"rb");
    if (!inFile) {
		printf("Cannot open input file %s\n",inFileName);
		goto lmt9_fail_sys;
    }
 
	/**
	 * Init standard MUSIC DTMF detector.
	 * Detects only 8 DTMF frequencies.
	 */
	lm_test(lm_dtmf_init(&d, Fs, XN) == 0, "lm_dtmf_init failed"); 

	frame_n = 0;
    sampleCount = 0;

    numWords = fread(intData, sizeof(int16_t), XN, inFile);

	// until end of file
	while(numWords == XN)
	{
		frame_n++;
		printf("\nframe [%u]:\n", frame_n);

		intToFloat(intData, inputData, numWords);

		sampleCount += XN;
		printf("samples_n = [%d, %d)\n", sampleCount - XN, sampleCount);

		res = lm_detect(&d, inputData, XN);
		if (res != 0) {
			printf("lm_detect failed with error %d\n", res);
			goto lmt9_fail_sys;
		}
		lm_test(res == 0, "lm_detect failed");

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
		
		numWords = fread(intData, sizeof(int16_t), XN, inFile );
	}

    printf("\nFinished. sampleCount = %d\n",sampleCount);
    fclose( inFile );
	
	/* Test compare. */
	lm_test(info1.peak_697 == d.info.peak_697, "peak error: 697 Hz (standard DTMF detector)");
	lm_test(info1.peak_770 == d.info.peak_770, "peak error: 770 Hz (standard DTMF detector)");
	lm_test(info1.peak_852 == d.info.peak_852, "peak error: 852 Hz (standard DTMF detector)");
	lm_test(info1.peak_941 == d.info.peak_941, "peak error: 941 Hz (standard DTMF detector)");
	lm_test(info1.peak_1209 == d.info.peak_1209, "peak error: 1209  Hz (standard DTMF detector)");
	lm_test(info1.peak_1336 == d.info.peak_1336, "peak error: 1336 Hz (standard DTMF detector)");
	lm_test(info1.peak_1477 == d.info.peak_1477, "peak error: 1477 Hz (standard DTMF detector)");
	lm_test(info1.peak_1633 == d.info.peak_1633, "peak error: 1633 Hz (standard DTMF detector)");
	
	lm_deinit(&d);

	return 0;

lmt9_fail_sys:
	fprintf(stderr, "Couldn't complete the test due to system error\n");
	return -1;
}
