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
 * Input is external file.
 */

#define Fs 8000
#define XN 16


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
	int dtmf = 0;
	uint16_t dtmf_row = 0;
	uint16_t dtmf_col = 0;


	/**
	 * Init standard MUSIC DTMF detector.
	 * Detects only 8 DTMF frequencies.
	 */
	lm_test(lm_dtmf_init(&d, Fs, XN) == 0, "lm_dtmf_init failed"); 

    if (argc == 2) {
        inFileName = argv[1];
    }
	
	printf("input file name = %s\n",inFileName);
 
    inFile = fopen(inFileName,"rb");
    if (!inFile) {
		printf("Cannot open input file %s\n",inFileName);
		goto lmt11_fail_sys;
    }

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

		dtmf = lm_dtmf_detect(&d, inputData, XN);
		if (dtmf < 0) {
			printf("lm_dtmf_detect failed with error %d\n", dtmf);
			goto lmt11_fail_sys;
		}
		lm_test(dtmf >= 0, "lm_dtmf_detect failed");

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

		if (dtmf > 0) {
			printf("DTMF decision=%d, DTMF digit=%c\n", dtmf, lm_dtmf_idx_2_char[dtmf]);
			dtmf_row = (uint16_t) (dtmf - 1) / 4;
			dtmf_col = (uint16_t) (dtmf - 1) % 4;  
			printf("row=%u, col=%u, freq1=%u, freq2=%u\n", dtmf_row, dtmf_col, lm_dtmf_row_freqs[dtmf_row], lm_dtmf_col_freqs[dtmf_col]); 
		}

		numWords = fread(intData, sizeof(int16_t), XN, inFile );
	}

    printf("\nFinished. sampleCount = %d\n",sampleCount);
    fclose(inFile);
	
	lm_dtmf_deinit(&d);

	return 0;

lmt11_fail_sys:
	fprintf(stderr, "Couldn't complete the test due to system error\n");
	return -1;
}
