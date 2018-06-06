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
 * Detect DTMF tones using lmusic.
 *
 * Input is external file.
 */

#define Fs 8000
#define XN 16


void usage(void)
{
	fprintf(stderr, "\nUsage:\n");
	fprintf(stderr, "lmtool <input file>\n\n");
	fprintf(stderr, "Example:\n\tlmtool fraction-131-140_16bit_s.raw\n\n");
}

typedef struct lmtool_detection_state_s {
	int dtmf_idx;
	size_t dtmf_start;
} lmtool_detection_state_t;

static void lmtool_update_detection_state(lmtool_detection_state_t *state, int dtmf_idx, size_t sample_n)
{
	if (dtmf_idx != state->dtmf_idx) {

		if (dtmf_idx == 0 || (state->dtmf_idx > 0)) {
			
			/**
			 * Record DTMF end.
			 */
			fprintf(stderr, "DTMF end. Digit=%c, sample=%zu, duration=%zums\n", lm_dtmf_idx_2_char[state->dtmf_idx], sample_n, (sample_n - state->dtmf_start) * 1000/Fs); 
			state->dtmf_start = 0;
		}

		if (dtmf_idx > 0) {

			/**
			 * Record DTMF start.
			 */
			fprintf(stderr, "DTMF start. Digit=%c, sample=%zu\n", lm_dtmf_idx_2_char[dtmf_idx], sample_n);
		
			state->dtmf_start = sample_n;
		}
		
		state->dtmf_idx = dtmf_idx;
	}
}

// convert 16 bit ints to doubles
void int2float(int16_t *input, double *output, int length)
{
    int i;
 
    for (i = 0; i < length; i++) {
        output[i] = (double)input[i] / ((double) 1.0);
    }
}

int main (int argc, char **argv)
{
    int16_t		data[XN];
    double		sample_vector[XN];
    size_t		word_n = 0;
    size_t		sample_n = 0;
    const char	*input_file_name = NULL;
    FILE		*fin = NULL;
	size_t		frame_n = 0;
	
	lm_detector_t				d = { 0 };
	int							dtmf_idx = 0;
	lmtool_detection_state_t	detection_state = { 0 };


    if (argc != 2) {
		goto lmtoolhelp;
	}

	input_file_name = argv[1];
	
	fprintf(stderr, "\nInput= %s\n", input_file_name);
 
    fin = fopen(input_file_name,"rb");
    if (!fin) {
		fprintf(stderr, "Cannot open input file %s\n", input_file_name);
		goto lmtoolerr;
    }
	
	/**
	 * Init standard MUSIC DTMF detector.
	 * Detects only 8 DTMF frequencies.
	 */

	if (lm_dtmf_init(&d, Fs, XN) != 0) {
		fprintf(stderr, "Cannot initialise lmusic descriptor.\n");
		goto lmtoolerr;
	}

    sample_n = 0;

    word_n = fread(data, sizeof(int16_t), XN, fin);

	while(word_n == XN)
	{
		frame_n++;

		int2float(data, sample_vector, word_n);

		dtmf_idx = lm_dtmf_detect(&d, sample_vector, XN);
		if (dtmf_idx < 0) {
			fprintf(stderr, "Error. lm_dtmf_detect failed with error %d\n", dtmf_idx);
			goto lmtoolerr;
		}

		/* Update detection state. */
		lmtool_update_detection_state(&detection_state, dtmf_idx, sample_n);
		
		sample_n += XN;

		/* Read next frame. */
		word_n = fread(data, sizeof(int16_t), XN, fin );
	}
	
	/* Run 0 through process to finish DTMF detection reports. */
	lmtool_update_detection_state(&detection_state, 0, sample_n);

    fprintf(stderr, "\nFinished. sampleCount = %zu\n", sample_n);
    fclose(fin);
	
	lm_dtmf_deinit(&d);

	return 0;

lmtoolhelp:
	usage();
	return EXIT_FAILURE;

lmtoolerr:
	return EXIT_FAILURE;
}
