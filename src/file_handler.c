/**
 * @file file_handler.c
 * @brief file handling
 * @author Suneth Samarasinghe (suneth@unsw.edu.au)

MIT License

Copyright (c) 2024 Suneth Samarasinghe (suneth@unsw.edu.au)

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.


******************************************************************************/


#include <stdio.h>
#include "file_handler.h"
#include "logger.h"
#include "misc.h"
#include "minimod.h"
#include "meth.h"
#include "error.h"

static const char *output_tsv;
static log_entry_t log_entry;

void set_output_file(const char *output_file) {
    output_tsv = output_file;
}


void write_output() {
    FILE *output_fp = fopen(output_tsv, "w");
    if (output_fp == NULL) {
        ERROR("could not open the output file %s", output_tsv);
        exit(EXIT_FAILURE);
    }
    print_stats(output_fp, 0);
    fclose(output_fp);
    return;
}

void read_file_contents(char *filepath) {

    double realtime0 = realtime();

    opt_t opt;
    init_opt(&opt); //initialise options to defaults
    
    //initialise the core data structure
    core_t* core = init_core(filepath, opt, realtime0);

    meth_freq(core);

    double realtime1 = realtime();

    write_output(output_tsv);

    double realtime2 = realtime();

    int stats_len = get_stats_len();
    
    log_entry.bamfile = filepath;
    log_entry.realtime_meth_freq = realtime1 - realtime0;
    log_entry.realtime_write_output = realtime2 - realtime1;
    log_entry.stats_len = stats_len;
    
    log_file_processed(&log_entry);

    //free the core data structure
    free_core(core,opt);
    return;
}

void read_files_from_stdin() {
    char *filepath = (char *)malloc(FILEPATH_LEN * sizeof(char));
    while (1) {
        int status = fscanf(stdin, "%s", filepath);
        if (ferror(stdin) || feof(stdin)) {
            break;
        }
        
        read_file_contents(filepath);
        
    }
    free(filepath);
}