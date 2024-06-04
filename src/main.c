/**
 * @file mod.h
 * @brief modification tags

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
#include <getopt.h>
#include <stdbool.h>
#include <stdlib.h>
#include <unistd.h>
#include "file_handler.h"
#include "logger.h"
#include "meth.h"
#include "error.h"
#include "ref.h"
#include "misc.h"

static char * reffile = NULL;
static char * outputfile = NULL;
static int is_bedmethyl = 0;
static int is_resuming = 0;
static char * dumpfile = NULL;

void initialize() {
    init_meth(reffile);
    if(is_resuming) {
        if(access(dumpfile, F_OK) == -1) {
            fprintf(stderr, "Dump file does not exist\n");
        }
        fprintf(stderr, "[realfreq] resuming, loading stats map from %s\n", dumpfile);
        load_stats_map(dumpfile);
        write_output(outputfile, is_bedmethyl);
    }
    load_ref(reffile);
    
}

void destroy() {
    destroy_meth();
}

void print_usage(FILE * fp) {
    fprintf(fp, "Usage: realfreq [options]\n");
    fprintf(fp, "Options:\n");
    fprintf(fp, "  -r, --reference FILE    reference file\n");
    fprintf(fp, "  -o, --output FILE       methylation frequency output file\n");
    fprintf(fp, "  -b, --bedmethyl         output in bedMethyl format\n");
    fprintf(fp, "  -d, --dump FILE         dump file\n");
    fprintf(fp, "  -s, --resume            resume from a dump file\n");
}

void read_file_contents(char *filepath) {

    double realtime0 = realtime();

    opt_t opt;
    init_opt(&opt); //initialise options to defaults
    
    //initialise the core data structure
    core_t* core = init_core(filepath, opt, realtime0);

    meth_freq(core);

    double realtime1 = realtime();

    dump_stats_map(dumpfile);
    write_output(outputfile, is_bedmethyl);

    double realtime2 = realtime();

    static log_entry_t log_entry;
    log_entry.bamfile = filepath;
    log_entry.realtime_meth_freq = realtime1 - realtime0;
    log_entry.realtime_write_output = realtime2 - realtime1;
    log_entry.stats_len = get_stats_len();
    log_file_processed(&log_entry);

    //free the core data structure
    free_core(core,opt);
    return;
}

int main(int argc, char* argv[]) {
    //parse the user args
    const char* optstring = "r:o:bd:s";
    int opt;
    while ((opt = getopt(argc, argv, optstring)) != -1) {
        switch (opt) {
            case 'r':
                reffile = optarg;
                break;
            case 'o':
                outputfile = optarg;
                break;
            case 'b':
                is_bedmethyl = 1;
                break;
            case 'd':
                dumpfile = optarg;
                break;
            case 's':
                is_resuming = 1;
                break;
            default:
                print_usage(stderr);
                exit(EXIT_FAILURE);
        }
    }

    if (reffile == NULL) {
        fprintf(stderr, "Reference file is not specified. Use -r or --reference flag\n");
        exit(EXIT_FAILURE);
    }

    if (outputfile == NULL) {
        fprintf(stderr, "Output file is not specified. Use -o or --output flag\n");
        exit(EXIT_FAILURE);
    }

    if (dumpfile == NULL) {
        fprintf(stderr, "Dump file is not specified. Use -d or --dump flag\n");
        exit(EXIT_FAILURE);
    }

    initialize();
    
    char *filepath = (char *)malloc(FILEPATH_LEN * sizeof(char));
    fprintf(stderr, "[realfreq] reading file path from stdin\n");
    while (1) {
        int status = fscanf(stdin, "%s", filepath);
        if (ferror(stdin)) {
            fprintf(stderr, "[realfreq] error reading from stdin\n");
            fprintf(stderr, "[realfreq] exiting...\n");
            break;
        }
        if (feof(stdin)) {
            fprintf(stderr, "[realfreq] end of stdin\n");
            fprintf(stderr, "[realfreq] exiting...\n");
            break;
        }
        
        fprintf(stderr, "[realfreq] processing file %s\n", filepath);
        read_file_contents(filepath);
        
    }
    free(filepath);

    destroy();
    return 0;
}
