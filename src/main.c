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
#include "file_handler.h"
#include "logger.h"
#include "meth.h"
#include "error.h"

static char * reffile = NULL;
static char * outputfile = NULL;
static int is_bedmethyl = 0;
static int is_resuming = 0;
static char * dumpfile = NULL;

void initialize() {
    init_meth(reffile, dumpfile, is_resuming);
    set_output_file(outputfile, is_bedmethyl);
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
    fprintf(fp, "  -s, --resume DUMP_FILE  resume from a dump file\n");
}

int main(int argc, char* argv[]) {
    //parse the user args
    const char* optstring = "r:o:b:s";
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
            case 's':
                is_resuming = 1;
                dumpfile = optarg;
            default:
                print_usage(stderr);
                exit(EXIT_FAILURE);
        }
    }

    if (reffile == NULL) {
        fprintf(stderr, "Reference file is not provided\n");
        exit(EXIT_FAILURE);
    }

    if (outputfile == NULL) {
        fprintf(stderr, "Output file is not provided\n");
        exit(EXIT_FAILURE);
    }

    if(is_resuming && dumpfile == NULL){
        fprintf(stderr, "Resuming but dump file is not provided\n");
        exit(EXIT_FAILURE);
    }

    initialize();
    read_files_from_stdin();
    destroy();
    return 0;
}
