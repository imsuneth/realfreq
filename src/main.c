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

static char * logfile = "realfreq_processed.log";
static char * reffile = NULL;
static char * outputfile = NULL;
static int is_bedmethyl = 0;

void initialize() {
    init_logger(logfile);
    init_mod(reffile);
}

void destroy() {
    destroy_mod();
    destroy_logger();
}

int main(int argc, char* argv[]) {
    //parse the user args
    const char* optstring = "yr:o:l:b";
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
            case 'l':
                logfile = optarg;
                break;
            case 'y':
                clear_log(true);
                break;
            default:
                fprintf(stderr, "Usage: %s [-r reference_file] [-o output_file] [-l processed_files_log]\n", argv[0]);
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

    set_output_file(outputfile, is_bedmethyl);

    initialize();
    read_files_from_stdin();
    destroy();
    return 0;
}
