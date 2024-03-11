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
#include "map.h"
#include "logger.h"
#include "misc.h"
#include "utils.h"
#include "mod.h"


void read_files_from_stdin() {
    char *filepath = (char *)malloc(FILEPATH_LEN * sizeof(char));
    while (1) {
        int status = fscanf(stdin, "%s", filepath);
        if (ferror(stdin) || feof(stdin)) {
            break;
        }
        read_file_contents(filepath);
        log_file_processed(filepath);
    }
}


void read_file_contents(char *filepath) {

    double realtime0 = realtime();

    opt_t opt;
    init_opt(&opt); //initialise options to defaults
    
    //initialise the core data structure
    core_t* core = init_core(filepath, opt, realtime0);

    meth_freq(core);

    //free the core data structure
    free_core(core,opt);

    return;
}

