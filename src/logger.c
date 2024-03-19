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
#include "logger.h"

static const char *processed_files_log;

void init_logger(const char *filepath) {
    processed_files_log = filepath;
}

void destroy_logger() {
    // free((void *)processed_files_log);
}

void clear_log(bool clear) {
    if (clear) {
        FILE *log_file = fopen(processed_files_log, "w");
        if (log_file) {
            fclose(log_file);
        }
    }
}

void log_file_processed(const char *bamfile, double realtime_meth_freq, double realtime_write_output) {
    FILE *log_file = fopen(processed_files_log, "a");
    if (log_file) {
        fprintf(log_file, "%s\t%.3f sec\t%.3f sec\n", bamfile, realtime_meth_freq, realtime_write_output);
        fclose(log_file);
    } else {
        fprintf(stderr, "Error: cannot open log file %s\n", processed_files_log);
    }
    
}

