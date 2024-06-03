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

void write_output(char *output_file, int is_bedmethyl) {
    FILE *output_fp = fopen(output_file, "w");
    if (output_fp == NULL) {
        ERROR("could not open the output file %s", output_file);
        exit(EXIT_FAILURE);
    }
    print_stats(output_fp, is_bedmethyl);
    fclose(output_fp);
    return;
}

