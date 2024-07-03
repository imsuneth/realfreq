/**
 * @file meth.c
 * @brief methylation calling and frequency calculation

MIT License

Copyright (c) 2024 Suneth Samarasinghe

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
#include "minimod.h"

void simple_meth_view(core_t* core);
void meth_freq(core_t* core);
void init_meth(const char * ref);
void destroy_meth();
void dump_stats_map(const char * stats_file);
void print_stats(FILE * output_file, int is_bedmethyl);
int get_stats_len();
void load_stats_map(const char * dump_file);
void dump_stats_map(const char * dump_file);

char* get_stats_range(int start, int end);
char* get_stats_contig(const char* contig);
char* get_stats_contig_range(const char *contig, int start, int end);
char* get_stats_contig_range_mod_code(const char *contig, int start, int end, char mod_code);
