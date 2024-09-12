/**
 * @file mod.h
 * @brief modification tags

MIT License

Copyright (c) 2023 Hasindu Gamaarachchi (hasindu@unsw.edu.au)

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

#ifndef PULSE_H
#define PULSE_H

#include "minimod.h"

uint16_t *get_mod_tag(bam1_t *record, char *tag, uint32_t *len_ptr);
const char *get_mm_tag_ptr(bam1_t *record);
uint8_t *get_ml_tag(bam1_t *record, uint32_t *len_ptr);
void modbases_single(core_t* core, db_t* db, int32_t i);
void update_freq_map(core_t* core, db_t* db);
void print_freq_tsv_header(opt_t opt);
void print_freq_output(opt_t opt);
void print_view_output(core_t* core, db_t* db);
void init_freq_map();
void destroy_freq_map();
char* get_stats_range(int start, int end);
char* get_stats_contig(const char* contig);
char* get_stats_contig_range(const char *contig, int start, int end);
char* get_stats_contig_range_mod_code(const char *contig, int start, int end, char mod_code);
void load_stats_map(const char * dump_file);
void dump_stats_map(const char * dump_file);
uint8_t parse_mod_codes(const char* mod_codes_str);
void parse_mod_threshes(const char* mod_codes_str, char* mod_thresh_str, uint8_t n_codes);
void print_mod_options(opt_t opt);

#endif
