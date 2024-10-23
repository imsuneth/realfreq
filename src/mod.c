/**
 * @file pulse.c
 * @brief modification tags

MIT License

Copyright (c) 2023 Hasindu Gamaarachchi
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
#define _GNU_SOURCE
#include "mod.h"
#include "misc.h"
#include "error.h"
#include "khash.h"
#include "ref.h"
#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <stdarg.h>
#include <ctype.h>
#include <stdbool.h>
#include <pthread.h>

static const int valid_bases[256] = {
    ['A'] = 1, ['C'] = 1, ['G'] = 1, ['T'] = 1, ['U'] = 1, ['N'] = 1,
    ['a'] = 1, ['c'] = 1, ['g'] = 1, ['t'] = 1, ['u'] = 1, ['n'] = 1
};

static const int valid_strands[256] = {
    ['+'] = 1,
    ['-'] = 1
};

static const int base_idx_lookup[256] = {
    ['A'] = 0,
    ['C'] = 1,
    ['G'] = 2,
    ['T'] = 3,
    ['U'] = 4,
    ['N'] = 5,
    ['a'] = 0,
    ['c'] = 1,
    ['g'] = 2,
    ['t'] = 3,
    ['u'] = 4,
    ['n'] = 5,
};

static const char base_complement_lookup[256] = {
    ['A'] = 'T',
    ['C'] = 'G',
    ['G'] = 'C',
    ['T'] = 'A',
    ['U'] = 'A',
    ['N'] = 'N',
    ['a'] = 't',
    ['c'] = 'g',
    ['g'] = 'c',
    ['t'] = 'a',
    ['u'] = 'a',
    ['n'] = 'n',
};

static const uint8_t valid_mod_codes[256] = {
    ['a'] = 1, ['b'] = 1, ['c'] = 1, ['e'] = 1, ['f'] = 1, ['g'] = 1, ['h'] = 1, ['m'] = 1, ['n'] = 1, ['o'] = 1, 
    ['A'] = 1, ['C'] = 1, ['G'] = 1, ['T'] = 1, ['U'] = 1, ['N'] = 1
};

// required mods for the mod freq calculation (255 means unset)
static uint8_t req_mods[256] = {
    ['a'] = 255, ['b'] = 255, ['c'] = 255, ['e'] = 255, ['f'] = 255, ['g'] = 255, ['h'] = 255, ['m'] = 255, ['n'] = 255, ['o'] = 255, 
    ['A'] = 255, ['C'] = 255, ['G'] = 255, ['T'] = 255, ['U'] = 255, ['N'] = 255
};

// threshold for the mod freq calculation (51 = 0.2 * 255)
static uint8_t req_threshes[256] = {
    ['a'] = 51, ['b'] = 51, ['c'] = 51, ['e'] = 51, ['f'] = 51, ['g'] = 51, ['h'] = 51, ['m'] = 51, ['n'] = 51, ['o'] = 51, 
    ['A'] = 51, ['C'] = 51, ['G'] = 51, ['T'] = 51, ['U'] = 51, ['N'] = 51
};

static inline int die(const char *format, ...) {
    va_list args;
    va_start(args, format);
    vfprintf(stderr, format, args);
    va_end(args);

    exit(EXIT_FAILURE);
}

const char *get_mm_tag_ptr(bam1_t *record){

    const char* tag = "MM";
    // get the mm
    uint8_t *data = bam_aux_get(record, tag);
    if(data == NULL){
        LOG_TRACE("%s tag not found in read %s",tag, bam_get_qname(record));
        return NULL;
    }

    const char *mm_str = bam_aux2Z(data);
    if(mm_str == NULL){
        LOG_TRACE("%s tag could not be decoded for %s. Is it type Z?",tag, bam_get_qname(record));
        return NULL;
    }

    return mm_str;
}

uint8_t *get_ml_tag(bam1_t *record, uint32_t *len_ptr){

    const char* tag = "ML";
    // get the mm
    uint8_t *data = bam_aux_get(record, tag);
    if(data == NULL){
        LOG_TRACE("%s tag not found in read %s",tag, bam_get_qname(record));
        return NULL;
    }

    // check if the type of the tag is array of integers
    const char aux_type = data[0];
    if (aux_type != 'B') {
        LOG_TRACE("%s tag is not of type B in read %s",tag, bam_get_qname(record));
        return NULL;
    }

    // get the array len
    uint32_t len = bam_auxB_len(data);
    if(len == 0){
        LOG_TRACE("%s tag array length is 0 in read %s",tag, bam_get_qname(record));
        return NULL;
    }

    // check if the array type is uint8_t
    const char aux_array_type = data[1];
    if (aux_array_type != 'C') {
        LOG_TRACE("%s array tag type '%c' is not of type C in read %s",tag, aux_array_type, bam_get_qname(record));
        return NULL;
    }

    // get the actual stuff
    uint8_t *array = (uint8_t *)malloc(sizeof(uint8_t)*len); //can be optimised for premalloced array as arg
    MALLOC_CHK(array);

    for(int i=0;i<len;i++){
        array[i] = bam_auxB2i(data,i);
    }

    //set the length
    *len_ptr = len;

    return array;

}

uint8_t parse_mod_codes(const char* mod_codes_str){
    uint8_t n_codes = 0, i=0;

    if(strlen(mod_codes_str)==0){
        ERROR("%s","Modification code(s) provided with -c option is empty");
        exit(EXIT_FAILURE);
    }

    char c;
    while((c = mod_codes_str[i]) != '\0'){
        if(valid_mod_codes[(int)c]==0){
            ERROR("Invalid modification code %c",c);
            exit(EXIT_FAILURE);
        }
        if(req_mods[(int)c]!=255){
            ERROR("Duplicate modification code %c",c);
            exit(EXIT_FAILURE);
        }
        req_mods[(int)c] = i;
        i++;
    }
    n_codes = i;

    return n_codes;
}

void parse_mod_threshes(const char* mod_codes_str, char* mod_thresh_str, uint8_t n_codes){
    uint8_t i=0;

    if(strlen(mod_thresh_str)==0){
        ERROR("%s","Modification threshold(s) provided with -m option is empty");
        exit(EXIT_FAILURE);
    }

    char *mod_thresh_str_copy = (char *)malloc(strlen(mod_thresh_str)+1);
    MALLOC_CHK(mod_thresh_str_copy);
    strcpy(mod_thresh_str_copy, mod_thresh_str);
    char* token = strtok(mod_thresh_str_copy, ",");
    i=0;
    while(token!=NULL){
        char *end;
        errno = 0;
        double d = strtod(token, &end);

        if (errno != 0 || end == token || *end != '\0') {
            ERROR("Invalid modification threshold %s",token);
            exit(EXIT_FAILURE);
        }
        
        if(d<0 || d>1){
            ERROR("Modification threshold should be in the range 0.0 to 1.0. You entered %f",d);
            exit(EXIT_FAILURE);
        }

        req_threshes[(int)mod_codes_str[i]] = d*255;
        token = strtok(NULL, ",");
        i++;
    }
    if(i!=n_codes){
        ERROR("Number of modification codes and thresholds do not match. Codes:%d, Thresholds:%d",n_codes,i);
        exit(EXIT_FAILURE);
    }
    free(mod_thresh_str_copy);
}

void print_mod_options(opt_t opt){
    for(int i=0;i<opt.n_mods;i++){
        INFO("Modification code: %c, Threshold: %f", opt.mod_codes_str[i], req_threshes[(int)opt.mod_codes_str[i]]/255.0);
    }
}

char* make_key(const char *chrom, int pos, char mod_code, char strand){
    int start_strlen = snprintf(NULL, 0, "%d", pos);
    int key_strlen = strlen(chrom) + start_strlen  + 6;
    
    char* key = (char *)malloc(key_strlen * sizeof(char));
    MALLOC_CHK(key);
    snprintf(key, key_strlen, "%s\t%d\t%c\t%c", chrom, pos, mod_code, strand);
    return key;
}

void decode_key(char *key, char **contig, int *pos, char *mod_code, char *strand){
    char* key_copy = (char *)malloc(strlen(key)+1);
    MALLOC_CHK(key_copy);
    strcpy(key_copy, key);
    char* token = strtok(key_copy, "\t");
    *contig = calloc(strlen(token)+1, sizeof(char));
    MALLOC_CHK(*contig);
    strcpy(*contig, token);
    *pos = atoi(strtok(NULL, "\t"));
    *mod_code = strtok(NULL, "\t")[0];
    *strand = strtok(NULL, "\t")[0];
    free(key_copy);
}

char* get_stats_contig_range_mod_code(const char *contig, int start, int end, char mod_code, khash_t(freqm) *freq_map) {
    freq_t ** freqs = (freq_t **)malloc(sizeof(freq_t *)*kh_size(freq_map));
    MALLOC_CHK(freqs);
    int len = 0;
    for (khiter_t k = kh_begin(freq_map); k != kh_end(freq_map); ++k) {
        if (kh_exist(freq_map, k)) {
            freq_t * freq = kh_value(freq_map, k);
            char *contig1 = NULL;
            int start1;
            char mod_code1;
            char strand1;
            char * key = (char *) kh_key(freq_map, k);
            decode_key(key, &contig1, &start1, &mod_code1, &strand1);
            if(strcmp(contig1, contig) == 0 && start1 >= start && start1 <= end && mod_code1 == mod_code){
                freqs[len] = freq;
                freq->key = key;
                len++;
            }
        }
    }
    if(len == 0){
        free(freqs);
        return NULL;
    }
    char *all_stats = (char *)malloc(sizeof(char)*len*100);
    MALLOC_CHK(all_stats);

    for(int i=0;i<len;i++){
        freq_t *freq = freqs[i];
        double freq_value = (double)freq->n_mod*100/freq->n_called;
        char *contig1 = NULL;
        int start1;
        char mod_code1;
        char strand1;
        char * key = freq->key;
        decode_key(key, &contig1, &start1, &mod_code1, &strand1);
        sprintf(all_stats, "%s%s\t%d\t%d\t%c\t%d\t%d\t%f\t%c\n", all_stats, contig1, start1, start1, strand1, freq->n_called, freq->n_mod, freq_value, mod_code1);
    }

    free(freqs);
    return all_stats;
}

char* get_stats_range(int start, int end, khash_t(freqm) *freq_map) {
    freq_t ** freqs = (freq_t **)malloc(sizeof(freq_t *)*kh_size(freq_map));
    MALLOC_CHK(freqs);
    int len = 0;
    for (khiter_t k = kh_begin(freq_map); k != kh_end(freq_map); ++k) {
        if (kh_exist(freq_map, k)) {
            freq_t * freq = kh_value(freq_map, k);
            char *contig1 = NULL;
            int start1;
            char mod_code1;
            char strand1;
            char * key = (char *) kh_key(freq_map, k);
            decode_key(key, &contig1, &start1, &mod_code1, &strand1);
            if(start1 >= start && start1 <= end){
                freqs[len] = freq;
                freq->key = key;
                len++;
            }
        }
    }
    if(len == 0){
        free(freqs);
        return NULL;
    }
    char *all_stats = (char *)malloc(sizeof(char)*len*100);
    MALLOC_CHK(all_stats);

    for(int i=0;i<len;i++){
        freq_t *freq = freqs[i];
        double freq_value = (double)freq->n_mod*100/freq->n_called;
        char *contig1 = NULL;
        int start1;
        char mod_code1;
        char strand1;
        char * key = freq->key;
        decode_key(key, &contig1, &start1, &mod_code1, &strand1);
        sprintf(all_stats, "%s%s\t%d\t%d\t%c\t%d\t%d\t%f\t%c\n", all_stats, contig1, start1, start1, strand1, freq->n_called, freq->n_mod, freq_value, mod_code1);
    }

    free(freqs);
    return all_stats;

}

char* get_stats_contig(const char* contig, khash_t(freqm) *freq_map) {
    freq_t ** freqs = (freq_t **)malloc(sizeof(freq_t *)*kh_size(freq_map));
    MALLOC_CHK(freqs);
    int len = 0;
    for (khiter_t k = kh_begin(freq_map); k != kh_end(freq_map); ++k) {
        if (kh_exist(freq_map, k)) {
            freq_t * freq = kh_value(freq_map, k);
            char *contig1 = NULL;
            int start1;
            char mod_code1;
            char strand1;
            char * key = (char *) kh_key(freq_map, k);
            decode_key(key, &contig1, &start1, &mod_code1, &strand1);
            if(strcmp(contig1, contig) == 0){
                freqs[len] = freq;
                freq->key = key;
                len++;
            }
        }
    }
    if(len == 0){
        free(freqs);
        return NULL;
    }
    char *all_stats = (char *)malloc(sizeof(char)*len*500);
    MALLOC_CHK(all_stats);

    for(int i=0;i<len;i++){
        freq_t *freq = freqs[i];
        double freq_value = (double)freq->n_mod*100/freq->n_called;
        char *contig1 = NULL;
        int start1;
        char mod_code1;
        char strand1;
        char * key = freq->key;
        decode_key(key, &contig1, &start1, &mod_code1, &strand1);
        sprintf(all_stats, "%s%s\t%d\t%d\t%c\t%d\t%d\t%f\t%c\n", all_stats, contig1, start1, start1, strand1, freq->n_called, freq->n_mod, freq_value, mod_code1);
    }

    free(freqs);
    return all_stats;

}

char* get_stats_contig_range(const char *contig, int start, int end, khash_t(freqm) *freq_map) {
    freq_t ** freqs = (freq_t **)malloc(sizeof(freq_t *)*kh_size(freq_map));
    MALLOC_CHK(freqs);
    int len = 0;
    for (khiter_t k = kh_begin(freq_map); k != kh_end(freq_map); ++k) {
        if (kh_exist(freq_map, k)) {
            freq_t * freq = kh_value(freq_map, k);
            char *contig1 = NULL;
            int start1;
            char mod_code1;
            char strand1;
            char * key = (char *) kh_key(freq_map, k);
            decode_key(key, &contig1, &start1, &mod_code1, &strand1);
            if(strcmp(contig1, contig) == 0 && start1 >= start && start1 <= end){
                freqs[len] = freq;
                freq->key = key;
                len++;
            }
        }
    }
    if(len == 0){
        free(freqs);
        return NULL;
    }
    char *all_stats = (char *)malloc(sizeof(char)*len*100);
    MALLOC_CHK(all_stats);

    for(int i=0;i<len;i++){
        freq_t *freq = freqs[i];
        double freq_value = (double)freq->n_mod*100/freq->n_called;
        char *contig1 = NULL;
        int start1;
        char mod_code1;
        char strand1;
        char * key = freq->key;
        decode_key(key, &contig1, &start1, &mod_code1, &strand1);
        sprintf(all_stats, "%s%s\t%d\t%d\t%c\t%d\t%d\t%f\t%c\n", all_stats, contig1, start1, start1, strand1, freq->n_called, freq->n_mod, freq_value, mod_code1);
    }

    free(freqs);
    return all_stats;
}


// static freq_t ** get_stats(khash_t(freqm)* freq_map, uint32_t *meth_freqs_len){
//     uint32_t len = 0;
//     freq_t ** freqs = (freq_t **)malloc(sizeof(freq_t *)*kh_size(freq_map));
//     MALLOC_CHK(freqs);
//     for (khiter_t k = kh_begin(freq_map); k != kh_end(freq_map); ++k) {
//         if (kh_exist(freq_map, k)) {
//             freq_t * freq = kh_value(freq_map, k);
//             freqs[len] = freq;
//             len++;
//         }
//     }
//     *meth_freqs_len = len;
//     return freqs;
// }

void process_tsv_file(const char *tsv_file, opt_t opt, khash_t(freqm) *freq_map){
    // columns: chromosome start end read_name log_lik_ratio log_lik_methylated log_lik_unmethylated num_calling_strands num_cpgs sequence
    FILE *fp = fopen(tsv_file, "r");
    if(fp == NULL){
        ERROR("Could not open file %s", tsv_file);
        exit(EXIT_FAILURE);
    }

    char *line = NULL;
    size_t len = 0;
    size_t read;
    //skip the header
    read = getline(&line, &len, fp);
    while ((read = getline(&line, &len, fp)) != -1) {
        char *contig = strtok(line, "\t");
        int start = atoi(strtok(NULL, "\t"));
        strtok(NULL, "\t"); // end
        strtok(NULL, "\t"); // read_name
        double log_lik_ratio = atof(strtok(NULL, "\t"));
        
        char *key = make_key(contig, start, 'm', '+');           
        khiter_t k = kh_get(freqm, freq_map, key);
        if (k == kh_end(freq_map)) { // not found, add to map
            freq_t * freq = (freq_t *)malloc(sizeof(freq_t));
            MALLOC_CHK(freq);
            freq->n_called = 1;
            freq->n_mod = log_lik_ratio >= 0 ? 1 : 0;
            int ret;
            k = kh_put(freqm, freq_map, key, &ret);
            kh_value(freq_map, k) = freq;
        } else { // found, update the values
            free(key);
            freq_t * freq = kh_value(freq_map, k);
            freq->n_called += 1;
            freq->n_mod += log_lik_ratio >= 0 ? 1 : 0;
        }
    }
    free(line);
    fclose(fp);
}

// freq_map is used by multiple threads, so need to lock it
void update_freq_map(core_t * core, db_t * db) {
    bam_hdr_t * hdr = core->bam_hdr;
    khash_t(freqm) * freq_map = core->freq_map;

    for(int i=0;i<db->n_bam_recs;i++){
        bam1_t *record = db->bam_recs[i];
        int8_t rev = bam_is_rev(record);
        int32_t tid = record->core.tid;
        assert(tid < hdr->n_targets);
        const char *tname = tid >= 0 ? hdr->target_name[tid] : "*";
        uint32_t seq_len = record->core.l_qseq;
        char strand = rev ? '-' : '+';

        modbase_t ** bases = db->modbases[i];

        for(int j=0;j<core->opt.n_mods;j++){
            char mod_code = core->opt.mod_codes_str[j];
            if(req_mods[(int)mod_code]==255){ // mod code not required
                continue;
            }
            for(int seq_i=0;seq_i<seq_len;seq_i++){
                modbase_t * base = &bases[j][seq_i];

                if(base->ref_pos == -1){ // no modification
                    continue;
                }

                uint8_t is_mod = 0, is_called = 0;
                uint8_t thresh = req_threshes[(int)mod_code];
                uint8_t mod_prob = base->mod_prob;
                
                if(mod_prob >= 255-thresh){ // modified with mod_code
                    is_called = 1;
                    is_mod = 1;
                } else if(mod_prob <= thresh){ // not modified with mod_code
                    is_called = 1;
                } else { // ambiguous
                    continue;
                }

                char *key = make_key(tname, base->ref_pos, mod_code, strand);           
                khiter_t k = kh_get(freqm, freq_map, key);
                if (k == kh_end(freq_map)) { // not found, add to map
                    freq_t * freq = (freq_t *)malloc(sizeof(freq_t));
                    MALLOC_CHK(freq);
                    freq->n_called = is_called;
                    freq->n_mod = is_mod;
                    int ret;
                    k = kh_put(freqm, freq_map, key, &ret);
                    kh_value(freq_map, k) = freq;
                } else { // found, update the values
                    free(key);
                    freq_t * freq = kh_value(freq_map, k);
                    freq->n_called += is_called;
                    freq->n_mod += is_mod;
                }
                
            }
        }
    }
}

void print_freq_output(opt_t opt, khash_t(freqm) *freq_map) {
    FILE * output_fp;
    if(opt.output_file == NULL) {
        output_fp = stdout;
    } else {
        output_fp = fopen(opt.output_file, "w");
        if(output_fp == NULL) {
            ERROR("Could not open output file %s", opt.output_file);
            exit(EXIT_FAILURE);
        }
    }
    if(!opt.bedmethyl_out) {
        fprintf(output_fp, "contig\tstart\tend\tstrand\tn_called\tn_mod\tfreq\tmod_code\n");
    }


    if(opt.bedmethyl_out) {
        // chrom, start, end, mod_code, n_called, strand, start, end, "255,0,0",  n_called, freq
        khint_t k;
        for (k = kh_begin(freq_map); k != kh_end(freq_map); ++k) {
            if (kh_exist(freq_map, k)) {
                freq_t* freq = kh_value(freq_map, k);
                double freq_value = (double)freq->n_mod*100/freq->n_called;
                char *contig = NULL;
                int ref_pos;
                char mod_code;
                char strand;
                char * key = (char *) kh_key(freq_map, k);
                decode_key(key, &contig, &ref_pos, &mod_code, &strand);
                int end = ref_pos+1;
                fprintf(output_fp, "%s\t%d\t%d\t%c\t%d\t%c\t%d\t%d\t255,0,0\t%d\t%f\n", contig, ref_pos, end, mod_code, freq->n_called, strand, ref_pos, end, freq->n_called, freq_value);
                free(contig);
            }
        }
        
    } else {
        // contig, start, end, strand, n_called, n_mod, freq, mod_code
        khint_t k;
        for (k = kh_begin(freq_map); k != kh_end(freq_map); ++k) {
            if (kh_exist(freq_map, k)) {
                freq_t* freq = kh_value(freq_map, k);
                double freq_value = (double)freq->n_mod/freq->n_called;
                char * contig = NULL;
                int ref_pos;
                char mod_code;
                char strand;
                char * key = (char *) kh_key(freq_map, k);
                decode_key(key, &contig, &ref_pos, &mod_code, &strand);

                fprintf(output_fp, "%s\t%d\t%d\t%c\t%d\t%d\t%f\t%c\n", contig, ref_pos, ref_pos, strand, freq->n_called, freq->n_mod, freq_value, mod_code);
                free(contig);
            }
        }
    }
    fclose(output_fp);
}

static void get_aln(int ** aln, bam_hdr_t *hdr, bam1_t *record){
    int32_t tid = record->core.tid;
    assert(tid < hdr->n_targets);
    int32_t pos = record->core.pos;
    int32_t end = bam_endpos(record);

    int8_t rev = bam_is_rev(record);

    uint32_t *cigar = bam_get_cigar(record);
    uint32_t n_cigar = record->core.n_cigar;

    int seq_len = record->core.l_qseq;
  
    int read_pos = 0;
    int ref_pos = pos;

    int * aligned_pairs = *aln;
    //fill the aligned_pairs array with -1
    for(int i=0;i<seq_len;i++){
        aligned_pairs[i] = -1;
    }

    for (uint32_t ci = 0; ci < n_cigar; ++ci) {
        uint32_t c = cigar[ci];
        if(rev) {
            c = cigar[n_cigar - ci - 1];
        }
        int cigar_len = bam_cigar_oplen(c);
        int cigar_op = bam_cigar_op(c);

        // Set the amount that the ref/read positions should be incremented
        // based on the cigar operation
        int read_inc = 0;
        int ref_inc = 0;

        // Process match between the read and the reference
        int8_t is_aligned = 0;
        if(cigar_op == BAM_CMATCH || cigar_op == BAM_CEQUAL || cigar_op == BAM_CDIFF) {
            is_aligned = 1;
            read_inc = 1;
            ref_inc = 1;
        } else if(cigar_op == BAM_CDEL) {
            ref_inc = 1;
        } else if(cigar_op == BAM_CREF_SKIP) {
            // end the current segment and start a new one
            //out.push_back(AlignedSegment());
            ref_inc = 1;
        } else if(cigar_op == BAM_CINS) {
            read_inc = 1;
        } else if(cigar_op == BAM_CSOFT_CLIP) {
            read_inc = 1;
        } else if(cigar_op == BAM_CHARD_CLIP) {
            read_inc = 0;
            ERROR("Hard clipping(%d) not supported. Use minimap2 with -Y to use soft clipping for suplimentary alignment.\n", cigar_op);
            exit(EXIT_FAILURE);
        } else {
            ERROR("Unhandled CIGAR OPT Cigar: %d\n", cigar_op);
            exit(EXIT_FAILURE);
        }

        // Iterate over the pairs of aligned bases
        for(int j = 0; j < cigar_len; ++j) {
            if(is_aligned) {
                ASSERT_MSG(read_pos < seq_len, "read_pos:%d seq_len:%d\n", read_pos, seq_len);
                int start = ref_pos;
                if(rev) {
                    start = pos + end - ref_pos - 1;
                }
                aligned_pairs[read_pos] = start;
            }

            // increment
            read_pos += read_inc;
            ref_pos += ref_inc;
        }
    }
}

static void get_bases(core_t * core, db_t *db, int32_t bam_i, const char *mm_string, uint8_t *ml, uint32_t ml_len, int *aln_pairs, bam_hdr_t *hdr, bam1_t *record) {
    const char *qname = bam_get_qname(record);
    int8_t rev = bam_is_rev(record);
    int32_t tid = record->core.tid;
    assert(tid < hdr->n_targets);
    const char *tname = (tid >= 0) ? hdr->target_name[tid] : "*";
    uint8_t *seq = bam_get_seq(record);
    uint32_t seq_len = record->core.l_qseq;

    // 5 int arrays to keep base pos of A, C, G, T, N bases.
    // A: 0, C: 1, G: 2, T: 3, U:4, N: 5
    // so that, nth base of A is at base_pos[0][n] and so on.
    int **bases_pos = db->bases_pos[bam_i];
    int bases_pos_lens[N_BASES] = {0};
    memset(db->mod_codes[bam_i], 0, core->opt.n_mods);

    int i;
    for(i=0;i<seq_len;i++){
        int base_char = seq_nt16_str[bam_seqi(seq, i)];
        int idx = base_idx_lookup[(int)base_char];
        bases_pos[idx][bases_pos_lens[idx]++] = i;

        modbase_t base;
        base.ref_pos = -1;
        base.mod_prob = 0;
        for(int j=0;j<core->opt.n_mods;j++){
            db->modbases[bam_i][j][i] = base;
        }
    }

    int mm_str_len = strlen(mm_string);
    i = 0;
    int ml_start_idx = 0;

    char modbase;
    // char mod_strand;  // commented for now. might need to revisit
    char * mod_codes = db->mod_codes[bam_i];
    int mod_codes_len;
    int * skip_counts = db->skip_counts[bam_i];
    int skip_counts_len;
    // char status_flag;

    while (i < mm_str_len) {
        // reset skip counts and mod codes
        skip_counts_len = 0;
        mod_codes_len = 0;

        // set default status flag to '.' (when not present or '.' in the MM string)
        // status_flag = '.';

        // get base
        if(i < mm_str_len) {
            ASSERT_MSG(valid_bases[(int)mm_string[i]], "Invalid base:%c\n", mm_string[i]);
            modbase = mm_string[i];
            i++;
        }

        // get strand
        if(i < mm_str_len) {
            ASSERT_MSG(valid_strands[(int)mm_string[i]], "Invalid strand:%c\n", mm_string[i]);
            // mod_strand = mm_string[i];
            i++;
        }

        // get base modification codes. can handle multiple codes giver as chars. TO-DO: handle when given as a ChEBI id
        int j = 0;
        while (i < mm_str_len && mm_string[i] != ',' && mm_string[i] != ';' && mm_string[i] != '?' && mm_string[i] != '.') {

            ASSERT_MSG(valid_mod_codes[(int)mm_string[i]], "Invalid base modification code:%c\n", mm_string[i]);

            if(j >= db->mod_codes_cap[bam_i]) {
                db->mod_codes_cap[bam_i] *= 2;
                db->mod_codes[bam_i] = (char *)realloc(db->mod_codes[bam_i], sizeof(char) * (db->mod_codes_cap[bam_i] + 1)); // +1 for null terminator
                MALLOC_CHK(db->mod_codes[bam_i]);
            }

            mod_codes[j] = mm_string[i];
            j++;

            i++;
        }
        mod_codes[j] = '\0';
        mod_codes_len = j;

        ASSERT_MSG(mod_codes_len>0 && mod_codes_len <= 16, "mod_codes_len:%d\n", mod_codes_len);

        // get modification status flag
        if(i < mm_str_len && ( mm_string[i] == '?' || mm_string[i] == '.' )) {
            // status_flag = mm_string[i];
            i++;
        } else { // if not present, set to '.'
            // status_flag = '.';
        }

        // get skip counts
        int k = 0;
        while (i < mm_str_len && mm_string[i] != ';') {

            // skip if a comma
            if(i < mm_str_len && mm_string[i] == ',') {
                i++;
                continue;
            }

            char skip_count_str[10];
            int l = 0;
            while (i < mm_str_len && mm_string[i] != ',' && mm_string[i] != ';') {
                skip_count_str[l] = mm_string[i];
                i++;
                l++;
                assert(l < 10); // if this fails, use dynamic allocation for skip_count_str
            }
            skip_count_str[l] = '\0';
            ASSERT_MSG(l > 0, "invalid skip count:%d.\n", l);
            sscanf(skip_count_str, "%d", &skip_counts[k]);
            ASSERT_MSG(skip_counts[k] >= 0, "skip count cannot be negative: %d.\n", skip_counts[k]);
            
            k++;
        }
        skip_counts_len = k;
        i++;

        if(skip_counts_len == 0) { // no skip counts, no modification
            continue;
        }

        int base_rank = -1;

        int ml_idx = ml_start_idx;
        for(int c=0; c<skip_counts_len; c++) {
            base_rank += skip_counts[c] + 1;
            char mb;
            int idx;
            int read_pos;

            if(rev) {
                mb = base_complement_lookup[(int)modbase];
            } else {
                mb = modbase;
            }

            idx = base_idx_lookup[(int)mb];

            // print_array(bases_pos_lens, 5, 'i');
            if(base_rank >= bases_pos_lens[idx]) {
                WARNING("%d th base of %c not found in SEQ. %c base count is %d read_id:%s seq_len:%d mod.base:%c mod_codes:%s\n", base_rank, mb, mb, bases_pos_lens[idx], qname, seq_len, modbase, mod_codes);
                continue;
            }
            ASSERT_MSG(base_rank < bases_pos_lens[idx], "%d th base of %c not found in SEQ. %c base count is %d read_id:%s seq_len:%d mod.base:%c mod_codes:%s\n", base_rank, mb, mb, bases_pos_lens[idx], qname, seq_len, modbase, mod_codes);
            
            if(rev) {
                read_pos = seq_len - bases_pos[idx][bases_pos_lens[idx] - base_rank - 1] - 1;
            } else {
                read_pos = bases_pos[idx][base_rank];
            }

            ASSERT_MSG(read_pos>=0 && read_pos < seq_len, "Base pos cannot exceed seq len. read_pos: %d seq_len: %d\n", read_pos, seq_len);

            int ref_pos = aln_pairs[read_pos];
            if(ref_pos == -1) {
                if(mod_codes_len > 0) {
                    ml_idx = ml_start_idx + c*mod_codes_len + mod_codes_len - 1;
                }
                continue;
            }

            int8_t is_cpg = 0;
            ref_t *ref = get_ref(tname);
            ASSERT_MSG(ref != NULL, "Contig %s not found in ref_map\n", tname);
            ASSERT_MSG(ref_pos >= 0 && ref_pos < ref->ref_seq_length, "ref_pos:%d ref_len:%d\n", ref_pos, ref->ref_seq_length);
            ASSERT_MSG(ref->ref_seq_length == hdr->target_len[tid], "ref_len:%d target_len:%d\n", ref->ref_seq_length, hdr->target_len[tid]);
            char * ref_seq = ref->forward;
            if ((!rev && ref_pos + 1 < ref->ref_seq_length && ref_seq[ref_pos] == 'C' && ref_seq[ref_pos + 1] == 'G') ||
                (rev && ref_pos > 0 && ref_seq[ref_pos] == 'G' && ref_seq[ref_pos - 1] == 'C')) {
                is_cpg = 1;
            }
            if(is_cpg == 0) {
                if(mod_codes_len > 0) {
                    ml_idx = ml_start_idx + c*mod_codes_len + mod_codes_len - 1;
                }
                continue;
            }
            
            // mod prob per each mod code. TO-DO: need to change when code is ChEBI id
            for(int m=0; m<mod_codes_len; m++) {
                char mod_code = mod_codes[m];
                
                ml_idx = ml_start_idx + c*mod_codes_len + m;

                if(req_mods[(int)mod_code]==255){ // mod code not required
                    continue;
                }

                ASSERT_MSG(ml_idx<ml_len, "ml_idx:%d ml_len:%d\n", ml_idx, ml_len);
                uint8_t mod_prob = ml[ml_idx];
                ASSERT_MSG(mod_prob <= 255 && mod_prob>=0, "mod_prob:%d\n", mod_prob);

                db->modbases[bam_i][req_mods[(int)mod_code]][read_pos].mod_prob = mod_prob;
                db->modbases[bam_i][req_mods[(int)mod_code]][read_pos].ref_pos = ref_pos;

            }

        }
        ml_start_idx = ml_idx + 1;

    }
}

// void print_view_output(core_t* core, db_t* db) {
//     fprintf(core->opt.output_fp, "ref_contig\tref_pos\tstrand\tread_id\tread_pos\tmod_code\tmod_prob\n");
//     for(int i=0;i<db->n_bam_recs;i++){
//         bam1_t *record = db->bam_recs[i];
//         bam_hdr_t * hdr = core->bam_hdr;
//         const char *qname = bam_get_qname(record);
//         int32_t tid = record->core.tid;
//         assert(tid < hdr->n_targets);
//         const char *tname = tid >= 0 ? hdr->target_name[tid] : "*";
//         int32_t seq_len = record->core.l_qseq;

//         int8_t rev = bam_is_rev(record);
//         char strand = rev ? '-' : '+';

//         modbase_t ** bases = db->modbases[i];
//         for(int j=0;j<core->opt.n_mods;j++){
//             char mod_code = core->opt.mod_codes_str[j];
//             if(req_mods[(int)mod_code]==255){ // mod code not required
//                 continue;
//             }
//             for(int seq_i=0;seq_i<seq_len;seq_i++){
//                 modbase_t base = bases[j][seq_i];
            
//                 if(base.ref_pos == -1){
//                     continue;
//                 }

//                 fprintf(core->opt.output_fp, "%s\t%d\t%c\t%s\t%d\t%c\t%f\n", tname, base.ref_pos, strand, qname, seq_i, mod_code, base.mod_prob/255.0);
            
//             }
//         }
//     }

// }

void dump_stats_map(const char * dump_file, khash_t(freqm) * freq_map){

    FILE *fp = fopen(dump_file, "wb");
    if (fp == NULL) {
        ERROR("Cannot open dump file %s\n", dump_file);
        exit(EXIT_FAILURE);
    }

    for (khiter_t k = kh_begin(freq_map); k != kh_end(freq_map); ++k) {
        if (kh_exist(freq_map, k)) {
            freq_t * freq = kh_value(freq_map, k);
            char *contig = NULL;
            int start;
            char mod_code;
            char strand;
            char * key = (char *) kh_key(freq_map, k);
            decode_key(key, &contig, &start, &mod_code, &strand);
            size_t chrom_len = strlen(contig);
            // fprintf(stderr, "writing chrom_len:%ld chrom:%s start:%d end:%d depth:%d n_mod:%d n_called:%d n_skipped:%d mod_code:%c mod_strand:%c ref_base:%c is_aln_cpg:%d\n", chrom_len, freq->contig, freq->start, freq->end, freq->depth, freq->n_mod, freq->n_called, freq->n_skipped, freq->mod_code, freq->strand, freq->ref_base, freq->is_aln_cpg);

            size_t r = fwrite(&chrom_len, sizeof(size_t), 1, fp);
            ASSERT_MSG(r==1, "Error writing chrom_len to dump file: %s r:%ld\n", dump_file, r);
            r = fwrite(contig, sizeof(char), chrom_len, fp);
            ASSERT_MSG(r==chrom_len, "Error writing chrom to dump file: %s r:%ld\n", dump_file, r);
            r = fwrite(&start, sizeof(int), 1, fp);
            ASSERT_MSG(r==1, "Error writing start to dump file: %s r:%ld\n", dump_file, r);
            r = fwrite(&freq->n_mod, sizeof(uint16_t), 1, fp);
            ASSERT_MSG(r==1, "Error writing n_mod to dump file: %s r:%ld\n", dump_file, r);
            r = fwrite(&freq->n_called, sizeof(uint16_t), 1, fp);
            ASSERT_MSG(r==1, "Error writing n_called to dump file: %s r:%ld\n", dump_file, r);
            r = fwrite(&mod_code, sizeof(char), 1, fp);
            ASSERT_MSG(r==1, "Error writing mod_code to dump file: %s r:%ld\n", dump_file, r);
            r = fwrite(&strand, sizeof(char), 1, fp);
            ASSERT_MSG(r==1, "Error writing mod_strand to dump file: %s r:%ld\n", dump_file, r);

            free(contig);
        }
    }

    fclose(fp);
}

void load_stats_map(const char * dump_file, khash_t(freqm) * freq_map){
    FILE *fp = fopen(dump_file, "rb");
    if (fp == NULL) {
        ERROR("Cannot open dump file %s\n", dump_file);
        exit(EXIT_FAILURE);
    }

    while (1) {
        
        size_t chrom_len;
        size_t r = fread(&chrom_len, sizeof(size_t), 1, fp);

        if(feof(fp)){
            fprintf(stderr, "[realfreq] reading dump file %s reached EOF\n", dump_file);
            break;
        }

        ASSERT_MSG(r==1, "Error reading chrom_len from dump file: %s r:%ld\n", dump_file, r);        

        freq_t * freq = (freq_t *)malloc(sizeof(freq_t));
        MALLOC_CHK(freq);
        char * contig = (char *)malloc((chrom_len+1)*sizeof(char));
        MALLOC_CHK(contig);
        int start;
        char mod_code;
        char strand;

        contig[chrom_len] = '\0';

        r = fread(contig, sizeof(char), chrom_len, fp);
        ASSERT_MSG(r==chrom_len, "Error reading chrom from dump file: %s r:%ld\n", dump_file, r);
        r = fread(&start, sizeof(int), 1, fp);
        ASSERT_MSG(r==1, "Error reading start from dump file: %s\n r:%ld", dump_file, r);
        r = fread(&freq->n_mod, sizeof(uint16_t), 1, fp);
        ASSERT_MSG(r==1, "Error reading n_mod from dump file: %s\n r:%ld", dump_file, r);
        r = fread(&freq->n_called, sizeof(uint16_t), 1, fp);
        ASSERT_MSG(r==1, "Error reading n_called from dump file: %s\n r:%ld", dump_file, r);
        r = fread(&mod_code, sizeof(char), 1, fp);
        ASSERT_MSG(r==1, "Error reading mod_code from dump file: %s\n r:%ld", dump_file, r);
        r = fread(&strand, sizeof(char), 1, fp);
        ASSERT_MSG(r==1, "Error reading mod_strand from dump file: %s\n r:%ld", dump_file, r);

        char *key = make_key(contig, start, mod_code, strand);
        int ret;
        khiter_t k = kh_put(freqm, freq_map, key, &ret);
        kh_value(freq_map, k) = freq;

        free(contig);
    }

    fclose(fp);
}




void modbases_single(core_t* core, db_t* db, int32_t i) {
    bam1_t *record = db->bam_recs[i];

    const char *mm = db->mm[i];
    uint32_t ml_len = db->ml_lens[i];
    uint8_t *ml = db->ml[i];

    bam_hdr_t *hdr = core->bam_hdr;

    get_aln(&(db->aln[i]), hdr, record);
    get_bases(core, db, i, mm, ml, ml_len, db->aln[i], hdr, record);
}
