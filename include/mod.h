#ifndef MOD_H
#define MOD_H

#include "utils.h"
#include "khash.h"

typedef struct {
    const char * chrom;
    int start;
    int end;
    int depth;
    int n_mod;
    int n_called;
    int n_skipped;
    double freq;
    char mod_code;
    char mod_strand;
    char strand;
} stat_t;

KHASH_MAP_INIT_STR(str, stat_t);

uint16_t *get_meth_tag(bam1_t *record, char *tag, uint32_t *len_ptr);
void meth_freq(core_t* core);

#endif // MOD_H