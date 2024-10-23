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
#include <pthread.h>
#include "mod.h"
#include "error.h"
#include "ref.h"
#include "misc.h"
#include "server.h"

#define FILEPATH_LEN 1024
pthread_t server_thread;

void print_help_msg(FILE * fp_help, opt_t opt) {
    fprintf(fp_help, "Usage: realfreq [options..] ref.fa\nrealfreq reads the input bam file path from stdin\n");
    fprintf(fp_help, "Options:\n");
    fprintf(fp_help,"   -b                         output in bedMethyl format [%s]\n", (opt.bedmethyl_out?"yes":"not set"));                                                    //0
    fprintf(fp_help,"   -c STR                     modification codes (ex. m , h or mh) [%s]\n", opt.mod_codes_str);                                                            //1
    fprintf(fp_help,"   -m FLOAT                   min modification threshold(s). Comma separated values for each modification code given in -c [%s]\n", opt.mod_threshes_str); //2
    fprintf(fp_help,"   -t INT                     number of processing threads [%d]\n",opt.num_thread);                                                                        //3
    fprintf(fp_help,"   -K INT                     batch size (max number of reads loaded at once) [%d]\n",opt.batch_size);                                                     //4
    fprintf(fp_help,"   -B FLOAT[K/M/G]            max number of bytes loaded at once [%.1fM]\n",opt.batch_size_bytes/(float)(1000*1000));                                      //5
    fprintf(fp_help,"   -h                         help\n");                                                                                                                    //6
    fprintf(fp_help,"   -p INT                     print progress every INT seconds (0: per batch) [%d]\n", opt.progress_interval);                                             //7
    fprintf(fp_help,"   -o FILE                    output file [%s]\n", opt.output_file==NULL?"stdout":opt.output_file);                                                        //8
    fprintf(fp_help,"   -d FILE                    dump file [%s]\n", opt.dump_file);                                                                                           //9
    fprintf(fp_help,"   -l FILE                    progress log file [%s]\n", opt.log_file);                                                                                    //10
    fprintf(fp_help,"   -r                         resume from dump file [%s]\n", (opt.is_resuming?"yes":"no"));                                                                //11
    fprintf(fp_help,"   -s PORT                    start server on PORT [%d]\n", opt.server_port);                                                                              //12
    fprintf(fp_help,"   -v INT                     verbosity level [%d]\n",(int)get_log_level());                                                                               //13
    fprintf(fp_help,"   -V                         print version\n");                                                                                                           //14
}

void start_realfreq(opt_t opt, khash_t(freqm)* freq_map) {
    char *bam_file = (char *)malloc(FILEPATH_LEN * sizeof(char));
    INFO("%s", "reading file path from stdin");

    while (1) {
        double realtime_prog = realtime();

        int ret = fscanf(stdin, "%s", bam_file);
        if (ferror(stdin)) {
            INFO("%s", "error reading from stdin");
            print_freq_output(opt, freq_map);
            break;
        }
        if (feof(stdin)) {
            INFO("%s", "end of stdin");
            print_freq_output(opt, freq_map);
            break;
        }
        if(ret!=1) {
            INFO("%s", "error reading from stdin");
            print_freq_output(opt, freq_map);
            break;
        }

        double realtime0 = realtime();
        if (strstr(bam_file, ".bam") != NULL) { //bam file
            INFO("processing file %s", bam_file);
            //initialise the core data structure
            core_t* core = init_core(bam_file, opt, realtime0);
            core->freq_map = freq_map;

            int32_t counter=0;

            //initialise a databatch
            db_t* db = init_db(core);

            ret_status_t status = {core->opt.batch_size,core->opt.batch_size_bytes};
            while (status.num_reads >= core->opt.batch_size || status.num_bytes>=core->opt.batch_size_bytes) {

                //load a databatch
                status = load_db(core, db);

                //process the data batch
                process_db(core, db);

                //write the output
                output_db(core, db);

                free_db_tmp(core, db);

                //print progress
                if(opt.progress_interval<=0 || realtime()-realtime_prog > opt.progress_interval){
                    fprintf(stderr, "[%s::%.3f*%.2f] %d Entries (%.1fM bytes) processed\t%d Entries (%.1fM bytes) skipped\n", __func__,
                            realtime() - realtime0, cputime() / (realtime() - realtime0),
                            (db->n_bam_recs), (db->sum_bytes)/(1000.0*1000.0),
                            db->skipped_reads,db->skipped_reads_bytes/(1000.0*1000.0));
                    realtime_prog = realtime();
                }

                //check if 90% of total reads are skipped
                if(core->skipped_reads>0.9*core->total_reads){
                    WARNING("%s","90% of the reads are skipped. Possible causes: unmapped bam, zero sequence lengths, or missing MM, ML tags (not performed base modification aware basecalling). Refer https://github.com/warp9seq/minimod for more information.");
                }

                if(opt.debug_break==counter){
                    break;
                }
                counter++;
            }

            // free the databatch
            free_db(core, db);

            fprintf(stderr, "[%s] total entries: %ld", __func__,(long)core->total_reads);
            fprintf(stderr,"\n[%s] total bytes: %.1f M",__func__,core->sum_bytes/(float)(1000*1000));
            fprintf(stderr,"\n[%s] total skipped entries: %ld",__func__,(long)core->skipped_reads);
            fprintf(stderr,"\n[%s] total skipped bytes: %.1f M",__func__,core->skipped_reads_bytes/(float)(1000*1000));
            fprintf(stderr,"\n[%s] total processed entries: %ld",__func__,(long)(core->total_reads-core->skipped_reads));
            fprintf(stderr,"\n[%s] total processed bytes: %.1f M",__func__,(core->sum_bytes-core->skipped_reads_bytes)/(float)(1000*1000));

            fprintf(stderr, "\n[%s] Data loading time: %.3f sec", __func__,core->load_db_time);
            fprintf(stderr, "\n[%s] Data processing time: %.3f sec", __func__,core->process_db_time);
            if((core->opt.flag&MINIMOD_PRF)|| core->opt.flag & MINIMOD_ACC){
                    fprintf(stderr, "\n[%s]     - Parse time: %.3f sec",__func__, core->parse_time);
                    fprintf(stderr, "\n[%s]     - Calc time: %.3f sec",__func__, core->calc_time);
            }
            fprintf(stderr, "\n[%s] Data output time: %.3f sec", __func__,core->output_time);
            fprintf(stderr,"\n");

            //free the core data structure
        free_core(core,opt);

        } else if (strstr(bam_file, ".tsv") != NULL) { //tsv file
            process_tsv_file(bam_file, opt, freq_map);
            print_freq_output(opt, freq_map);
            dump_stats_map(opt.dump_file, freq_map);
            INFO("processed file %s", bam_file);
            continue;
        }
        

        double realtime1 = realtime();

        print_freq_output(opt, freq_map);
        dump_stats_map(opt.dump_file, freq_map);

        double realtime2 = realtime();

        //write to log file
        if(opt.log_file!=NULL){
            FILE *fp = fopen(opt.log_file, "a");
            if (fp == NULL) {
                ERROR("Cannot open log file %s for writing", opt.log_file);
                exit(EXIT_FAILURE);
            }
            fprintf(fp, "%s\t%.3f sec\t%.3f sec\t%d entries\n", bam_file, realtime1 - realtime0, realtime2 - realtime1, kh_size(freq_map));
            fclose(fp);
        }
        
        
    }
    free(bam_file);
}

int main(int argc, char* argv[]) {

    double realtime0 = realtime();

    char *ref_file = NULL;

    FILE *fp_help = stderr;

    opt_t opt;
    init_opt(&opt); //initialise options to defaults
    opt.subtool = MOD_FREQ;
    //parse the user args
    const char* optstring = "bc:m:t:K:B:hp:o:d:l:rs:v:V";
    int32_t c = -1;
    //parse the user args
    while ((c = getopt(argc, argv, optstring)) != -1) {
        if (c == 'b') {
            opt.bedmethyl_out = 1;
        } else if (c == 'c') {
            opt.mod_codes_str = optarg;
        } else if (c == 'm') {
            opt.mod_threshes_str = optarg;
        } else if (c == 't') {
            opt.num_thread = atoi(optarg);
            if (opt.num_thread < 1) {
                ERROR("Number of threads should larger than 0. You entered %d", opt.num_thread);
                exit(EXIT_FAILURE);
            }
        } else if (c == 'K') {
            opt.batch_size = atoi(optarg);
            if (opt.batch_size < 1) {
                ERROR("Batch size should larger than 0. You entered %d",opt.batch_size);
                exit(EXIT_FAILURE);
            }
        } else if (c == 'B') {
            opt.batch_size_bytes = mm_parse_num(optarg);
            if(opt.batch_size_bytes<=0){
                ERROR("%s","Maximum number of bytes should be larger than 0.");
                exit(EXIT_FAILURE);
            }
        } else if (c == 'h') {
            print_help_msg(fp_help, opt);
            exit(EXIT_SUCCESS);
        } else if (c == 'p') {
            if (atoi(optarg) < 0) {
                ERROR("Progress interval should be 0 or positive. You entered %d", atoi(optarg));
                exit(EXIT_FAILURE);
            }
            opt.progress_interval = atoi(optarg);
        } else if (c == 'o') {
            FILE *fp = fopen(optarg, "w");
            if (fp == NULL) {
                ERROR("Cannot open output file %s for writing", optarg);
                exit(EXIT_FAILURE);
            }
            fclose(fp);
            opt.output_file = optarg;
        } else if (c == 'd') {
            opt.dump_file = optarg;
        } else if (c == 'l') {
            FILE *fp = fopen(optarg, "a");
            if (fp == NULL) {
                ERROR("Cannot open log file %s for writing", optarg);
                exit(EXIT_FAILURE);
            }
            fclose(fp);
            opt.log_file = optarg;
        } else if (c == 'r') {
            opt.is_resuming = 1;
        } else if (c == 's') {
            opt.server_port = atoi(optarg);
        } else if (c == 'v') {
            int v = atoi(optarg);
            set_log_level((enum log_level_opt)v);
        } else if (c == 'V') {
            fprintf(stderr, "realfreq version %s\n", VERSION);
            exit(EXIT_SUCCESS);
        } else {
            ERROR("Unknown option %c", c);
            exit(EXIT_FAILURE);
        }
    }

    if(opt.mod_codes_str==NULL){
        INFO("%s", "Modification codes not provided. Using default modification code m");
        opt.mod_codes_str = "m";
    }
    opt.n_mods = parse_mod_codes(opt.mod_codes_str);
    
    if(opt.mod_threshes_str==NULL){
        INFO("%s", "Modification threshold not provided. Using default threshold 0.8");
    } else {
        parse_mod_threshes(opt.mod_codes_str, opt.mod_threshes_str, opt.n_mods);
    }

    print_mod_options(opt);

    // No arguments given
    if (argc - optind != 1 || fp_help == stdout) {
        WARNING("%s","Missing arguments");
        print_help_msg(fp_help, opt);
        if(fp_help == stdout){
            exit(EXIT_SUCCESS);
        }
        exit(EXIT_FAILURE);
    }

    ref_file = argv[optind];

    if (ref_file == NULL) {
        WARNING("%s","Reference file not provided");
        print_help_msg(fp_help, opt);
        if(fp_help == stdout){
            exit(EXIT_SUCCESS);
        }
        exit(EXIT_FAILURE);
    }

    if (opt.output_file == NULL) {
        ERROR("%s", "Output file is not specified. Use -o or --output flag\n");
        exit(EXIT_FAILURE);
    }

    if (opt.dump_file == NULL) {
        ERROR("%s", "Dump file is not specified. Use -d or --dump flag\n");
        exit(EXIT_FAILURE);
    }

    //load the reference genome
    fprintf(stderr, "[%s] Loading reference genome %s\n", __func__, ref_file);
    load_ref(ref_file);
    fprintf(stderr, "[%s] Reference genome loaded in %.3f sec\n", __func__, realtime()-realtime0);

    //initialise the freq_map
    khash_t(freqm)* freq_map = kh_init(freqm);

    if(opt.is_resuming){
        load_stats_map(opt.dump_file, freq_map);
    }

    // start the server
    server_args_t server_args = {opt.server_port, freq_map};
    if(opt.server_port!=-1) {
        INFO("Starting server on port %d\n", opt.server_port);
        pthread_create(&server_thread, NULL, start_server, (void *) &server_args);
    }
    
    start_realfreq(opt, freq_map);

    if(opt.server_port!=-1) {
        pthread_cancel(server_thread);
        INFO("Server on port %d stopped\n", opt.server_port);
        pthread_join(server_thread, NULL);
    }

    // destroy the freq_map
    khint_t k;
    for (k = kh_begin(freq_map); k != kh_end(freq_map); k++) {
        if (kh_exist(freq_map, k)) {
            char *key = (char *) kh_key(freq_map, k);
            freq_t* freq = kh_value(freq_map, k);
            free(freq);
            free(key);
        }
    }
    kh_destroy(freqm, freq_map);

    destroy_ref();

    INFO("%s", "exiting\n");
    return 0;
}
