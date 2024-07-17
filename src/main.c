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
#include "logger.h"
#include "meth.h"
#include "error.h"
#include "ref.h"
#include "misc.h"
#include "server.h"

#define FILEPATH_LEN 1024
static int server_port = -1;
pthread_t server_thread;

void destroy() {
    destroy_meth();
}

void print_usage(FILE * fp) {
    fprintf(fp, "Usage: realfreq [options] -r FILE -o FILE\nrealfreq reads the input bam file path from stdin\n");
    fprintf(fp, "Options:\n");
    fprintf(fp, "  -h, --help                   print this help message\n");
    fprintf(fp, "  -r FILE, --reference FILE    reference file\n");
    fprintf(fp, "  -o FILE, --output FILE       methylation frequency output file\n");
    fprintf(fp, "  -m FLOAT, --mod-thresh FLOAT min modification threshold inclusive [0.2]\n");
    fprintf(fp, "  -b, --bedmethyl              output in bedMethyl format\n");
    fprintf(fp, "  -d FILE, --dump FILE         dump file\n");
    fprintf(fp, "  -s, --resume                 resume from a dump file\n");
    fprintf(fp, "  -c, --server PORT            start the server on PORT\n");
    
}

void start_realfreq(opt_t opt) {
    char *filepath = (char *)malloc(FILEPATH_LEN * sizeof(char));
    INFO("%s", "reading file path from stdin\n");
    while (1) {
        int status = fscanf(stdin, "%s", filepath);
        if (ferror(stdin)) {
            INFO("%s", "error reading from stdin\n");
            break;
        }
        if (feof(stdin)) {
            INFO("%s", "end of stdin\n");
            break;
        }
        if(status!=1) {
            INFO("%s", "error reading from stdin\n");
            break;
        }
        
        INFO("processing file %s\n", filepath);
        double realtime0 = realtime();
    
        //initialise the core data structure
        core_t* core = init_core(filepath, opt, realtime0);

        meth_freq(core);

        double realtime1 = realtime();

        dump_stats_map(core->opt.dump_file);
        write_output(core->opt.out_file, core->opt.bedmethyl_out, opt.mod_code);

        double realtime2 = realtime();

        static log_entry_t log_entry;
        log_entry.bamfile = filepath;
        log_entry.realtime_meth_freq = realtime1 - realtime0;
        log_entry.realtime_write_output = realtime2 - realtime1;
        log_entry.stats_len = get_stats_len();
        log_file_processed(&log_entry);

        //free the core data structure
        free_core(core,opt);
        
    }
    free(filepath);
}

int main(int argc, char* argv[]) {
    opt_t opt;
    init_opt(&opt); //initialise options to defaults
    //parse the user args
    const char* optstring = "hr:o:bd:sc:";
    int o;
    while ((o = getopt(argc, argv, optstring)) != -1) {
        switch (o) {
            case 'h':
                print_usage(stdout);
                exit(EXIT_SUCCESS);
                break;
            case 'r':
                opt.ref_file = optarg;
                break;
            case 'o':
                opt.out_file = optarg;
                break;
            case 'b':
                opt.bedmethyl_out = 1;
                break;
            case 'd':
                opt.dump_file = optarg;
                break;
            case 's':
                opt.is_resuming = 1;
                break;
            case 'c':
                server_port = atoi(optarg);
                break;
            default:
                print_usage(stderr);
                exit(EXIT_FAILURE);
                
        }
    }

    if (opt.ref_file == NULL) {
        ERROR("%s", "Reference file is not specified. Use -r or --reference flag\n");
        exit(EXIT_FAILURE);
    }

    if (opt.out_file == NULL) {
        ERROR("%s", "Output file is not specified. Use -o or --output flag\n");
        exit(EXIT_FAILURE);
    }

    if (opt.dump_file == NULL) {
        ERROR("%s", "Dump file is not specified. Use -d or --dump flag\n");
        exit(EXIT_FAILURE);
    }

    init_meth(opt);

    fprintf(stderr, "[realfreq] loading reference from %s\n", opt.ref_file);
    load_ref(opt.ref_file);

    if(server_port!=-1) {
        pthread_create(&server_thread, NULL, start_server, (void *) &server_port);
    }
    
    start_realfreq(opt);

    if(server_port!=-1) {
        pthread_join(server_thread, NULL);
    }
    
    destroy();

    INFO("%s", "exiting\n");
    return 0;
}
