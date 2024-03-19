#include <stdio.h>
#include <getopt.h>
#include <stdbool.h>
#include <stdlib.h>
#include "file_handler.h"
#include "logger.h"
#include "mod.h"
#include "error.h"

static const char * logfile = "realfreq_processed.log";

void initialize() {
    init_logger(logfile);
    mod_init();
}

void destroy() {
    destroy_mod();
    destroy_logger();
}

int main(int argc, char* argv[]) {
    //parse the user args
    const char* optstring = "yo:l:";
    int opt;
    while ((opt = getopt(argc, argv, optstring)) != -1) {
        switch (opt) {
            case 'o':
                set_output_file(optarg);
                break;
            case 'l':
                logfile = (const char *)malloc(sizeof(char)*strlen(optarg));
                MALLOC_CHK(logfile);
                strcpy(logfile, optarg);
                break;
            case 'y':
                clear_log(true);
                break;
            default:
                fprintf(stderr, "Usage: %s [-o output_file] [-l processed_files_log]\n", argv[0]);
                exit(EXIT_FAILURE);
        }
    }

    initialize();
    read_files_from_stdin();
    destroy();
    return 0;
}
