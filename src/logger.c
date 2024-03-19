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

