#include <stdio.h>
#include "logger.h"

void log_file_processed(const char *path) {
    FILE *log_file = fopen(LOG_FILE_PATH, "a");
    if (log_file) {
        fprintf(log_file, "Processed file: %s\n", path);
        fclose(log_file);
    }
}

