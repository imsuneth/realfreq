#ifndef LOGGER_H
#define LOGGER_H

#include <stdbool.h>

void set_processed_files_log(const char *filepath);
void lear_log(bool clear);
void log_file_processed(const char *bamfile, double realtime_meth_freq, double realtime_write_output);

#endif /* LOGGER_H */
