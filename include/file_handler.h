#ifndef FILE_HANDLER_H
#define FILE_HANDLER_H

#include <stdio.h>

#define FILEPATH_LEN 500

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <htslib/hts.h>
#include <htslib/sam.h>

void read_files_from_stdin(const char *output_file);

#endif /* FILE_HANDLER_H */
