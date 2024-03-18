#include <stdio.h>
#include <getopt.h>
#include <stdbool.h>
#include "file_handler.h"
#include "logger.h"
#include "mod.h"

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
                set_processed_files_log(optarg);
                break;
            case 'y':
                clear_log(true);
                break;
            default:
                fprintf(stderr, "Usage: %s [-o output_file] [-l processed_files_log]\n", argv[0]);
                exit(EXIT_FAILURE);
        }
    }

    init_maps();
    read_files_from_stdin();
    destroy_maps();
    return 0;
}
