#include <stdio.h>
#include <getopt.h>
#include "file_handler.h"
#include "mod.h"

int main(int argc, char* argv[]) {
    
    const char* output_file;

    //parse the user args
    const char* optstring = "o:";
    int opt;
    while ((opt = getopt(argc, argv, optstring)) != -1) {
        switch (opt) {
            case 'o':
                output_file = optarg;
                break;
            default:
                fprintf(stderr, "Usage: %s [-o output_file]\n", argv[0]);
                exit(EXIT_FAILURE);
        }
    }

    init_maps();
    read_files_from_stdin(output_file);
    destroy_maps();
    return 0;
}
