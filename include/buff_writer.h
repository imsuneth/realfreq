
#ifndef BUFF_WRITER_H
#define BUFF_WRITER_H

#include <stdio.h>

typedef struct {
    char * buffer;
    size_t size;
    size_t capacity;
    char * filename;
    FILE * fp;
} buff_writer_t;

buff_writer_t * buff_writer_init(size_t capacity, char * filename);
void buff_writer_pushf(buff_writer_t * writer, const char * format, ...);
void buff_writer_destroy(buff_writer_t * writer);
void buff_writer_flush(buff_writer_t * writer);
void buff_writer_clear_file(buff_writer_t * writer);

#endif