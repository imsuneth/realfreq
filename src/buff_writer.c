/**
 * @file buff_write.c
 * @brief write to a file using a buffer

MIT License

Copyright (c) 2024 Suneth Samarasinghe

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

#include "buff_writer.h"
#include "error.h"
#include <stdlib.h>
#include <stdarg.h>

FILE * fp;

buff_writer_t * buff_writer_init(size_t capacity, char * filename) {

    buff_writer_t * writer = (buff_writer_t *)malloc(sizeof(buff_writer_t));
    MALLOC_CHK(writer);

    writer->filename = filename;

    if(writer->filename == NULL) {
        fp = stdout;
    } else {
        writer->fp = fopen(filename, "w");
        if(writer->fp == NULL) {
            ERROR("Cannot open file %s", filename);
        }
    }

    writer->buffer = (char *)malloc(sizeof(char)*capacity);
    MALLOC_CHK(writer->buffer);

    writer->size = 0;
    writer->capacity = capacity;
    return writer;
}

void buff_writer_clear_file(buff_writer_t * writer) {
    if(writer->fp != NULL) {
        fclose(writer->fp);
    }
    writer->fp = fopen(writer->filename, "w");
    if(writer->fp == NULL) {
        ERROR("Cannot open file %s", writer->filename);
    }
}

void buff_writer_pushf(buff_writer_t * writer, const char * format, ...) {
    va_list args;
    va_start(args, format);
    size_t len = vsnprintf(NULL, 0, format, args);
    va_end(args);

    if(writer->size + len >= writer->capacity) { // if the buffer is full, write to file
        fwrite(writer->buffer, sizeof(char), writer->size, writer->fp);
        writer->size = 0;
    }

    va_start(args, format);
    vsnprintf(writer->buffer + writer->size, len + 1, format, args);
    va_end(args);

    writer->size += len;
}

void buff_writer_flush(buff_writer_t * writer) {
    if(writer->size > 0) {
        fwrite(writer->buffer, sizeof(char), writer->size, writer->fp);
        writer->size = 0;
    }
}

void buff_writer_destroy(buff_writer_t * writer) {
    fwrite(writer->buffer, sizeof(char), writer->size, writer->fp);
    free(writer->buffer);
    free(writer);
    if(writer->filename != NULL && writer->fp != NULL) {
        fclose(writer->fp);
    }
}