/**
 * @file meth.c
 * @brief methylation calling and frequency calculation

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

#include "error.h"
#include "mod.h"
#include "minimod.h"
#include "server.h"
#include "misc.h"
#include <stdio.h> 
#include <netdb.h> 
#include <netinet/in.h> 
#include <stdlib.h> 
#include <string.h> 
#include <sys/socket.h> 
#include <sys/types.h> 
#include <unistd.h> // read(), write(), close()
#include <pthread.h>
#define MAX 80
#define SA struct sockaddr

char* help = "Available query commands:\n \
            help\n \
            \tshow this help message\n \
            get_contig:<contig>\n \
            \tquery by contig\n \
            get_range:<start_pos>:<end_pos>\n \
            \tquery data between start and end positions (both inclusive)\n \
            get_contig_range:<contig>:<start_pos>:<end_pos>\n \
            \tquery by contig and between start and end positions\n \
            get_contig_range_mod:<contig>:<start_pos>:<end_pos>:<mod_code>\n \
            \tquery by contig and between start and end positions and by mod code\n";

khash_t(freqm)* map;
/* 
test commands:
    nc localhost 8080 <<< help
    nc localhost 8080 <<< get_contig:chr1
    nc localhost 8080 <<< get_range:1:100
    nc localhost 8080 <<< get_contig_range:chr22:18850302:49514860
    nc localhost 8080 <<< get_contig_range_mod:chr22:18850302:49514860:m
*/

void * handle_request(void * clientfd) 
{ 
    char buff[MAX]; 
    memset(buff, 0, MAX);
    int client_fd = *((int *)clientfd);

    // read the message from client
    if (read(client_fd, buff, sizeof(buff)) == 0) {
        INFO("%s", "client disconnected\n");
        pthread_exit(0);
    }

    char* tok;
    char* request = strtok(buff, "\r\n");
    char* response = "Invalid query\n";

    printf("query: %s\n", request); 

    if(strcmp(request, "help") == 0) {
        response = help;
    } else if (request == NULL) {
        response = "Invalid query\n";
    } else {
        tok = strtok(request, ":");
        if(tok == NULL) {
            response = "Invalid query\n";
        } else if(strcmp(tok, "get_contig_range") == 0) {
            char * req_contig = strtok(NULL, ":");
            char * req_start_pos = strtok(NULL, ":");
            char * req_end_pos = strtok(NULL, ":");

            if(req_contig == NULL || req_start_pos == NULL || req_end_pos == NULL) {
                response = "Invalid query. get_contig_range:<contig>:<start_pos>:<end_pos>\n";
            } else {
                response = get_stats_contig_range(req_contig, atoi(req_start_pos), atoi(req_end_pos), map);
                if (response == NULL) {
                    response = "No matching data found\n";
                }
            }

        } else if(strcmp(tok, "get_range") == 0) {
            char * req_start_pos = strtok(NULL, ":");
            char * req_end_pos = strtok(NULL, ":");

            if(req_start_pos == NULL || req_end_pos == NULL) {
                response = "Invalid query. get_range:<start_pos>:<end_pos>\n";
            } else {
                response = get_stats_range(atoi(req_start_pos), atoi(req_end_pos), map);
                if (response == NULL) {
                    response = "No matching data found\n";
                }
            }

        } else if(strcmp(tok, "get_contig") == 0) {
            char * req_contig = strtok(NULL, ":");

            if(req_contig == NULL) {
                response = "Invalid query get_contig:<contig>\n";
            } else {
                response = get_stats_contig(req_contig, map);
                if (response == NULL) {
                    response = "No matching data found\n";
                }
            }

        } else if (strcmp(tok, "get_contig_range_mod") == 0) {
            char * req_contig = strtok(NULL, ":");
            char * req_start_pos = strtok(NULL, ":");
            char * req_end_pos = strtok(NULL, ":");
            char * req_mod_code = strtok(NULL, ":");

            if(req_contig == NULL || req_start_pos == NULL || req_end_pos == NULL || req_mod_code == NULL) {
                response = "Invalid query get_contig_range_mod:<contig>:<start_pos>:<end_pos>:<mod_code>\n";
            } else {
                response = get_stats_contig_range_mod_code(req_contig, atoi(req_start_pos), atoi(req_end_pos), req_mod_code[0], map);
                if (response == NULL) {
                    response = "No matching data found\n";
                }
            }
        } else {
            response = "Invalid query\n";
        }
        
    }

    // and send that buffer to client 
    if (write(client_fd, response, strlen(response)*sizeof(char)) < 0) {
        ERROR("%s", "error writing to socket\n");
    }
    
    close(client_fd);
    pthread_exit(0);
}

void * start_server(void* arg) {
    server_args_t *args = (server_args_t *)arg;
    int port = args->port;
    map = args->freq_map;
    int sockfd;
    struct sockaddr_in servaddr;

    // socket create and verification
    double t1 = realtime();
    while (realtime() - t1 < TIMEOUT) {
        sockfd = socket(AF_INET, SOCK_STREAM, 0);
        if (sockfd == -1) {
            ERROR("%s", "socket creation failed...retrying in 5s\n");
            sleep(5);
        }
        else {
            break;
        }
    }
    
    memset(&servaddr, 0, sizeof(servaddr));

    // assign IP, PORT
    servaddr.sin_family = AF_INET;
    servaddr.sin_addr.s_addr = htonl(INADDR_ANY);
    servaddr.sin_port = htons(port);

    // Binding newly created socket to given IP and verification
    t1 = realtime();
    while(realtime() - t1 < TIMEOUT) {
        if ((bind(sockfd, (struct sockaddr*)&servaddr, sizeof(servaddr))) != 0) {
            ERROR("socket bind to port %d failed...retrying in 5s\n", port);
            sleep(5);
        }
        else {
            break;
        }
    }
    

    // Now server is ready to listen and verification
    t1 = realtime();
    while(realtime() - t1 < TIMEOUT) {
        if ((listen(sockfd, 5)) != 0) {
            ERROR("listen on port %dfailed...\n", port);
            sleep(5);
        }
        else{
            INFO("server listening on port %d..\n", port);
            break;
        }
            
    }

    // Accept the data packet from client and verification
    while(1) {
        // Accept the data packet from client and verification
        struct sockaddr_in cli;
        u_int32_t client_fd;
        u_int32_t len;
        len = sizeof(cli);
        client_fd = accept(sockfd, (struct sockaddr*)&cli, &len);
        if (client_fd < 0) {
            ERROR("server acccept failed for client %d...\n", client_fd);
            sleep(5);
        }
        else {
            pthread_t client_thread;
            pthread_create(&client_thread, NULL, handle_request, (void*) &client_fd);
            pthread_detach(client_thread);
        }
        
    }

    close(sockfd);
    pthread_exit(0);

    return NULL;
}
