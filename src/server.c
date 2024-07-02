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
#include "meth.h"
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
#define PORT 8080 
#define SA struct sockaddr

char* help = "Available commands:\n \
            \tget:<contig>:<position>:<mod_code>\n \
            \texit\n";

// test with:
// telnet 127.0.0.1 8080
// get:chr22:19981716:m
// exit

// Function designed for chat between client and server. 
void handle_request(int connfd) 
{ 
    char buff[MAX]; 
    int n; 
    // infinite loop for chat 
    for (;;) { 
        // check if the connfd
        if (connfd < 0) {
            ERROR("%s", "server acccept failed...\n");
            exit(EXIT_FAILURE);
        }

        bzero(buff, MAX); 
   
        // read the message from client and copy it in buffer 
        read(connfd, buff, sizeof(buff)); 


        printf("Request: %s\n", buff); 

        
        char* tok;
        char* request = strtok(buff, "\r\n");
        char* response = "Invalid request\n";

        if(strcmp(request, "exit") == 0) {
            return;
        } else if(strcmp(request, "help") == 0) {
            response = help;
        } else if (request == NULL) {
            response = "Invalid request\n";
        } else {
            tok = strtok(request, ":");
            if(tok == NULL) {
                response = "Invalid request\n";
            } else if(strcmp(tok, "get") == 0) {
                char * res_chr = strtok(NULL, ":");
                char * res_pos = strtok(NULL, ":");
                char * res_mod_code = strtok(NULL, ":");

                if(res_chr == NULL || res_pos == NULL || res_mod_code == NULL) {
                    response = "Invalid request get:<contig>:<position>:<mod_code>\n";
                } else {
                    printf("Response: get chr: %s, pos: %s, mod_code: %c\n", res_chr, res_pos, res_mod_code[0]);
                    response = get_all_stats_at_pos(res_chr, atoi(res_pos), res_mod_code[0]);
                    if (response == NULL) {
                        response = "No data found\n";
                    }
                }

            } else if(strcmp(tok, "list") == 0) {
                response = "to be implemented\n";
            }
        }

        

        bzero(buff, MAX); 
        n = 0;
        while ((buff[n++] = response[n]) != '\0') 
            ;
   
        // and send that buffer to client 
        write(connfd, buff, sizeof(buff)); 
        
        INFO("Message sent to client: %s", buff);

    } 
}


void start_server(){
    int sockfd;
    struct sockaddr_in servaddr;

    // socket create and verification
    while (1)
    {
        sockfd = socket(AF_INET, SOCK_STREAM, 0);
        if (sockfd == -1) {
            ERROR("%s", "socket creation failed...retrying in 5s\n");
            sleep(5);
        }
        else {
            INFO("%s", "Socket successfully created..\n");
            break;
        }
    }
    
    
    bzero(&servaddr, sizeof(servaddr));

    // assign IP, PORT
    servaddr.sin_family = AF_INET;
    servaddr.sin_addr.s_addr = htonl(INADDR_ANY);
    servaddr.sin_port = htons(PORT);

    // Binding newly created socket to given IP and verification
    while(1){
        if ((bind(sockfd, (struct sockaddr*)&servaddr, sizeof(servaddr))) != 0) {
            ERROR("%s", "socket bind failed...retrying in 5s\n");
            sleep(5);
        }
        else {
            INFO("%s", "Socket successfully binded..\n");
            break;
        }
    }
    

    // Now server is ready to listen and verification
    while(1){
        if ((listen(sockfd, 5)) != 0) {
            ERROR("%s", "Listen failed...\n");
            sleep(5);
        }
        else{
            INFO("%s", "Server listening..\n");
            break;
        }
            
    }
    
    struct sockaddr_in cli;
    int connfd, len;
    len = sizeof(cli);
    
    // Accept the data packet from client and verification
    connfd = accept(sockfd, (struct sockaddr*)&cli, &len);
    if (connfd < 0) {
        ERROR("%s", "server acccept failed...\n");
        exit(EXIT_FAILURE);
    }
    else
        INFO("%s", "server acccept the client...\n");
        
    // Function for chatting between client and server
    handle_request(connfd);

    // After chatting close the socket
    close(sockfd);
    INFO("%s", "Server closed\n");
}

