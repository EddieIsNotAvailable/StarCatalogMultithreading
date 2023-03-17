// MIT License
// 
// Copyright (c) 2023 Trevor Bakker 
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <getopt.h>
#include <stdint.h>
#include "utility.h"
#include "star.h"
#include "float.h"
#include <pthread.h>

#define NUM_STARS 30000 
#define MAX_LINE 1024
#define DELIMITER " \t\n"

struct Star star_array[ NUM_STARS ];
uint8_t   (*distance_calculated)[NUM_STARS];

pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;

double min_global = FLT_MAX;
double max_global = FLT_MIN;

int thread_ct=4;

typedef struct {
  uint32_t start;
  uint32_t end;
  uint64_t count;
  double mean;
  double min;
  double max;
} ThreadArg;

void * work(void * arg) {
  ThreadArg *this = (ThreadArg*) arg;
  uint32_t i,j;
  int end = this->end;
  int start = this->start;
  double sum=0;
  double test;
  for (i = start; i < end; i++)
  {
    for (j = start; j < end; j++)
    {
      pthread_mutex_lock(&mutex);
      int a = distance_calculated[i][j];
      pthread_mutex_unlock(&mutex);

      if( i!=j && a == 0 )
      {
        double distance = calculateAngularDistance( star_array[i].RightAscension, star_array[j].Declination,
                                                    star_array[j].RightAscension, star_array[j].Declination ) ;
        pthread_mutex_lock(&mutex);
        distance_calculated[i][j] = 1;
        distance_calculated[j][i] = 1;
        pthread_mutex_unlock(&mutex);
        this->count++;

        if( this->min > distance )
        {
          this->min = distance;
        }

        if( this->max < distance )
        {
          this->max = distance;
        }
        //this->mean = this->mean + (distance - this->mean)/this->count;
        //this->mean += distance;
        //sum += distance;

        test = this->mean + (distance - this->mean)/this->count;

        if(isnan(test)) {
          printf("Thread w start %d at i: %d, j: %d Resulted in nan -  mean: %f dist: %f count: %lu\n",this->start, i,j,this->mean,distance,this->count);
          printf("Array info there: i.ra= %f, i.des = %f, j.ra = %f, j.d = %f\n",star_array[i].RightAscension,star_array[i].Declination,star_array[j].RightAscension,star_array[j].Declination);
        }
        else this->mean = test;
        

      }
    }
    
    
  }
  //*(ThreadArg*) arg = *this;
  //pthread_exit(arg);
  pthread_exit(NULL);
}


double runThreads() {
  pthread_t threads[thread_ct];
  ThreadArg *args[thread_ct];
  int i;

  int load = NUM_STARS/thread_ct;
  int start;


  for(i=0; i<thread_ct; i++)
  {
    start = load * i;
    args[i] = malloc(sizeof(ThreadArg));
    args[i]->count=0;
    args[i]->start = start;
    args[i]->end = start + load;
    args[i]->mean = 0;
    args[i]->min = FLT_MAX;
    args[i]->max = FLT_MIN;
    
    if(pthread_create(&threads[i], NULL, &work, (void *) args[i]) != 0)
      perror("Thread creation failed");
  }

  uint64_t count_sum=0;
  double mean_sum=0,min,max;
  for(i=0; i<thread_ct; i++)
  {
    pthread_join(threads[i], NULL);
    //count_sum += args[i]->count;
    min = args[i]->min;
    max = args[i]->max;
    if(min_global>min) min_global = min;
    if(max_global<max) max_global = max;

    printf("Thread[%d] - Count: %lu Start: %d End: %d Mean: %f Min: %f Max: %f\n",
    i, args[i]->count, args[i]->start, args[i]->end, args[i]->mean, min, max);
    //mean_sum = mean_sum + (args[i]->mean - mean_sum) / count_sum;
    mean_sum += args[i]->mean;
  }
  printf("Mean Sum: %f\n",mean_sum);


  //return mean_sum/count_sum;
  return mean_sum/thread_ct;
}


int main( int argc, char * argv[] )
{
  struct timespec start, finish;
  clock_gettime(CLOCK_MONOTONIC, &start);

  FILE *fp;
  uint32_t star_count = 0;

  distance_calculated = malloc(sizeof(uint8_t[NUM_STARS][NUM_STARS]));

  if( distance_calculated == NULL )
  {
    uint64_t num_stars = NUM_STARS;
    uint64_t size = num_stars * num_stars * sizeof(uint8_t);
    printf("Could not allocate %ld bytes\n", size);
    exit( EXIT_FAILURE );
  }

  memset(distance_calculated,0,sizeof(uint8_t)*NUM_STARS*NUM_STARS);

  fp = fopen( "data/tycho-trimmed.csv", "r" );

  if( fp == NULL )
  {
    printf("ERROR: Unable to open the file data/tycho-trimmed.csv\n");
    exit(1);
  }

  char line[MAX_LINE];
  while (fgets(line, 1024, fp))
  {
    uint32_t column = 0;

    char* tok;
    for (tok = strtok(line, " ");
            tok && *tok;
            tok = strtok(NULL, " "))
    {
       switch( column )
       {
          case 0:
              star_array[star_count].ID = atoi(tok);
              break;
       
          case 1:
              star_array[star_count].RightAscension = atof(tok);
              break;
       
          case 2:
              star_array[star_count].Declination = atof(tok);
              break;

          default: 
             printf("ERROR: line %d had more than 3 columns\n", star_count );
             exit(1);
             break;
       }
       column++;
    }
    star_count++;
  }
  printf("%d records read\n", star_count );

  printf("Threads: %d\n",thread_ct);
  
  // Find the average angular distance in the most inefficient way possible
  double distance =  runThreads();
  pthread_mutex_destroy(&mutex);
  printf("Average distance found is %lf\n", distance );
  printf("Minimum distance found is %lf\n", min_global );
  printf("Maximum distance found is %lf\n", max_global );
  clock_gettime(CLOCK_MONOTONIC, &finish);
 
  double elapsed = (finish.tv_sec - start.tv_sec);
  elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
  printf("Seconds taken: %f\n",elapsed);
  return 0;
}
