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

struct Star star_array[ NUM_STARS ];
uint8_t   (*distance_calculated)[NUM_STARS];

double  min_global  = FLT_MAX;
double  max_global  = FLT_MIN;

int thread_ct=1;
int load = NUM_STARS;

typedef struct
{
  uint32_t start;
  uint32_t end;
  uint64_t count;
  double dist_sum;
  double min;
  double max;
} ThreadArg;

void * determineAverageAngularDistance( void * arg )
{
	ThreadArg *this = (ThreadArg*) arg;
	uint32_t i, j, end = this->end;
	double distance,count=0,dist_sum=0;

	for (i = this->start; i < end; i++)
	{
	  for (j = i+1; j < NUM_STARS; j++)
	  {
		if( i!=j && distance_calculated[i][j] == 0 )
		{
			distance = calculateAngularDistance( star_array[i].RightAscension, star_array[i].Declination,
														star_array[j].RightAscension, star_array[j].Declination ) ;
			
			distance_calculated[i][j] = 1;
			distance_calculated[j][i] = 1;
			
			count++;
			dist_sum += distance;

			if( this->min > distance ) this->min = distance;
			if( this->max < distance ) this->max = distance;
			
		}
	  }
	}
	this->dist_sum = dist_sum;
	this->count = count;
	pthread_exit(NULL);
}

double runThreads()
{
	pthread_t threads[thread_ct];
	ThreadArg *args[thread_ct];
	uint32_t start=0;
	int i;
	for(i=0; i< thread_ct; i++)
	{
		args[i] = malloc(sizeof(ThreadArg));
		start = i * load;
		args[i]->start = start;
		args[i]->end = start + load;
		args[i]->min = FLT_MAX;
		args[i]->max = FLT_MIN;
		pthread_create(&threads[i], NULL, determineAverageAngularDistance, (void *) args[i]);
	}
	
	double dist_sum=0,count_sum=0;
	for(i=0;i<thread_ct;i++)
	{
		pthread_join(threads[i], NULL);
		if (args[i]->min<min_global) min_global = args[i]->min;
		if (args[i]->max>max_global) max_global = args[i]->max;
		count_sum += args[i]->count;
		dist_sum += args[i]->dist_sum;
		free(args[i]);
	}
	free(distance_calculated);

	return dist_sum / count_sum;
}


int main( int argc, char * argv[] )
{
	struct timespec start, end;
	clock_gettime(CLOCK_MONOTONIC, &start);
	FILE *fp;
	uint32_t star_count=0,n,column;

	distance_calculated = malloc(sizeof(uint8_t[NUM_STARS][NUM_STARS]));
	memset(distance_calculated,0,sizeof(uint8_t[NUM_STARS][NUM_STARS]));

	for( n = 1; n < argc; n++ )          
	{
		if( strcmp(argv[n], "-t")==0)
		{
			thread_ct = atoi(argv[n+1]);
			load = load/thread_ct;
		}
	}

	fp = fopen( "data/tycho-trimmed.csv", "r" );

	if( fp == NULL )
	{
		printf("ERROR: Unable to open the file data/tycho-trimmed.csv\n");
		exit(1);
	}

	char line[MAX_LINE];
	while (fgets(line, 1024, fp))
	{
		column=0;
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

	double dist = runThreads();

	clock_gettime(CLOCK_MONOTONIC, &end);

	printf("Average distance found is %f\n", dist );
	printf("Minimum distance found is %lf\n", min_global );
	printf("Maximum distance found is %lf\n", max_global );
	printf("Seconds taken: %lf\n", (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1000000000.0);
	return 0;
}
