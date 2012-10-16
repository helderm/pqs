/*
 * qsort.c
 *
 *  Created on: Dec 6, 2011
 *      Author: helder
 */

#include <stdio.h>
#include <stdlib.h>

typedef unsigned short uint16_t;
typedef unsigned long uint64_t;

typedef struct{
	uint64_t size;
	uint16_t *element;
} array_t;

//#define DEBUG

int initArray(array_t **a, int size);
int compare (const void * a, const void * b);

int main(int argc, char **argv){
	array_t *array = NULL;
	int sizeArray = atol(argv[1]);

	if(argc != 2 || sizeArray <= 0){
		printf("Invalid args!\n Usage: ./qsort <array_size>\n");
		return 1;
	}

	if(initArray(&array, sizeArray) == 1)
		return 1;

	qsort(array->element, array->size, sizeof(uint16_t), compare);

#ifdef DEBUG
	int i;
	for(i=0; i<array->size; i++){
		printf("%s():%d - el[%d] = [%d]\n", __FUNCTION__, __LINE__, i, array->element[i]);
	}
#endif

	return 0;
}

int initArray(array_t **a, int size){
	array_t *array = NULL;
	int i;

	if((array = (array_t *)malloc(sizeof(array_t))) == NULL){
		printf("Malloc error!\n");
		return 1;
	}

	array->size = size;
	if((array->element = (uint16_t *)malloc(sizeof(uint16_t) * size)) == NULL){
		printf("Malloc error!\n");
		return 1;
	}

	srand(time(NULL));

	//fill array with random numbers
	for(i=0; i<size; i++){
		array->element[i] = rand() % 65534; //max value that uint16_t support
	}

	*a = array;
	return 0;
}

// qsort compare function
int compare (const void * a, const void * b){
	return ( *(uint16_t*)a - *(uint16_t*)b );
}
