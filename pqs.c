/*
 * pqs.c
 *
 *  Created on: Nov 15, 2011
 *      Author: helder
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "mpi.h"

typedef unsigned short uint16_t;
typedef unsigned long uint64_t;

typedef struct{
	uint64_t size;
	uint16_t *element;
} array_t;

#define PQS_ERROR -1
#define PQS_SUCCESS 0

#define PQS_ROOT_RANK 0

//#define DEBUG

int initArray(array_t **a, uint64_t size);
int distributeArray(array_t **ap, array_t *a, uint64_t sizeArray, int numProcs);
int gatherSamples(array_t **allSamples, array_t *arrayPart, int numProcs, int rank, uint64_t sizeArray);
int broadcastPivots(array_t **piv, array_t *allSamples, int numPivots, int rank);
int sendReceivePartitions(array_t **partitions, uint64_t **displacements, array_t *arrayPart, array_t *pivots, int numProcs, int rank);
int mergePartitions(array_t **orderedArray, array_t *p, uint64_t *displacements, int numProcs, int rank);

int chooseSamples(array_t **chosenSamples, array_t *list, int numProcs, uint64_t sizeArray);
int choosePivots(array_t **chosenPivots, array_t *list, int numProcs);
int binarySearch(array_t *array, uint16_t key, uint64_t start);
int compare (const void * a, const void * b);

int main(int argc, char **argv){

	array_t *array = NULL;
	array_t *arrayPart = NULL;
	array_t *allSamples = NULL;
	array_t *pivots = NULL;

	array_t *partitions = NULL;
	uint64_t *displacements = NULL;

	array_t *orderedArray = NULL;

	int rank = 0, numProcs = 0;
	uint64_t i;
	uint64_t sizeArray = atol(argv[1]);

	if(argc != 2 || sizeArray <= 0){
		printf("Invalid args!\n Usage: ./pqs <array_size>\n");
		return PQS_ERROR;
	}

	MPI_Init(&argc, &argv);

	//get my rank
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	//get the size of the comm
	MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

	//the size of the array must be greater than and divisible by the number of process, for the sake of simplicity
	if(sizeArray < numProcs || sizeArray % numProcs != 0){
		printf("The size of the array must be greater than and divisible by the number of processes!\n Usage: ./pqs <array_size>\n");
		MPI_Finalize();
		return PQS_ERROR;
	}

	//fill the array with random numbers
	if(rank == PQS_ROOT_RANK){
		if(initArray(&array, sizeArray) == PQS_ERROR)
			return PQS_ERROR;
	}

	//distribute the array between the running process
	distributeArray(&arrayPart, array, sizeArray, numProcs);

	if(rank == PQS_ROOT_RANK){
		free(array->element); free(array); array = NULL;
	}

	//each process will sort its part of the array using qsort
	qsort(arrayPart->element, arrayPart->size, sizeof(uint16_t), compare);

#ifdef DEBUG
	for(i=0; i<arrayPart->size; i++){
		printf("%s():%d - I'm [%d], el[%lu] = [%hu]\n", __FUNCTION__, __LINE__, rank, i, arrayPart->element[i]);
	}
#endif

	//each process will select numProcs samples and send to the root process
	gatherSamples(&allSamples, arrayPart, numProcs, rank, sizeArray);

#ifdef DEBUG
	if(rank == PQS_ROOT_RANK){
		for(i=0; i<allSamples->size; i++){
			printf("%s():%d - I'm [%d], allSamples[%lu] = [%hu]\n", __FUNCTION__, __LINE__, rank, i, allSamples->element[i]);
		}
	}
#endif

	broadcastPivots(&pivots, allSamples, numProcs, rank);

#ifdef DEBUG
	if(rank == PQS_ROOT_RANK){
		for(i=0; i<pivots->size; i++){
			printf("%s():%d - I'm [%d], pivots[%lu] = [%hu]\n", __FUNCTION__, __LINE__, rank, i, pivots->element[i]);
		}
	}
#endif

	if(rank == PQS_ROOT_RANK){
		free(allSamples->element); free(allSamples); allSamples = NULL;
	}

	//send my partitions to each process and receive partitions from all processes
	sendReceivePartitions(&partitions, &displacements, arrayPart, pivots, numProcs, rank);

	free(pivots->element); free(pivots); pivots = NULL;

#ifdef DEBUG
	for(i=0; i<partitions->size; i++){
		printf("%s():%d - I'm [%d], partitions->el[%lu] = [%hu]\n",
				__FUNCTION__, __LINE__, rank, i, partitions->element[i]);
	}

	for(i=0; i<numProcs; i++){
		printf("%s():%d - I'm [%d], displacements[%lu] = [%lu]\n",
				__FUNCTION__, __LINE__, rank, i, displacements[i]);
	}
#endif

	//free'ing arrayPart, since 'partitions' already has all elements that we need
	free(arrayPart->element); free(arrayPart); arrayPart = NULL;

	mergePartitions(&orderedArray, partitions, displacements, numProcs, rank);

	//print output
	for(i=0; i<orderedArray->size; i++){
		printf("%s():%d - I'm [%d], orderedArray->el[%lu] = [%hu]\n",
			__FUNCTION__, __LINE__, rank, i, orderedArray->element[i]);
	}

	MPI_Finalize();

	return PQS_SUCCESS;
}

int initArray(array_t **a, uint64_t size){
	array_t *array = NULL;
	int i;

	if((array = (array_t *)malloc(sizeof(array_t))) == NULL){
		printf("Malloc error!\n");
		return PQS_ERROR;
	}

	array->size = size;
	if((array->element = (uint16_t *)malloc(sizeof(uint16_t) * size)) == NULL){
		printf("Malloc error!\n");
		return PQS_ERROR;
	}

	srand(time(NULL));

	//fill array with random numbers
	for(i=0; i<size; i++){
		array->element[i] = rand() % 65534; //max value that uint16_t support
	}

	*a = array;
	return PQS_SUCCESS;
}

int distributeArray(array_t **ap, array_t *a, uint64_t sizeArray, int numProcs){
	array_t *arrayPart = NULL;
	int numElementsPerProcess = 0;

	//number of elements that each process will receive
	numElementsPerProcess = sizeArray/numProcs;

	//malloc the buffer var that will hold part of the array
	if((arrayPart = (array_t *)malloc(sizeof(array_t))) == NULL){
		printf("Malloc error!\n");
		return PQS_ERROR;
	}

	arrayPart->size = numElementsPerProcess;
	if((arrayPart->element = (uint16_t *)malloc(sizeof(uint16_t) * numElementsPerProcess)) == NULL){
		printf("Malloc error!\n");
		return PQS_ERROR;
	}

	MPI_Scatter((a!=NULL? a->element : NULL), numElementsPerProcess, MPI_UNSIGNED_SHORT, arrayPart->element, numElementsPerProcess, MPI_UNSIGNED_SHORT, 0, MPI_COMM_WORLD);

	*ap = arrayPart;
	return PQS_SUCCESS;
}

int gatherSamples(array_t **as, array_t *arrayPart, int numProcs, int rank, uint64_t sizeArray){
	array_t *samples = NULL;
	array_t *allSamples = NULL;

	srand(time(NULL));

	if(chooseSamples(&samples, arrayPart, numProcs, sizeArray) != PQS_SUCCESS)
		return PQS_ERROR;

#ifdef DEBUG
	uint64_t i;
	for(i=0; i<samples->size; i++){
		printf("%s():%d - I'm [%d], chosenSamples[%lu] = [%hu]\n", __FUNCTION__, __LINE__, rank, i, samples->element[i]);
	}
#endif

	//malloc receiving buffer
	if(rank == PQS_ROOT_RANK){

		if((allSamples = (array_t*) malloc(sizeof(array_t))) == NULL){
			printf("Malloc error!\n");
			return PQS_ERROR;
		}

		allSamples->size = samples->size  * numProcs;

		if((allSamples->element = (uint16_t*) malloc(sizeof(uint16_t) * allSamples->size)) == NULL){
			printf("Malloc error!\n");
			return PQS_ERROR;
		}
	}

	MPI_Gather(samples->element, samples->size, MPI_UNSIGNED_SHORT, rank==0?allSamples->element:NULL, samples->size, MPI_UNSIGNED_SHORT, 0, MPI_COMM_WORLD);

	free(samples->element);
	free(samples);

	//sort the allSamples list
	if(rank == PQS_ROOT_RANK){
		qsort(allSamples->element, allSamples->size, sizeof(uint16_t), compare);
	}

	 *as = allSamples;
	return PQS_SUCCESS;

}

int broadcastPivots(array_t **piv, array_t *allSamples, int numProcs, int rank){
	array_t *pivots = NULL;

	//choose pivots from allSamples
	if(choosePivots(&pivots, allSamples, numProcs) != PQS_SUCCESS)
		return PQS_ERROR;

	MPI_Bcast(pivots->element, numProcs-1, MPI_UNSIGNED_SHORT, 0, MPI_COMM_WORLD);

	pivots->size = numProcs-1;

	*piv = pivots;

	return PQS_SUCCESS;
}

// qsort compare function
int compare (const void * a, const void * b){
	return ( *(uint16_t*)a - *(uint16_t*)b );
}

int chooseSamples(array_t **chosenSamples, array_t *list, int numProcs, uint64_t sizeArray){
	array_t *samples = NULL;
	int i, /*aux,*/ index;

	//malloc the buffer var that will hold part of the array
	if((samples = (array_t *)malloc(sizeof(array_t))) == NULL){
		printf("Malloc error!\n");
		return PQS_ERROR;
	}

	samples->size = 0;

	if((samples->element = (uint16_t*) malloc(sizeof(uint16_t) * numProcs)) == NULL){
		printf("Malloc error!\n");
		return PQS_ERROR;
	}

	//select samples from list
	double aux = ceil((double)sizeArray / (double)(numProcs * numProcs));

	for(i=0; i<numProcs; i++){
		index = i * aux;

		if(index > list->size-1)
			break;

		samples->element[i] = list->element[index];
		samples->size++;
	}

	*chosenSamples = samples;

	return PQS_SUCCESS;
}

int choosePivots(array_t **chosenPivots, array_t *list, int numProcs){
	array_t *pivots = NULL;
	int i, aux, index;

	//malloc the buffer var that will hold part of the array
	if((pivots = (array_t *)malloc(sizeof(array_t))) == NULL){
		printf("Malloc error!\n");
		return PQS_ERROR;
	}

	pivots->size = 0;

	if((pivots->element = (uint16_t*) malloc(sizeof(uint16_t) * numProcs-1)) == NULL){
		printf("Malloc error!\n");
		return PQS_ERROR;
	}

	//if there is no list, this is not the root process. only malloc it to use it later through mpi calls
	if(list == NULL){
		*chosenPivots = pivots;

		return PQS_SUCCESS;
	}

	//select samples from list
	aux = floor(numProcs/2) - 1;
	for(i=1; i<=numProcs-1; i++){
		index = (i * numProcs) + aux;

		pivots->element[i-1] = list->element[index];
		pivots->size++;
	}

	*chosenPivots = pivots;

	return PQS_SUCCESS;
}

int sendReceivePartitions(array_t **partitions, uint64_t **displacements, array_t *arrayPart, array_t *pivots, int numProcs, int rank){
	uint64_t *sendCounts, *recvCounts, *sendDispls, *recvDispls;
	array_t *recvPartitions = NULL;

	uint64_t init, index, i;
	uint64_t sizeAllPartitions = 0;

	//malloc memory to mpi_alltoallv args
	sendCounts = (uint64_t*) malloc(numProcs * sizeof(uint64_t));
	recvCounts = (uint64_t*) malloc(numProcs * sizeof(uint64_t));
	sendDispls = (uint64_t*) malloc(numProcs * sizeof(uint64_t));
	recvDispls = (uint64_t*) malloc(numProcs * sizeof(uint64_t));

	if(sendCounts == NULL || recvCounts == NULL ||
			sendDispls == NULL|| recvDispls == NULL){

		printf("%s():%d - I'm [%d], malloc error!\n", __FUNCTION__, __LINE__, rank);
		return PQS_ERROR;
	}

	//traverse sample list
	init = 0;
	for(i=0; i<pivots->size; i++){
		index = binarySearch(arrayPart, pivots->element[i], init);
		if(index == PQS_ERROR){
			printf("%s():%d - I'm [%d], error while executing binarySearch()!\n", __FUNCTION__, __LINE__, rank);
			return PQS_ERROR;
		}

#ifdef DEBUG
		//printf("%s():%d - I'm [%d], index = [%lu], init = [%lu]\n", __FUNCTION__, __LINE__, rank, index, init);
#endif

		sendCounts[i] = index - init + 1;
		sendDispls[i] = (i!=0)? sendDispls[i-1] + sendCounts[i-1]: 0;

#ifdef DEBUG
		printf("%s():%d - I'm [%d], sendCount[%lu] = [%lu], sendDispls[%lu] = [%lu]\n",
				__FUNCTION__, __LINE__, rank, i, sendCounts[i], i, sendDispls[i]);
#endif

		init = index + 1; //preset for next pivot
	}

	//set the sendCounts and sendDispls for the last process
	sendCounts[i] = arrayPart->size - (init+1>=arrayPart->size? arrayPart->size: init);
	sendDispls[i] = sendDispls[i-1] + sendCounts[i-1];

#ifdef DEBUG
	printf("%s():%d - I'm [%d], sendCount[%lu] = [%lu], sendDispls[%lu] = [%lu]\n",
			__FUNCTION__, __LINE__, rank, i, sendCounts[i], i, sendDispls[i]);
#endif

	//distribute the counts for each process
	MPI_Alltoall(sendCounts, 1, MPI_UNSIGNED_LONG, recvCounts, 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);


#ifdef DEBUG
	//for(i=0; i<numProcs; i++){
	//	printf("%s():%d - I'm [%d], recvCounts[%lu] = [%lu]\n",
	//			__FUNCTION__, __LINE__, rank, i, recvCounts[i]);
	//}
#endif

	//set recv args to mpi_alltoallv
	for(i=0; i<numProcs; i++){
		sizeAllPartitions += recvCounts[i];
		recvDispls[i] = (i>0? recvDispls[i-1] + recvCounts[i-1]: 0);
	}
	recvPartitions = (array_t*) malloc(sizeof(array_t));
	recvPartitions->size = sizeAllPartitions;
	recvPartitions->element = (uint16_t*)malloc(sizeAllPartitions * sizeof(uint16_t));


#ifdef DEBUG
	/*printf("%s():%d - I'm [%d], sizeAllPartitions = [%lu]\n",
			__FUNCTION__, __LINE__, rank, sizeAllPartitions);

	for(i=0; i<numProcs; i++){
		printf("%s():%d - I'm [%d], recvDispls[%lu] = [%lu]\n",
				__FUNCTION__, __LINE__, rank, i, recvDispls[i]);
	}*/
#endif

	//TODO: for some reason, this fucntion requires that the count and displacements arrays be an int *
	//			but i cant understand why. For really large arrays, this might be a problem
	MPI_Alltoallv(arrayPart->element, sendCounts, sendDispls, MPI_UNSIGNED_SHORT,
			recvPartitions->element, recvCounts, recvDispls, MPI_UNSIGNED_SHORT,
			MPI_COMM_WORLD);

	//TODO: REMOVE this, just to debug
	//qsort(recvPartitions->element, recvPartitions->size, sizeof(uint16_t), compare);

	//cleanup
	free(sendCounts);
	free(recvCounts);
	free(sendDispls);

	*partitions = recvPartitions;
	*displacements = recvDispls;

	return PQS_SUCCESS;

}

int mergePartitions(array_t **ordArray, array_t *p, uint64_t *displacements, int numProcs, int rank){
	array_t **partitionList = NULL;
	array_t *partition = NULL;
	array_t *orderedArray = NULL;

	uint16_t *lesserElement = NULL;

	int i, parIndex = 0;
	uint64_t j;

	if((partitionList = (array_t**) malloc(numProcs * sizeof(array_t*))) == NULL){
		MPI_Abort(MPI_COMM_WORLD, PQS_ERROR);
		return PQS_ERROR;
	}

	//partitionList->size = numProcs;

	//divide the full partition list into 'numProcs' lists
	for(i=0; i<numProcs; i++){

		if((partition = (array_t *) malloc(sizeof(array_t))) == NULL){
			MPI_Abort(MPI_COMM_WORLD, PQS_ERROR);
			return PQS_ERROR;
		}

		partition->size = (i+1<numProcs? displacements[i+1] : p->size) - displacements[i];
		partition->element = &p->element[displacements[i]];

		partitionList[i] = partition;

#ifdef DEBUG
		printf("%s():%d - I'm [%d], partition[%d]->el = [%hu], size = [%lu]\n",
						__FUNCTION__, __LINE__, rank, i, partition->size!=0? partition->element[0] : 0, partition->size);
#endif
	}

	//malloc memory to the ordered array
	if((orderedArray = (array_t*) malloc(sizeof(array_t))) == NULL){
		MPI_Abort(MPI_COMM_WORLD, PQS_ERROR);
		return PQS_ERROR;
	}

	orderedArray->size = p->size;
	if((orderedArray->element = (uint16_t*) malloc(p->size * sizeof(uint16_t))) == NULL){
		MPI_Abort(MPI_COMM_WORLD, PQS_ERROR);
		return PQS_ERROR;
	}

	//create the new ordered list by comparing the partitions
	for(j=0; j<orderedArray->size; j++){

		for(i=0; i<numProcs; i++){

			if((lesserElement == NULL && partitionList[i]->size > 0) || (partitionList[i]->size > 0 && partitionList[i]->element[0] < *lesserElement)){

				lesserElement = &partitionList[i]->element[0];
				parIndex = i;
			}
		}

		partitionList[parIndex]->size--;
		partitionList[parIndex]->element++;

		orderedArray->element[j] = *lesserElement;

		lesserElement = NULL;
	}

#ifdef DEBUG
	for(j=0; j<orderedArray->size; j++){
		printf("%s():%d - I'm [%d], orderedArray->el[%lu] = [%hu]\n",
			__FUNCTION__, __LINE__, rank, j, orderedArray->element[j]);
	}


#endif

	//cleanup
	for(i=0; i<numProcs; i++){
		free(partitionList[i]);
	}
	free(partitionList);

	*ordArray = orderedArray;

	return PQS_SUCCESS;
}


//TODO: our array can have equal values, MAYBE we need to check that before returning the index!
int binarySearch(array_t *array, uint16_t key, uint64_t start){
	uint64_t low, high, mid;

	low = start; high = array->size-1;
	while(low <= high) {
		mid = ceil(low+high)/2;

		if(key < array->element[mid]){
			//if the array doesnt have the exact value, return the closest one
			if(key >= array->element[mid-1])
				return mid-1;
			else
				high = mid-1;

		}else if(key > array->element[mid]){

			//if the array doesnt have the exact value, return the closest one
			if(key < array->element[mid+1])
				return mid;
			else
				low = mid+1;
		}
		else return mid; /* found */
	}

#ifdef DEBUG
	printf("%s():%d - Key not found happened, el[0] = [%hu]\n", __FUNCTION__, __LINE__, array->element[0]);
#endif

	//key not found, it is greater or lesser than all elements of the array
	if(key < array->element[start])
		return start;
	else
		return array->size-1;

	//return PQS_ERROR;
}
