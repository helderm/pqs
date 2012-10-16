pqs
===

A parallel implementation of Quicksort using MPI

How to use
----------

- Install and configure MPICH

sudo apt-get install mpich2

- Build

mpicc -o pqs pqs.c -lm

gcc -o qsort qsort.c

- Debug Build

mpicc -lm -DDEBUG -o pqs pqs.c

gcc -DDEBUG -o qsort qsort.c

- Run

time mpiexec -n <#process> ./pqs <size_array>

time qsort <size_array>
