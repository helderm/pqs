pqs
===

A parallel implementation of Quicksort using MPI

How to use
===

1) Install and configure MPICH
sudo apt-get install mpich2

2) Build
mpicc -o pqs pqs.c -lm
gcc -o qsort qsort.c

3) Debug Build
mpicc -lm -DDEBUG -o pqs pqs.c
gcc -DDEBUG -o qsort qsort.c

3) Run
time mpiexec -n <#process> ./pqs <size_array>
time qsort <size_array>
