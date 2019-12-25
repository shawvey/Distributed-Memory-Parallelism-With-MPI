# Distributed-Memory-Parallelism-With-MPI
this assignment uses MPI “Single Program Multiple Data (SPMD)” distributed memory parallelism to run stencil code from one core up to all cores of 2 nodes.


## Building and running the code

1. you have to run the following command to configure environment.

   ```shell
   source env.sh
   ```

2. A simple `Makefile` is provided to build the code using the GCC compiler.  Just
   type `make` to build the code.  A binary named `stencil` is produced on a
   successful build.

3. There are *three* input problems tested, representing different grid sizes.  The
   inputs are all set as command line arguments to the executable in the following
   format:

   ```shell
   ./stencil nx ny niters
   ```

4. The inputs required are:

| nx   | ny   | niters | Command                   |
| ---- | ---- | ------ | ------------------------- |
| 1024 | 1024 | 100    | `./stencil 1024 1024 100` |
| 4096 | 4096 | 100    | `./stencil 4096 4096 100` |
| 8000 | 8000 | 100    | `./stencil 8000 8000 100` |


## Checking the results

The program will have executed correctly if the output image matches the
provided reference output images with a small tolerance of +/- 1.  A Python
check script is provided to check the values. 

    python check.py --ref-stencil-file stencil_1024_1024_100.pgm --stencil-file stencil.pgm

If any errors are found, the script can be rerun with the addition of the
`--verbose` flag to provide information about where the errors occur in the
grid.

The reference input files for the different problems are named:

| nx   | ny   | niters | Reference file              |
| ---- | ---- | ------ | --------------------------- |
| 1024 | 1024 | 100    | `stencil_1024_1024_100.pgm` |
| 4096 | 4096 | 100    | `stencil_4096_4096_100.pgm` |
| 8000 | 8000 | 100    | `stencil_8000_8000_100.pgm` |


## Final results

|           | Before optimisation | After optimisation |
| :-------- | :------------------ | :------------------ |
| 1024*1024 | 0.110686 s          | 0.028604 s          |
| 4096*4096 | 2.878536 s          | 0.180593 s          |
| 8000*8000 | 10.348350 s         | 1.241563 s          |
