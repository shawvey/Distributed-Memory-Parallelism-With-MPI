stencil: stencil.c
	mpicc -Ofast -mtune=native -std=c99 -Wall $^ -o $@


