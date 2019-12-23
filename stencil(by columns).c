#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <string.h>
#include "mpi.h"


//Define output file name                                                                                                       jm
#define OUTPUT_FILE "stencil.pgm"
// Define MASTER variable
#define MASTER 0


int calc_ncols_from_rank(int rank, int size, int NCOLS);
void stencil(const int nx, const int ny, const int width, const int height,
             float* image, float* tmp_image);
void init_image(const int nx, const int ny, const int width, const int height,
                float* image, float* tmp_image);
void output_image(const char* file_name, const int nx, const int ny,
                  const int width, const int height, float* image);
double wtime(void);

int main(int argc, char* argv[])
{
  int ii,jj;             /* row and column indices for the grid */
  int rank;              /* the rank of this process */
  int left;              /* the rank of the process to the left */
  int right;             /* the rank of the process to the right */
  int size;              /* number of processes in the communicator */
  int tag = 0;           /* scope for adding extra information to a message */
  MPI_Status status;     /* struct used by MPI_Recv */
  int local_nrows;       /* number of rows apportioned to this rank */
  int local_ncols;       /* number of columns apportioned to this rank */
  float *w;             /* local temperature grid at time t     */
  float *sendbuf;       /* buffer to hold values to send */
  float *recvbuf;       /* buffer to hold received values */
  float *sendbuf_v;
  float *recvbuf_v;
  float *v;    /* temporary rank array */

  // Check usage
  if (argc != 4) {
    fprintf(stderr, "Usage: %s nx ny niters\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  int nx = atoi(argv[1]);       /* numbers of rows */
  int ny = atoi(argv[2]);       /* numbers of columns */
  int niters = atoi(argv[3]);   /* iteration times */
  // Initiliase problem dimensions from command line arguments
  int height = nx + 2;//rows+2=height
  int width = ny + 2; //columns+2=width

  

  // Allocate the image
  float* image = malloc(sizeof(float) * width * height);
  float* tmp_image = malloc(sizeof(float) * width * height);


  // Set the input image
  init_image(nx, ny, height, width, image, tmp_image);
  

  MPI_Init( &argc, &argv );
  MPI_Comm_size( MPI_COMM_WORLD, &size );
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  double mintic;
  double maxtoc;
  
  left = (rank == MASTER) ? (rank + size - 1) : (rank - 1);
  right = (rank + 1) % size;


  // Allocate the spaces
  local_nrows = nx + 2;
  local_ncols = calc_ncols_from_rank(rank, size, width);
  if (local_ncols < 1) {
    fprintf(stderr,"Error: too many processes:- local_ncols < 1\n");
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }


  if(rank == MASTER || rank == (size-1)){
    w = (float*)malloc(sizeof(float*) * local_nrows * (local_ncols+1));
    v = (float*)malloc(sizeof(float*) * local_nrows * (local_ncols+1));
    /* The last rank has the most columns apportioned.
       printbuf must be big enough to hold this number */ 
  }
  else{
    w = (float*)malloc(sizeof(float*) * (local_nrows) * (local_ncols+2));
    v = (float*)malloc(sizeof(float*) * (local_nrows) * (local_ncols+2));
  }

  sendbuf = (float*)malloc(sizeof(float) * local_nrows);
  recvbuf = (float*)malloc(sizeof(float) * local_nrows);
  sendbuf_v = (float*)malloc(sizeof(float) * local_nrows);
  recvbuf_v = (float*)malloc(sizeof(float) * local_nrows);
 
  //first rank
  if(rank == MASTER){
    for(ii=0;ii<local_nrows;ii++) { 
      for(jj=0; jj<local_ncols + 1; jj++) {       
        if(jj == local_ncols){
          w[ii * (local_ncols + 1) + jj] = 0.0f;
          v[ii * (local_ncols + 1) + jj] = 0.0f;
        }

        else{
          w[ii * (local_ncols + 1) + jj] = image[ii * width + jj]; 
          v[ii * (local_ncols + 1) + jj] = tmp_image[ii * width + jj];
        }
      /* halo cells */ }
    }
  }

  //last rank
  else if(rank == size-1){
    for(ii=0;ii<local_nrows;ii++) { 
        for(jj=0; jj<local_ncols + 1; jj++) { 
          if (jj == 0){
            w[ii * (local_ncols + 1) + jj] = 0.0f;
            v[ii * (local_ncols + 1) + jj] = 0.0f;
          }
          else{
            w[ii * (local_ncols + 1) + jj] = image[ii * width + (width-local_ncols-1) + jj]; 
            v[ii * (local_ncols + 1) + jj] = tmp_image[ii * width + (width-local_ncols-1) + jj];
          }
          /* halo cells */ }
    }
  }

  //给中间的ranks赋值
  else{
    for(ii=0;ii<local_nrows;ii++) {
        for(jj=0; jj<local_ncols + 2; jj++) { 
          if (jj == 0 || jj == local_ncols + 1){
            w[ii * (local_ncols + 2) + jj] = 0.0f;
            v[ii * (local_ncols + 2) + jj] = 0.0f;
          }
          else{
            w[ii * (local_ncols + 2) + jj] = image[ (ii * width) + (rank * local_ncols)-1 + jj]; 
            v[ii * (local_ncols + 2) + jj] = tmp_image[ (ii * width) + (rank * local_ncols)-1 + jj];
        }
      }
    }
  }





  /*
  ** halo exchange for the local grids w:
  ** - first send to the left and receive from the right,
  ** - then send to the right and receive from the left.
  ** for each direction:
  ** - first, pack the send buffer using values from the grid
  ** - exchange using MPI_Sendrecv()
  ** - unpack values from the recieve buffer into the grid
  */

  /* send to the left, receive from right */




  double tic = wtime();
  MPI_Reduce(&tic, &mintic, 1, MPI_DOUBLE, MPI_MIN, MASTER,
           MPI_COMM_WORLD);

  for (int t = 0; t < niters; ++t) {
    if(rank == MASTER){

      for(ii=0; ii < local_nrows; ii++){
          sendbuf[ii] = w[ii * (local_ncols + 1)]; //传0的那列过去
          sendbuf_v[ii] = w[ii * (local_ncols + 1) + local_ncols-1]; 
        }
      MPI_Sendrecv(sendbuf, local_nrows, MPI_FLOAT, left, tag,
         recvbuf, local_nrows, MPI_FLOAT, right, tag,
         MPI_COMM_WORLD, &status);

      MPI_Sendrecv(sendbuf_v, local_nrows, MPI_FLOAT, right, tag,
         recvbuf_v, local_nrows, MPI_FLOAT, left, tag,
         MPI_COMM_WORLD, &status);

      for(ii=0; ii < local_nrows; ii++){
          w[ii * (local_ncols + 1) + local_ncols] = recvbuf[ii];
          w[ii * (local_ncols + 1)] = recvbuf_v[ii];
        }

      stencil(nx, local_ncols-1, local_ncols + 1, local_nrows, w, v);

      for(ii=0; ii < local_nrows; ii++){
          sendbuf[ii] = v[ii * (local_ncols + 1)]; //传0的那列过去
          sendbuf_v[ii] = v[ii * (local_ncols + 1) + local_ncols-1]; 
        }
      MPI_Sendrecv(sendbuf, local_nrows, MPI_FLOAT, left, tag,
         recvbuf, local_nrows, MPI_FLOAT, right, tag,
         MPI_COMM_WORLD, &status);

      MPI_Sendrecv(sendbuf_v, local_nrows, MPI_FLOAT, right, tag,
         recvbuf_v, local_nrows, MPI_FLOAT, left, tag,
         MPI_COMM_WORLD, &status);

      for(ii=0; ii < local_nrows; ii++){
          v[ii * (local_ncols + 1) + local_ncols] = recvbuf[ii];
          v[ii * (local_ncols + 1)] = recvbuf_v[ii];
        }

      stencil(nx, local_ncols-1, local_ncols + 1, local_nrows, v, w);
    }

    else if(rank ==size-1){

      for(ii=0; ii < local_nrows; ii++){
        sendbuf[ii] = w[ii * (local_ncols + 1) + 1];
        sendbuf_v[ii] = w[ii * (local_ncols + 1) + local_ncols];//传0那列
      }

      MPI_Sendrecv(sendbuf, local_nrows, MPI_FLOAT, left, tag,
         recvbuf, local_nrows, MPI_FLOAT, right, tag,
         MPI_COMM_WORLD, &status);

    
      MPI_Sendrecv(sendbuf_v, local_nrows, MPI_FLOAT, right, tag,
         recvbuf_v, local_nrows, MPI_FLOAT, left, tag,
         MPI_COMM_WORLD, &status);

      for(ii=0; ii < local_nrows; ii++){
        w[ii * (local_ncols + 1) + local_ncols] = recvbuf[ii];
        w[ii * (local_ncols + 1)] = recvbuf_v[ii];
      }

      stencil(nx, local_ncols-1, local_ncols+1, local_nrows, w, v);

      for(ii=0; ii < local_nrows; ii++){
        sendbuf[ii] = v[ii * (local_ncols + 1) + 1];
        sendbuf_v[ii] = v[ii * (local_ncols + 1) + local_ncols];//传0那列
      }

      MPI_Sendrecv(sendbuf, local_nrows, MPI_FLOAT, left, tag,
         recvbuf, local_nrows, MPI_FLOAT, right, tag,
         MPI_COMM_WORLD, &status);

    
      MPI_Sendrecv(sendbuf_v, local_nrows, MPI_FLOAT, right, tag,
         recvbuf_v, local_nrows, MPI_FLOAT, left, tag,
         MPI_COMM_WORLD, &status);

      for(ii=0; ii < local_nrows; ii++){
        v[ii * (local_ncols + 1) + local_ncols] = recvbuf[ii];
        v[ii * (local_ncols + 1)] = recvbuf_v[ii];
      }
      stencil(nx, local_ncols-1, local_ncols+1, local_nrows, v, w);
      
    }

    else{
      for(ii=0; ii < local_nrows; ii++){
          sendbuf[ii] = w[ii * (local_ncols + 2) + 1];
          sendbuf_v[ii] = w[ii * (local_ncols + 2) + local_ncols];
        }
      MPI_Sendrecv(sendbuf, local_nrows, MPI_FLOAT, left, tag,
         recvbuf, local_nrows, MPI_FLOAT, right, tag,
         MPI_COMM_WORLD, &status);


      MPI_Sendrecv(sendbuf_v, local_nrows, MPI_FLOAT, right, tag,
         recvbuf_v, local_nrows, MPI_FLOAT, left, tag,
         MPI_COMM_WORLD, &status);

      for(ii=0; ii < local_nrows; ii++){
          w[ii * (local_ncols + 2) + local_ncols + 1] = recvbuf[ii];
          w[ii * (local_ncols + 2)] = recvbuf_v[ii];
        }

      stencil(nx, local_ncols, local_ncols+2, local_nrows, w, v);

      for(ii=0; ii < local_nrows; ii++){
          sendbuf[ii] = v[ii * (local_ncols + 2) + 1];
          sendbuf_v[ii] = v[ii * (local_ncols + 2) + local_ncols];
        }
      MPI_Sendrecv(sendbuf, local_nrows, MPI_FLOAT, left, tag,
         recvbuf, local_nrows, MPI_FLOAT, right, tag,
         MPI_COMM_WORLD, &status);


      MPI_Sendrecv(sendbuf_v, local_nrows, MPI_FLOAT, right, tag,
         recvbuf_v, local_nrows, MPI_FLOAT, left, tag,
         MPI_COMM_WORLD, &status);

      for(ii=0; ii < local_nrows; ii++){
          v[ii * (local_ncols + 2) + local_ncols + 1] = recvbuf[ii];
          v[ii * (local_ncols + 2)] = recvbuf_v[ii];
        }

      stencil(nx, local_ncols, local_ncols+2, local_nrows, v, w);
    }
  }

  double toc = wtime();
  MPI_Reduce(&toc, &maxtoc, 1, MPI_DOUBLE, MPI_MAX, MASTER,
           MPI_COMM_WORLD);

  if(rank == MASTER){
    printf("------------------------------------\n");
    printf(" runtime: %lf s\n",  maxtoc - mintic);
    printf("------------------------------------\n");
  }



  if(rank == MASTER){

    for(ii=1; ii < nx+1 ; ii++){
      for(jj = 1; jj < local_ncols; jj++){
        image[ii * width + jj] = w[ii * (local_ncols + 1) + jj];
      }
    }
    MPI_Sendrecv(image, width*height, MPI_FLOAT, right, tag, image, width*height, MPI_FLOAT, left, tag,
         MPI_COMM_WORLD, &status); 
    output_image(OUTPUT_FILE, nx, ny, width, height, image);
  }

  else if(rank == size-1){
    MPI_Recv(image, width*height, MPI_FLOAT, left, tag,
         MPI_COMM_WORLD, &status);

    for(ii=1; ii < nx+1; ii++){
        for(jj= 1; jj < local_ncols; jj++){
          image[ii * width + (width-local_ncols-1) + jj] = w[ii * (local_ncols + 1) + jj] ; 
          }
        }
      MPI_Send(image, width*height, MPI_FLOAT, right, tag, MPI_COMM_WORLD);
  }
  else{
    MPI_Recv(image, width*height, MPI_FLOAT, left, tag,
         MPI_COMM_WORLD, &status);
    for(ii=1; ii < nx+1; ii++){
        for(jj= 1; jj < local_ncols + 1; jj++){
          image[(ii * width) + (rank * local_ncols-1) + jj] = w[ii * (local_ncols + 2) + jj] ; 
        }
      }
      MPI_Send(image, width*height, MPI_FLOAT, right, tag, MPI_COMM_WORLD);
  }

  
  MPI_Finalize();

  free(image);
  free(tmp_image);
  free(w);
  free(v);
  free(sendbuf);
  free(recvbuf);
  free(sendbuf_v);
  free(recvbuf_v);
}


// calculate the ncols for each rank
int calc_ncols_from_rank(int rank, int size, int NCOLS)
{
  int ncols;
  ncols = NCOLS / size;       // integer division 
  if ((NCOLS % size) != 0) {  // if there is a remainder 
    if (rank == size - 1)
      ncols += NCOLS % size;  // add remainder to last rank 
  }
  return ncols;
}



void stencil(const int nx, const int ny, const int width, const int height,
             float* restrict image, float* restrict tmp_image)
{
  for (int i = 1; i < nx + 1; ++i) {
    for (int j = 1; j < ny + 1; ++j) {
      tmp_image[j + i * width] =  image[j     + i       * width] * 3.0f / 5.0f;
      tmp_image[j + i * width] += image[j     + (i - 1) * width] * 0.5f / 5.0f;
      tmp_image[j + i * width] += image[j     + (i + 1) * width] * 0.5f / 5.0f;
      tmp_image[j + i * width] += image[j - 1 + i       * width] * 0.5f / 5.0f;
      tmp_image[j + i * width] += image[j + 1 + i       * width] * 0.5f / 5.0f;
      }
    }
}


// Create the input image
void init_image(const int nx, const int ny, const int width, const int height,
                float* image, float* tmp_image)
{
  // Zero everything
  for (int j = 0; j < ny + 2; ++j) {
    for (int i = 0; i < nx + 2; ++i) {
      image[j + i * width] = 0.0f;
      tmp_image[j + i * width] = 0.0f;
    }
  }

  const int tile_size = 64;
  // checkerboard pattern
  for (int jb = 0; jb < ny; jb += tile_size) {
    for (int ib = 0; ib < nx; ib += tile_size) {
      if ((ib + jb) % (tile_size * 2)) {
        const int jlim = (jb + tile_size > ny) ? ny : jb + tile_size;
        const int ilim = (ib + tile_size > nx) ? nx : ib + tile_size;
        for (int j = jb + 1; j < jlim + 1; ++j) {
          for (int i = ib + 1; i < ilim + 1; ++i) {
            image[j + i * width] = 100.0;
          }
        }
      }
    }
  }
}

// Routine to output the image in Netpbm grayscale binary image format
void output_image(const char* file_name, const int nx, const int ny,
                  const int width, const int height, float* image)
{
  // Open output file
  FILE* fp = fopen(file_name, "w");
  if (!fp) {
    fprintf(stderr, "Error: Could not open %s\n", OUTPUT_FILE);
    exit(EXIT_FAILURE);
  }

  // Ouptut image header
  fprintf(fp, "P5 %d %d 255\n", nx, ny);

  // Calculate maximum value of image
  // This is used to rescale the values
  // to a range of 0-255 for output
  double maximum = 0.0;
  for (int j = 1; j < ny + 1; ++j) {
    for (int i = 1; i < nx + 1; ++i) {
      if (image[j + i * width] > maximum) maximum = image[j + i * width];
    }
  }

  // Output image, converting to numbers 0-255
  for (int j = 1; j < ny + 1; ++j) {
    for (int i = 1; i < nx + 1; ++i) {
      fputc((char)(255.0f * image[j + i * width] / maximum), fp);
    }
  }

  // Close the file
  fclose(fp);
}

// Get the current time in seconds since the Epoch
double wtime(void)
{
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec + tv.tv_usec * 1e-6;
}


