#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <string.h>
#include "mpi.h"


//Define output file name
#define OUTPUT_FILE "stencil.pgm"
// Define MASTER variable
#define MASTER 0


int calc_nrows_from_rank(int rank, int size, int height);
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
  int up;              /* the rank of the process to the left */
  int down;             /* the rank of the process to the right */
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
  float *v; 



  // Check usage
  if (argc != 4) {
    fprintf(stderr, "Usage: %s nx ny niters\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  int nx = atoi(argv[1]);       /* numbers of rows */
  int ny = atoi(argv[2]);       /* numbers of columns */
  int niters = atoi(argv[3]);   /* iteration times */
  // Initiliase problem dimensions from command line arguments
  int height = nx + 2;
  int width = ny + 2; 



  // Allocate the image
  float* image = malloc(sizeof(float) * width * height);
 
  float* tmp_image = malloc(sizeof(float) * width * height);


  // Set the input image
  init_image(nx, ny, height, width, image, tmp_image);

  MPI_Init( &argc, &argv );
  MPI_Comm_size( MPI_COMM_WORLD, &size );
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  double maxtoc;
  double mintic;


  up = (rank == MASTER) ? (rank + size - 1) : (rank - 1);
  down = (rank + 1) % size;

  // we pad the outer edge of the image to avoid out of range address issues in
  // stencil


  // Allocate the spaces
  local_ncols = ny + 2;
  local_nrows = calc_nrows_from_rank(rank, size, height);
  if (local_ncols < 1) {
    fprintf(stderr,"Error: too many processes:- local_ncols < 1\n");
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }




  sendbuf = (float*)malloc(sizeof(float) * width);
  recvbuf = (float*)malloc(sizeof(float) * width);
  sendbuf_v = (float*)malloc(sizeof(float) * width);
  recvbuf_v = (float*)malloc(sizeof(float) * width);

  if(rank==MASTER || rank==(size-1)){
    w = (float*)malloc(sizeof(float*) * (local_nrows+1) *local_ncols);
    v = (float*)malloc(sizeof(float*) * (local_nrows+1) *local_ncols);
  }
  else{
    w = (float*)malloc(sizeof(float*) * (local_nrows+2) *local_ncols);
    v = (float*)malloc(sizeof(float*) * (local_nrows+2) *local_ncols);
  }


  if(rank == MASTER){
    for(ii=0;ii<local_nrows+1;ii++) {
      if(ii==local_nrows){
        for(jj=0; jj<width; jj++) {
          w[ii * width + jj] = 0.0f; 
        }
      }
      else{
        for(jj=0; jj<width; jj++) {
          w[ii * width + jj] = image[ii * width + jj];
      }
      }
      }
    }


  else if(rank == size-1){
    for(ii=0;ii<local_nrows+1;ii++) { 
        if (ii == 0){
          for(jj=0; jj<width; jj++) { 
            w[ii * width + jj] = 0.0f; 
          }
        }
        else{
          for(jj=0; jj<width; jj++) {
            w[ii * width+ jj] = image[(height-local_nrows-1+ii)*width+jj]; 
          }
        }
    }
  }

  else{
    for(ii=0;ii<local_nrows+2;ii++) { 
        if (ii == 0 || ii==local_nrows+1){
          for(jj=0; jj<width; jj++) { 
            w[ii * width + jj] = 0.0f; 
          }
        }
        else{
          for(jj=0; jj<width; jj++) { 
            w[ii * width+ jj] = image[(rank*local_nrows-1+ii)*width+jj]; 

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


  
  double tic = wtime();
  MPI_Reduce(&tic, &mintic, 1, MPI_DOUBLE, MPI_MIN, MASTER,
           MPI_COMM_WORLD);

  if(rank == MASTER){
    for (int t = 0; t < niters; ++t) {
      for(jj=0; jj<width; jj++) {
        sendbuf_v[jj] = w[(local_nrows-1)*width+jj];
      }
      MPI_Send(sendbuf_v, width, MPI_FLOAT, down, tag, MPI_COMM_WORLD);
      MPI_Recv(recvbuf, width, MPI_FLOAT, down, tag,
             MPI_COMM_WORLD, &status);



        for(jj=0; jj < width; jj++){
            w[local_nrows*width+jj] = recvbuf[jj];
        }
        stencil(local_nrows-1, ny, width,local_nrows+1,w, v);
        for(jj=0; jj<width; jj++) {
          sendbuf_v[jj] = v[(local_nrows-1)*width+jj];
        }
        MPI_Send(sendbuf_v, width, MPI_FLOAT, down, tag, MPI_COMM_WORLD);
        MPI_Recv(recvbuf, width, MPI_FLOAT, down, tag,
             MPI_COMM_WORLD, &status);



        for(jj=0; jj < width; jj++){
          v[local_nrows*width+jj] = recvbuf[jj];
      }
      stencil(local_nrows-1, ny, width,local_nrows+1,v, w);
  }
  }

  else if(rank ==size-1){
    for (int t = 0; t < niters; ++t) {
      for(jj=0; jj<width; jj++) {
        sendbuf[jj] = w[width+ jj];
      }

      MPI_Send(sendbuf, width, MPI_FLOAT, up, tag, MPI_COMM_WORLD);
      MPI_Recv(recvbuf_v, width, MPI_FLOAT, up, tag,
           MPI_COMM_WORLD, &status);


      for(jj=0; jj<width; jj++){
        w[jj] = recvbuf_v[jj];
      }

      stencil(local_nrows-1, ny, width,local_nrows+1,w, v);

      for(jj=0; jj<width; jj++) {
        sendbuf[jj] = v[width+ jj];
      }

      MPI_Send(sendbuf, width, MPI_FLOAT, up, tag, MPI_COMM_WORLD);
      MPI_Recv(recvbuf_v, width, MPI_FLOAT, up, tag,
           MPI_COMM_WORLD, &status);


      for(jj=0; jj<width; jj++){
        v[jj] = recvbuf_v[jj];
      }

      stencil(local_nrows-1, ny, width,local_nrows+1,v, w);
    }
  }

  else{
    for (int t = 0; t < niters; ++t) {
      for(jj=0; jj<width; jj++) {
        sendbuf[jj] = w[width+jj];
        sendbuf_v[jj] = w[local_nrows*width+jj];//传0那列
      }

      MPI_Sendrecv(sendbuf, width, MPI_FLOAT, up, tag,
           recvbuf, width, MPI_FLOAT, down, tag,
           MPI_COMM_WORLD, &status);
      MPI_Sendrecv(sendbuf_v, width, MPI_FLOAT, down, tag,
           recvbuf_v, width, MPI_FLOAT, up, tag,
           MPI_COMM_WORLD, &status);

      for(jj=0; jj<width; jj++){
        w[(local_nrows+1)*width+jj] = recvbuf[jj];
        w[jj] = recvbuf_v[jj];
      }

      stencil(local_nrows, ny, width,local_nrows+2,w, v);

      for(jj=0; jj<width; jj++) {
        sendbuf[jj] = v[width+jj];
        sendbuf_v[jj] = v[local_nrows*width+jj];//传0那列
      }

      MPI_Sendrecv(sendbuf, width, MPI_FLOAT, up, tag,
           recvbuf, width, MPI_FLOAT, down, tag,
           MPI_COMM_WORLD, &status);
      MPI_Sendrecv(sendbuf_v, width, MPI_FLOAT, down, tag,
           recvbuf_v, width, MPI_FLOAT, up, tag,
           MPI_COMM_WORLD, &status);

      for(jj=0; jj<width; jj++){
        v[(local_nrows+1)*width+jj] = recvbuf[jj];
        v[jj] = recvbuf_v[jj];
      }

      stencil(local_nrows, ny, width,local_nrows+2,v, w);
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

    for(ii=1; ii < local_nrows ; ii++){
      for(jj = 1; jj < width; jj++){
        image[ii * width + jj] = w[ii * width + jj];
      }
    }

    MPI_Sendrecv(image, width*height, MPI_FLOAT, down, tag, image, width*height, MPI_FLOAT, up, tag,
         MPI_COMM_WORLD, &status);
    output_image(OUTPUT_FILE, nx,ny,width,height,image);
  }

  else if(rank == size-1){
    MPI_Recv(image, width*height, MPI_FLOAT, up, tag,
         MPI_COMM_WORLD,&status);

    for(ii=1; ii < local_nrows+1; ii++){
        for(jj= 1; jj < width; jj++){
          image[(height-local_nrows-1+ii)*width+jj] = w[ii * width+ jj];
        }
      }
    MPI_Send(image, width*height, MPI_FLOAT, down, tag, MPI_COMM_WORLD);
  }

  else{
    MPI_Recv(image, width*height, MPI_FLOAT, up, tag,
         MPI_COMM_WORLD,&status);
    for(ii=1; ii < local_nrows+2; ii++){
        for(jj= 1; jj < width+1; jj++){
          image[(rank*local_nrows-1+ii)*width+jj]=w[ii * width+ jj];
        }
      }
    MPI_Send(image, width*height, MPI_FLOAT, down, tag, MPI_COMM_WORLD);
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
int calc_nrows_from_rank(int rank, int size, int height)
{
  int nrows;
  nrows = height / size;       // integer division
  if ((height % size) != 0) {  // if there is a remainder
    if (rank == size - 1)
      nrows += height % size;  // add remainder to last rank
  }
  return nrows;
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
            image[j + i * width] = 100.0f;
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
