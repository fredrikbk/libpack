#include <stdio.h>
#include <stdlib.h>

#include "mpi.h"

#include "ddtbench.h"

#define itag 0

void timing_basic_ping_pong_nelements( int DIM1, int loop, char* testname, MPI_Comm local_communicator) {

  float* array;
  int myrank;
  int base, typesize, bytes, i;
  char method[50];

  array = malloc( DIM1 * sizeof(float));

  MPI_Comm_rank( local_communicator, &myrank );

  base = myrank * DIM1 + 1;
  utilities_fill_unique_array_1D_float( &array[0], DIM1, base );

  if ( myrank == 0 ) {
    snprintf(&method[0], 50, "reference");

    MPI_Type_size( MPI_FLOAT, &typesize );
    bytes = typesize * DIM1;

    timing_init( testname, &method[0], bytes );
  }

  for( i=0 ; i<loop ; i++ ){
    if ( myrank == 0 ) {
      MPI_Send( &array[0], DIM1, MPI_FLOAT, 1, itag, local_communicator );
      MPI_Recv( &array[0], DIM1, MPI_FLOAT, 1, itag, local_communicator, MPI_STATUS_IGNORE );
      timing_record(3);
    } else {
      MPI_Recv( &array[0], DIM1, MPI_FLOAT, 0, itag, local_communicator, MPI_STATUS_IGNORE );
      MPI_Send( &array[0], DIM1, MPI_FLOAT, 0, itag, local_communicator );
    }
  }

  if ( myrank == 0 ) {
    timing_print( 1 );
  }

  free(array);
}

void timing_basic_alltoall_nelements( int DIM1, int procs, int loop, char* testname, MPI_Comm local_communicator) {
      
  float* send_array;
  float* recv_array;
  
  int myrank;
  int base, typesize, bytes, i;
  char method[50];

  send_array = malloc( DIM1 * procs * sizeof(float));
  recv_array = malloc( DIM1 * procs * sizeof(float));

  MPI_Comm_rank( local_communicator, &myrank );

  base = myrank * DIM1 + 1;
  utilities_fill_unique_array_1D_float( &send_array[0], DIM1, base );

  if ( myrank == 0 ) {
    snprintf(method, 50, "reference");
        
    MPI_Type_size( MPI_FLOAT, &typesize );
    bytes = typesize * DIM1 * procs;

    timing_init( testname, &method[0], bytes );
  }

  for( i=0 ; i<loop ; i++ ) {
    MPI_Alltoall(&send_array[0], DIM1, MPI_FLOAT, &recv_array[0], DIM1, MPI_FLOAT, local_communicator );
    MPI_Alltoall(&recv_array[0], DIM1, MPI_FLOAT, &send_array[0], DIM1, MPI_FLOAT, local_communicator );
    if ( myrank == 0 ) {
      timing_record(3);
    }
  }

  if ( myrank == 0 ) {
    timing_print( 1 );
  }

  free(send_array);
  free(recv_array);
}     
