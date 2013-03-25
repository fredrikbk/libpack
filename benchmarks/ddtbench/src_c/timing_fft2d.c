#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mpi.h"

#include "ddtbench.h"
      
void timing_fft2d_ddt( int DIM1, int procs, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator ) {

  double* matrix;
  double* recv_array;

  int myrank;
  int i, j , base, typesize, bytes;

//! variables for the datatype construction
  MPI_Datatype dtype_vector_t, dtype_resize_t, dtype_scatter_t, dtype_gather_t, dtype_complex_t;
  MPI_Aint lb, extent;

  char method[50];

//! just some statements to prevent compiler warnings of unused variables
//! those parameter are included for future features
  *correct_flag = 0;
  *ptypesize = 0;
//  typesize = filehandle_debug

  matrix = malloc( DIM1 * DIM1/procs * 2 * sizeof(double) );
  recv_array = malloc( DIM1 * DIM1/procs * 2 * sizeof(double) );

  MPI_Comm_rank( local_communicator, &myrank );

//! ================= initialize the arrays =================

  base = myrank * DIM1 * DIM1/procs * 2 + 1;
  utilities_fill_unique_array_2D_double( &matrix[0], 2*DIM1, DIM1/procs, base );
 
  if ( myrank == 0 ) {
    snprintf( &method[0], 50, "mpi_ddt" );

    MPI_Type_size( MPI_DOUBLE, &typesize );
    bytes = DIM1/procs * DIM1 * 2 * typesize;
 
    timing_init( testname, &method[0], bytes );
  }
     
  for( i=0 ; i<outer_loop ; i++ ) {

    MPI_Type_contiguous( 2, MPI_DOUBLE, &dtype_complex_t );
    MPI_Type_vector( DIM1/procs, 1, DIM1, dtype_complex_t, &dtype_vector_t );
    lb = 0;
    extent = 2 * sizeof(double);
    MPI_Type_create_resized( dtype_vector_t, lb, extent, &dtype_resize_t );
    MPI_Type_contiguous( DIM1/procs, dtype_resize_t, &dtype_scatter_t );
    MPI_Type_commit( &dtype_scatter_t  );

    MPI_Type_free( &dtype_vector_t );
    MPI_Type_free( &dtype_resize_t );

    MPI_Type_vector( DIM1/procs, DIM1/procs, DIM1, dtype_complex_t, &dtype_vector_t );
    lb = 0;
    extent = DIM1/procs * 2 * sizeof(double);
    MPI_Type_create_resized( dtype_vector_t, lb, extent, &dtype_gather_t );
    MPI_Type_commit( &dtype_gather_t );

    MPI_Type_free( &dtype_vector_t );
    MPI_Type_free( &dtype_complex_t );

    if ( myrank == 0 ) {
      timing_record(1);
    }

    for( j=0 ; j<inner_loop ; j++ ) {

      MPI_Alltoall( &matrix[0], 1, dtype_gather_t, &recv_array[0], 1, dtype_scatter_t, local_communicator );
      MPI_Alltoall( &recv_array[0], 1, dtype_gather_t, &matrix[0], 1, dtype_scatter_t, local_communicator );
      if ( myrank == 0 ) {
        timing_record(3);
      }

    }  //! inner loop

    MPI_Type_free( &dtype_gather_t );
    MPI_Type_free( &dtype_scatter_t );

    if ( myrank == 0 ) {
      timing_record(5);
    }

  }  //! outer loop

  if ( myrank == 0 ) {
    timing_print( 1 );
  }

  free(recv_array);
  free(matrix);
}

void timing_fft2d_manual( int DIM1, int procs, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator ) {

  double* matrix;
  double* recv_buffer;
  double* buffer;

  int myrank;
  int i, j, k, l, base, typesize, bytes;

  char method[50];

//! just some statements to prevent compiler warnings of unused variables
//! those parameter are included for future features
  *correct_flag = 0;
  *ptypesize = 0;
//  typesize = filehandle_debug;

  matrix = malloc( 2 * DIM1 * DIM1/procs * sizeof(double) );
  recv_buffer = malloc( 2 * DIM1 * DIM1/procs * sizeof(double) );

  MPI_Comm_rank( local_communicator, &myrank );

//! ================= initialize the arrays =================

  base = myrank * DIM1 * DIM1/procs * 2 + 1;
  utilities_fill_unique_array_2D_double( &matrix[0], 2*DIM1, DIM1/procs, base );

  if ( myrank == 0 ) {
    snprintf( &method[0], 50, "manual" );

    MPI_Type_size( MPI_DOUBLE, &typesize );
    bytes = 2 * DIM1/procs * DIM1 * typesize;

    timing_init( testname, &method[0], bytes );
  }
     
  for( i=0 ; i<outer_loop ; i++ ) {

    buffer = malloc( 2 * DIM1 * DIM1/procs * sizeof(double) );

    if ( myrank == 0 ) {
      timing_record(1);
    }

    for( j=0 ; j<inner_loop ; j++ ) {
//! pack the data
      for( k=0 ; k < procs ; k++ ) {
        for( l=0 ; l<DIM1/procs ; l++ ) {
          memcpy( &buffer[k * DIM1/procs * DIM1/procs + l * DIM1/procs], &matrix[k * DIM1/procs + l * DIM1], 2 * DIM1/procs * sizeof(double) );
        }
      }

      if ( myrank == 0 ) {
        timing_record(2);
      }

      MPI_Alltoall( &buffer[0], 2*DIM1/procs*DIM1/procs, MPI_DOUBLE, &recv_buffer[0], 2*DIM1/procs*DIM1/procs, MPI_DOUBLE, local_communicator );

//! unpack the data
      for( k=0 ; k<DIM1/procs ; k++ ) {
        for( l=0 ; l<DIM1 ; l++ ) {
          matrix[l+k*DIM1] = recv_buffer[k + l * DIM1/procs];
        }
      }

//! pack the data            
      for( k=0 ; k < procs ; k++ ) {
        for( l=0 ; l<DIM1/procs ; l++ ) {
          memcpy( &buffer[k * DIM1/procs * DIM1/procs + l * DIM1/procs], &matrix[k * DIM1/procs + l * DIM1], 2 * DIM1/procs * sizeof(double) );
        }
      }

      MPI_Alltoall( &buffer[0], 2*DIM1/procs*DIM1/procs, MPI_DOUBLE, recv_buffer, 2*DIM1/procs*DIM1/procs, MPI_DOUBLE, local_communicator );

      if ( myrank == 0 ) {
        timing_record(3);
      }

//! unpack the data
      for( k=0 ; k<DIM1/procs ; k++ ) {
        for( l=0 ; l<DIM1 ; l++ ) {
          matrix[l+k*DIM1] = recv_buffer[k + l * DIM1/procs];
        }
      }

      if ( myrank == 0 ) {
        timing_record(4);
      }

    } //! inner loop

    free( buffer );

    if ( myrank == 0 ) {
      timing_record(5);
    }

  } //! outer loop

  if ( myrank == 0 ) {
    timing_print( 1 );
  }

  free(recv_buffer);
  free(matrix);
}

void timing_fft2d_mpi_pack_ddt( int DIM1, int procs, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator ) {

  double* matrix;
  double* recv_buffer;
  double* buffer;      

  int myrank;
  int i, j, base, typesize, bytes, pos;

//! variables for the datatype construction
  MPI_Datatype dtype_vector_t, dtype_resize_t, dtype_scatter_t, dtype_gather_t, dtype_complex_t;
  MPI_Aint lb, extent;

  char method[50];

//! just some statements to prevent compiler warnings of unused variables
//! those parameter are included for future features
  *correct_flag = 0;
  *ptypesize = 0;
//  typesize = filehandle_debug

  matrix = malloc( DIM1 * DIM1/procs * 2 * sizeof(double) );
  recv_buffer = malloc( DIM1 * DIM1/procs * 2 * sizeof(double) );

  MPI_Comm_rank( local_communicator, &myrank );

//! ================= initialize the arrays =================

  base = myrank * DIM1 * DIM1/procs * 2 + 1;
  utilities_fill_unique_array_2D_double( &matrix[0], 2*DIM1, DIM1/procs, base );
 
  if ( myrank == 0 ) {
    snprintf( &method[0], 50, "mpi_pack_ddt" );

    MPI_Type_size( MPI_DOUBLE, &typesize );
    bytes = DIM1/procs * DIM1 * 2 * typesize;

    timing_init( testname, &method[0], bytes );
  }
     
  for( i=0 ; i<outer_loop ; i++ ) {

    buffer = malloc( 2* DIM1 * DIM1/procs * sizeof(double) );
    MPI_Type_size( MPI_DOUBLE, &typesize );
    bytes = 2 * DIM1/procs * DIM1/procs * typesize;

    MPI_Type_contiguous( 2, MPI_DOUBLE, &dtype_complex_t );
    MPI_Type_vector( DIM1/procs, 1, DIM1, dtype_complex_t, &dtype_vector_t );
    lb = 0;
    extent = 2 * sizeof(double);
    MPI_Type_create_resized( dtype_vector_t, lb, extent, &dtype_resize_t );
    MPI_Type_contiguous( DIM1/procs, dtype_resize_t, &dtype_scatter_t );
    MPI_Type_commit( &dtype_scatter_t );

    MPI_Type_free( &dtype_vector_t );
    MPI_Type_free( &dtype_resize_t );

    MPI_Type_vector( DIM1/procs, DIM1/procs, DIM1, dtype_complex_t, &dtype_vector_t );
    lb = 0;
    extent = DIM1/procs * 2 * sizeof(double);
    MPI_Type_create_resized( dtype_vector_t, lb, extent, &dtype_gather_t );
    MPI_Type_commit( &dtype_gather_t );

    MPI_Type_free( &dtype_vector_t );
    MPI_Type_free( &dtype_complex_t );

    if ( myrank == 0 ) {
      timing_record(1);
    }

    for( j=0 ; j<inner_loop ; j++ ) {
//! pack the data
      pos = 0;
      MPI_Pack( &matrix[0], 1, dtype_gather_t, &buffer[0], bytes, &pos, local_communicator );
      if ( myrank == 0 ) {
        timing_record(2);
      }

      MPI_Alltoall( &buffer[0], bytes, MPI_PACKED, &recv_buffer[0], bytes, MPI_PACKED, local_communicator );

//! unpack the data
      pos = 0;
      MPI_Unpack( &recv_buffer[0], bytes, &pos, &matrix[0], 1, dtype_scatter_t, local_communicator );

//! pack the data
      pos = 0;
      MPI_Pack( &matrix[0], 1, dtype_gather_t, &buffer[0], bytes, &pos, local_communicator );

      MPI_Alltoall( &buffer[0], bytes, MPI_PACKED, &recv_buffer[0], bytes, MPI_PACKED, local_communicator );

      if ( myrank == 0 ) {
        timing_record(3);
      }
//! unpack the data
      pos = 0;
      MPI_Unpack( &recv_buffer[0], bytes, &pos, &matrix[0], 1, dtype_scatter_t, local_communicator );
      if ( myrank == 0 ) {
        timing_record(4);
      }

   } //! inner loop

   MPI_Type_free( &dtype_gather_t );
   MPI_Type_free( &dtype_scatter_t );

   free( buffer );

   if ( myrank == 0 ) {
     timing_record(5);
   }

 } //! outer loop

 if ( myrank == 0 ) {
   timing_print( 1 );
 }

 free(recv_buffer);
 free(matrix);
}
