#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mpi.h"

#include "ddtbench.h"

#define itag 0

static inline int idx2D(int x, int y, int DIM1) {
  return x+y*DIM1;
}

void timing_specfem3D_oc_ddt( int DIM1, int icount, int* list, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator ) {
      
  float *array;

  int i, j, bytes, base, typesize;
  int myrank;

//! variables for the MPI derived datatypes
  int* displacement;
  MPI_Datatype dtype_indexed_t;

  char method[50];

//! just some statements to prevent compiler warnings of unused variables
//! those parameter are included for future features
  *correct_flag = 0;
  *ptypesize = 0;
//  typesize = filehandle_debug

//conversion from fortran to c
  displacement = malloc( icount * outer_loop * sizeof(float) );
  for ( i=0 ; i<outer_loop ; i++ ) {
    for( j=0 ; j < icount ; j++ ) {
      displacement[idx2D(j,i,icount)] = list[idx2D(j,i,icount)] - 1;
    }
  }

  array = malloc( DIM1 * sizeof(float) );

  MPI_Comm_rank( local_communicator, &myrank );

//! ================= initialize the array ==================

  base = myrank * DIM1 + 1;
  utilities_fill_unique_array_1D_float( &array[0], DIM1, base);

  if ( myrank == 0 ) {
    snprintf(method, 50, "mpi_ddt");

    MPI_Type_size( MPI_FLOAT, &typesize );
    bytes = icount * typesize;

    timing_init( testname, &method[0], bytes );
  }

  for( i=0 ; i<outer_loop ; i++ ) {

//! ========== building the MPI derived datatype ============

    MPI_Type_create_indexed_block(icount, 1, &displacement[idx2D(0,i,icount)], MPI_FLOAT, &dtype_indexed_t );
    MPI_Type_commit( &dtype_indexed_t );

    if ( myrank == 0 ) {
      timing_record(1);
    }

    for( j=0 ; j<inner_loop ; j++ ) {
//! =============== ping pong communication =================

      if ( myrank == 0 ) {
//! send the data from rank 0 to rank 1          
        MPI_Send( &array[0], 1, dtype_indexed_t, 1, itag, local_communicator );
//! receive the data from rank 1 back
        MPI_Recv( &array[0], 1, dtype_indexed_t, 1, itag, local_communicator, MPI_STATUS_IGNORE );
        timing_record(3);
//! now for rank 1
      } else {
        MPI_Recv( &array[0], 1, dtype_indexed_t, 0, itag, local_communicator, MPI_STATUS_IGNORE );
//! send the received data back
        MPI_Send( &array[0], 1, dtype_indexed_t, 0, itag, local_communicator );
      }

    } //! inner loop

//! ======================= clean up ========================

    MPI_Type_free( &dtype_indexed_t );

    if ( myrank == 0 ) {
      timing_record(5);
    }

  } //! outer loop

  if ( myrank == 0 ) {
    timing_print( 1 );
  }

  free(array);
  free(displacement);
}

void timing_specfem3D_oc_manual( int DIM1, int icount, int* list, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator ) {

  float* array;
  float* buffer;

  int i, j, k, bytes, base, typesize;
  int myrank;

  int* displacement;

  char method[50];

//! just some statements to prevent compiler warnings of unused variables
//! those parameter are included for future features
  *correct_flag = 0;
  *ptypesize = 0;
//  typesize = filehandle_debug

//conversion from fortran to c
  displacement = malloc( icount * outer_loop * sizeof(float) );
  for ( i=0 ; i<outer_loop ; i++ ) {
    for( j=0 ; j < icount ; j++ ) {
      displacement[idx2D(j,i,icount)] = list[idx2D(j,i,icount)] - 1;
    }
  }

  array = malloc( DIM1 * sizeof(float) );

  MPI_Comm_rank( local_communicator, &myrank );

//! ================= initialize the array ==================

  base = myrank * DIM1 + 1;
  utilities_fill_unique_array_1D_float( array, DIM1, base);

  if ( myrank == 0 ) {
    snprintf(method, 50, "manual");

    MPI_Type_size( MPI_FLOAT, &typesize );
    bytes = icount * typesize;

    timing_init( testname, &method[0], bytes );
  }

  for( i=0 ; i<outer_loop ; i++ ) {

    buffer = malloc( icount *sizeof(float) );

    if ( myrank == 0 ) {
      timing_record(1);
    }

    for( j=0 ; j<inner_loop ; j++ ) {
//! =============== ping pong communication =================

      if ( myrank == 0 ) {
//! pack of the data
        for( k=0 ; k < icount ; k++ ) {
          buffer[k] = array[displacement[idx2D(k,i,icount)]];
        }
        timing_record(2);
//! send the data from rank 0 to rank 1          
        MPI_Send( buffer, icount, MPI_FLOAT, 1, itag, local_communicator );
//! receive the data from rank 1 back
        MPI_Recv( buffer, icount, MPI_FLOAT, 1, itag, local_communicator, MPI_STATUS_IGNORE );
        timing_record(3);
//! unpack of the data
        for( k=1 ; k<icount ; k++ ) {
          array[displacement[idx2D(k,i,icount)]] = buffer[k];
        }
        timing_record(4);
//! now for rank 1
      } else {
        MPI_Recv( buffer, icount, MPI_FLOAT, 0, itag, local_communicator, MPI_STATUS_IGNORE );
        for( k=0 ; k < icount ; k++ ) {
          array[displacement[idx2D(k,i,icount)]] = buffer[k];
        }
        for( k=0 ; k < icount ; k++ ) {
          buffer[k] = array[displacement[idx2D(k,i,icount)]];
        }
//! send the received data back
        MPI_Send( buffer, icount, MPI_FLOAT, 0, itag, local_communicator );
      }

    } //! inner loop

//! ======================= clean up ========================

    free( buffer );

    if ( myrank == 0 ) {
      timing_record(5);
    }

  } //! outer loop

  if ( myrank == 0 ) {
    timing_print( 1 );
  }

  free(array);
  free(displacement);
}

void timing_specfem3D_oc_mpi_pack_ddt( int DIM1, int icount, int* list, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator ) {

  float* array;
  float* buffer;

  int i, j, bytes, base, typesize, pos;
  int myrank;

//! variables for the MPI derived datatypes
  int* displacement;
  MPI_Datatype dtype_indexed_t;

  char method[50];

//! just some statements to prevent compiler warnings of unused variables
//! those parameter are included for future features
  correct_flag = 0;
  ptypesize = 0;
//  typesize = filehandle_debug

//conversion from fortran to c
  displacement = malloc( icount * outer_loop * sizeof(float) );
  for ( i=0 ; i<outer_loop ; i++ ) {
    for( j=0 ; j < icount ; j++ ) {
      displacement[idx2D(j,i,icount)] = list[idx2D(j,i,icount)] - 1;
    }
  }

  array = malloc( DIM1 * sizeof(float) );

  MPI_Comm_rank( local_communicator, &myrank );

//! ================= initialize the array ==================

  base = myrank * DIM1 + 1;
  utilities_fill_unique_array_1D_float( &array[0], DIM1, base);

  if ( myrank == 0 ) {
    snprintf( method, 50, "mpi_pack_ddt" );

    MPI_Type_size( MPI_FLOAT, &typesize );
    bytes = icount * typesize;

    timing_init( testname, &method[0], bytes );
  }

  for( i=0 ; i<outer_loop ; i++ ) {

    buffer = malloc( icount * sizeof(float) );

    MPI_Type_size( MPI_FLOAT, &typesize );
    bytes = icount * typesize;

//! ========== building the MPI derived datatype ============
    MPI_Type_create_indexed_block(icount, 1, &displacement[idx2D(0,i,icount)], MPI_FLOAT, &dtype_indexed_t );
    MPI_Type_commit( &dtype_indexed_t );

    if ( myrank == 0 ) {
      timing_record(1);
    }

    for( j=0 ; j<inner_loop ; j++ ) {
//! =============== ping pong communication =================

      if ( myrank == 0 ) {
//! pack of the data
        pos = 0;
        MPI_Pack( &array[0], 1, dtype_indexed_t, &buffer[0], bytes, &pos, local_communicator );
        timing_record(2);
//! send the data from rank 0 to rank 1          
        MPI_Send( &buffer[0], pos, MPI_PACKED, 1, itag, local_communicator );
//! receive the data from rank 1 back
        MPI_Recv( &buffer[0], bytes, MPI_PACKED, 1, itag, local_communicator, MPI_STATUS_IGNORE );
        timing_record(3);
//! unpack of the data
        pos = 0;
        MPI_Unpack( &buffer[0], bytes, &pos, &array[0], 1, dtype_indexed_t, local_communicator );
        timing_record(4);
//! now for rank 1
      } else {
        MPI_Recv( &buffer[0], bytes, MPI_PACKED, 0, itag, local_communicator, MPI_STATUS_IGNORE );
        pos = 0;
        MPI_Unpack( &buffer[0], bytes, &pos, &array[0], 1, dtype_indexed_t, local_communicator );
        pos = 0;
        MPI_Pack( &array[0], 1, dtype_indexed_t, &buffer[0], bytes, &pos, local_communicator );
//! send the received data back
        MPI_Send( &buffer[0], pos, MPI_PACKED, 0, itag, local_communicator );
      }

    } //! inner loop

//! ======================= clean up ========================

    free( buffer );

    MPI_Type_free( &dtype_indexed_t );

    if ( myrank == 0 ) {
      timing_record(5);
    }

  } //! outer loop

  if ( myrank == 0 ) {
    timing_print( 1 );
  }

  free(array);
  free(displacement);
}

void timing_specfem3D_cm_ddt( int DIM2_cm, int DIM2_ic, int icount_cm, int icount_ic, int* list_cm, int* list_ic, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator ) {
      
  float* array_cm;
  float* array_ic;
      
  int i, j, base, bytes, typesize;
  int myrank, maximum;

  char method[50];

//! variables for the MPI derived datatypes
  MPI_Aint struct_displacement[2];
  MPI_Datatype dtype_temp_t[2], dtype_indexed_t;
  int blocklength[2];

  int* temp_displacement_ic;
  int* temp_displacement_cm;
  int* displacement;

//! just some statements to prevent compiler warnings of unused variables
//! those parameter are included for future features
  *correct_flag = 0;
  *ptypesize = 0;
 //typesize = filehandle_debug
 
  array_cm = malloc( 3 * DIM2_cm * sizeof(float) );
  array_ic = malloc( 3 * DIM2_ic * sizeof(float) );

// fortran to c conversion
  temp_displacement_cm = malloc( icount_cm * outer_loop * sizeof(int) );
  temp_displacement_ic = malloc( icount_ic * outer_loop * sizeof(int) );
  for( i=0 ; i<outer_loop ; i++ ) {
    for( j=0 ; j<icount_cm ; j++ ) {
      temp_displacement_cm[idx2D(j,i,icount_cm)] = list_cm[idx2D(j,i,icount_cm)] - 1;
    }
    for( j=0 ; j<icount_ic ; j++ ) {
      temp_displacement_ic[idx2D(j,i,icount_ic)] = list_ic[idx2D(j,i,icount_ic)] - 1;
    }
  }

  MPI_Comm_rank( local_communicator, &myrank );

//! ================= initialize the arrays =================
  base = myrank * (DIM2_cm+DIM2_ic) * 3 + 1;
  utilities_fill_unique_array_2D_float( &array_cm[0], 3, DIM2_cm, base);
  base = base + DIM2_cm * 3;
  utilities_fill_unique_array_2D_float( &array_ic[0], 3, DIM2_ic, base);

  if ( myrank == 0 ) {
    snprintf( method, 50, "mpi_ddt" );

    MPI_Type_size( MPI_FLOAT, &typesize );
    bytes = (icount_cm + icount_ic) * 3 * typesize;
 
    timing_init( testname, &method[0], bytes );
  }

  for( i=0 ; i<outer_loop ; i++ ) {
//! ========== building the MPI derived datatype ============
    if (icount_cm < icount_ic) {
      maximum = icount_ic;
    } else {
      maximum = icount_cm;
    }
    displacement = malloc( maximum * sizeof(int) );

    MPI_Get_address( &array_cm[0], &struct_displacement[0] );
    MPI_Get_address( &array_ic[0], &struct_displacement[1] );
 
    blocklength[0] = 1;
    blocklength[1] = 1;

    for( j = 0 ; j < icount_cm ; j++ ) {
      displacement[j] = temp_displacement_cm[idx2D(j,i,icount_cm)] * 3;
    }
    MPI_Type_create_indexed_block( icount_cm, 3, &displacement[0], MPI_FLOAT, &dtype_temp_t[0] );
    for( j = 0 ; j < icount_ic ; j++ ) {
      displacement[j] = temp_displacement_ic[idx2D(j,i,icount_ic)] * 3;
    }
    MPI_Type_create_indexed_block( icount_ic, 3, &displacement[0], MPI_FLOAT, &dtype_temp_t[1] );

    MPI_Type_create_struct( 2, &blocklength[0], &struct_displacement[0], &dtype_temp_t[0], &dtype_indexed_t );
    MPI_Type_commit( &dtype_indexed_t );
    MPI_Type_free( &dtype_temp_t[0] );
    MPI_Type_free( &dtype_temp_t[1] );

    free( displacement );

    if ( myrank == 0 ) {
      timing_record(1);
    }

//! =============== ping pong communication =================

    for( j=0 ; j<inner_loop ; j++ ) {
      if ( myrank == 0 ) {
//! send the data from rank 0 to rank 1          
        MPI_Send( MPI_BOTTOM, 1, dtype_indexed_t, 1, itag, local_communicator );
//! receive the data from rank 1 back
        MPI_Recv( MPI_BOTTOM, 1, dtype_indexed_t, 1, itag, local_communicator, MPI_STATUS_IGNORE );
        timing_record(3);
//! now for rank 1
      } else {
//! receive the data from rank 0
        MPI_Recv( MPI_BOTTOM, 1, dtype_indexed_t, 0, itag, local_communicator, MPI_STATUS_IGNORE );
//! send the data back to rank 0
        MPI_Send( MPI_BOTTOM, 1, dtype_indexed_t, 0, itag, local_communicator );
      }
    }//! inner loop

//! ======================= clean up ========================

    MPI_Type_free( &dtype_indexed_t );

    if ( myrank == 0 ) {
      timing_record(5);
    }

  } //! outer loop

  if ( myrank == 0 ) {
    timing_print( 1 );
  }

  free(array_cm);
  free(array_ic);
  free(temp_displacement_cm);
  free(temp_displacement_ic);
}

void timing_specfem3D_cm_manual( int DIM2_cm, int DIM2_ic, int icount_cm, int icount_ic, int* list_cm, int* list_ic, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator ) {

  float* array_cm;
  float* array_ic;

  float* buffer;

  int i, j, base, bytes, typesize, counter, isize, k;
  int myrank;

  int* temp_displacement_ic;
  int* temp_displacement_cm;

  char method[50];

//! just some statements to prevent compiler warnings of unused variables
//! those parameter are included for future features
  *correct_flag = 0;
  *ptypesize = 0;
//  typesize = filehandle_debug

  array_cm = malloc( 3 * DIM2_cm * sizeof(float) );
  array_ic = malloc( 3 * DIM2_ic * sizeof(float) );

// fortran to c conversion
  temp_displacement_cm = malloc( icount_cm * outer_loop * sizeof(int) );
  temp_displacement_ic = malloc( icount_ic * outer_loop * sizeof(int) );
  for( i=0 ; i<outer_loop ; i++ ) {
    for( j=0 ; j<icount_cm ; j++ ) {
      temp_displacement_cm[idx2D(j,i,icount_cm)] = list_cm[idx2D(j,i,icount_cm)] - 1;
    }
    for( j=0 ; j<icount_ic ; j++ ) {
      temp_displacement_ic[idx2D(j,i,icount_ic)] = list_ic[idx2D(j,i,icount_ic)] - 1;
    }
  }

  MPI_Comm_rank( local_communicator, &myrank );

//! ================= initialize the arrays =================
  base = myrank * (DIM2_cm+DIM2_ic) * 3 + 1;
  utilities_fill_unique_array_2D_float( array_cm, 3, DIM2_cm, base );
  base = base + DIM2_cm * 3;
  utilities_fill_unique_array_2D_float( array_ic, 3, DIM2_ic, base );

  if ( myrank == 0 ) {
    snprintf(method, 50, "manual");

    MPI_Type_size( MPI_FLOAT, &typesize );
    bytes = (icount_cm + icount_ic) * 3 * typesize;

    timing_init( testname, &method[0], bytes );
  }

  for( i=0 ; i<outer_loop ; i++ ) {
    isize = (icount_cm+icount_ic) * 3;
    buffer = malloc( isize * sizeof(float) );
        
    if ( myrank == 0 ) {
      timing_record(1);
    }

//! =============== ping pong communication =================

    for( j=0 ; j<inner_loop ; j++ ) {
      if ( myrank == 0 ) {
//! pack the buffer
        counter = 0;
        for( k=0 ; k<icount_cm ; k++ ) {
          buffer[counter++] = array_cm[idx2D(0,temp_displacement_cm[idx2D(k,i,icount_cm)],3)];
          buffer[counter++] = array_cm[idx2D(1,temp_displacement_cm[idx2D(k,i,icount_cm)],3)];
          buffer[counter++] = array_cm[idx2D(2,temp_displacement_cm[idx2D(k,i,icount_cm)],3)];
        }
        for( k=0 ; k<icount_ic ; k++ ) { 
          buffer[counter++] = array_ic[idx2D(0,temp_displacement_ic[idx2D(k,i,icount_ic)],3)];
          buffer[counter++] = array_ic[idx2D(1,temp_displacement_ic[idx2D(k,i,icount_ic)],3)];
          buffer[counter++] = array_ic[idx2D(2,temp_displacement_ic[idx2D(k,i,icount_ic)],3)];
        }
        timing_record(2);
//! send the data from rank 0 to rank 1          
        MPI_Send( &buffer[0], isize, MPI_FLOAT, 1, itag, local_communicator );
//! receive the data from rank 1 back
        MPI_Recv( &buffer[0], isize, MPI_FLOAT, 1, itag, local_communicator, MPI_STATUS_IGNORE );
        timing_record(3);
//! unpack the data
        counter = 0;
        for( k=0 ; k<icount_cm ; k++ ) {
          array_cm[idx2D(0,temp_displacement_cm[idx2D(k,i,icount_cm)],3)] = buffer[counter++];
          array_cm[idx2D(1,temp_displacement_cm[idx2D(k,i,icount_cm)],3)] = buffer[counter++];
          array_cm[idx2D(2,temp_displacement_cm[idx2D(k,i,icount_cm)],3)] = buffer[counter++];
        }
        for( k=0 ; k<icount_ic ; k++ ) {
          array_ic[idx2D(0,temp_displacement_ic[idx2D(k,i,icount_ic)],3)] = buffer[counter++];
          array_ic[idx2D(1,temp_displacement_ic[idx2D(k,i,icount_ic)],3)] = buffer[counter++];
          array_ic[idx2D(2,temp_displacement_ic[idx2D(k,i,icount_ic)],3)] = buffer[counter++];
        }
        timing_record(4);
//! now for rank 1
      } else {
//! receive the data from rank 0
        MPI_Recv( &buffer[0], isize, MPI_FLOAT, 0, itag, local_communicator, MPI_STATUS_IGNORE );
//! unpack the data
        counter = 0;
        for( k=0 ; k<icount_cm ; k++ ) {
          array_cm[idx2D(0,temp_displacement_cm[idx2D(k,i,icount_cm)],3)] = buffer[counter++];
          array_cm[idx2D(1,temp_displacement_cm[idx2D(k,i,icount_cm)],3)] = buffer[counter++];
          array_cm[idx2D(2,temp_displacement_cm[idx2D(k,i,icount_cm)],3)] = buffer[counter++];
        }
        for( k=0 ; k<icount_ic ; k++ ) {
          array_ic[idx2D(0,temp_displacement_ic[idx2D(k,i,icount_ic)],3)] = buffer[counter++];
          array_ic[idx2D(1,temp_displacement_ic[idx2D(k,i,icount_ic)],3)] = buffer[counter++];
          array_ic[idx2D(2,temp_displacement_ic[idx2D(k,i,icount_ic)],3)] = buffer[counter++];
        }
//! pack the data
        counter = 0;
        for( k=0 ; k<icount_cm ; k++ ) {
          buffer[counter++] = array_cm[idx2D(0,temp_displacement_cm[idx2D(k,i,icount_cm)],3)];
          buffer[counter++] = array_cm[idx2D(1,temp_displacement_cm[idx2D(k,i,icount_cm)],3)];
          buffer[counter++] = array_cm[idx2D(2,temp_displacement_cm[idx2D(k,i,icount_cm)],3)];
        }
        for( k=0 ; k<icount_ic ; k++ ) { 
          buffer[counter++] = array_ic[idx2D(0,temp_displacement_ic[idx2D(k,i,icount_ic)],3)];
          buffer[counter++] = array_ic[idx2D(1,temp_displacement_ic[idx2D(k,i,icount_ic)],3)];
          buffer[counter++] = array_ic[idx2D(2,temp_displacement_ic[idx2D(k,i,icount_ic)],3)];
        }
//! send the data back to rank 0
        MPI_Send( buffer, isize, MPI_FLOAT, 0, itag, local_communicator );
      }
    } //! inner loop

//! ======================= clean up ========================
    free( buffer );

    if ( myrank == 0 ) {
      timing_record(5);
    }

  } //! outer loop

  if ( myrank == 0 ) {
    timing_print( 1 );
  }

  free(temp_displacement_ic);
  free(temp_displacement_cm);
  free(array_ic);
  free(array_cm);
}

void timing_specfem3D_cm_mpi_pack_ddt( int DIM2_cm, int DIM2_ic, int icount_cm, int icount_ic, int* list_cm, int* list_ic, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator ) {

  float* array_cm;
  float* array_ic;

  float* buffer;

  int i, j, base, bytes, typesize, isize;
  int myrank, maximum, pos;

  char method[50];

//! variables for the MPI derived datatypes
  MPI_Aint struct_displacement[2];
  MPI_Datatype dtype_temp_t[2], dtype_indexed_t;
  int blocklength[2];

  int* displacement;
  int* temp_displacement_ic;
  int* temp_displacement_cm;

//! just some statements to prevent compiler warnings of unused variables
//! those parameter are included for future features
  *correct_flag = 0;
  *ptypesize = 0;
//  typesize = filehandle_debug

  array_cm = malloc( 3 * DIM2_cm * sizeof(float) );
  array_ic = malloc( 3 * DIM2_ic * sizeof(float) );

// fortran to c conversion
  temp_displacement_cm = malloc( icount_cm * outer_loop * sizeof(int) );
  temp_displacement_ic = malloc( icount_ic * outer_loop * sizeof(int) );
  for( i=0 ; i<outer_loop ; i++ ) {
    for( j=0 ; j<icount_cm ; j++ ) {
      temp_displacement_cm[idx2D(j,i,icount_cm)] = list_cm[idx2D(j,i,icount_cm)] - 1;
    }
    for( j=0 ; j<icount_ic ; j++ ) {
      temp_displacement_ic[idx2D(j,i,icount_ic)] = list_ic[idx2D(j,i,icount_ic)] - 1;
    }
  }

  MPI_Comm_rank( local_communicator, &myrank );

//! ================= initialize the arrays =================
  base = myrank * (DIM2_cm+DIM2_ic) * 3 + 1;
  utilities_fill_unique_array_2D_float( &array_cm[0], 3, DIM2_cm, base );
  base = base + DIM2_cm * 3;
  utilities_fill_unique_array_2D_float( &array_ic[0], 3, DIM2_ic, base );

  if ( myrank == 0 ) {
    snprintf(method, 50, "mpi_pack_ddt" );

    MPI_Type_size( MPI_FLOAT, &typesize );
    bytes = (icount_cm + icount_ic) * 3 * typesize;

    timing_init( testname, &method[0], bytes );
  }

  for( i=0 ; i < outer_loop ; i++ ) {
    isize = (DIM2_cm+DIM2_ic) * 3;
    MPI_Type_size( MPI_FLOAT, &typesize );
    bytes = isize * typesize;

    buffer = malloc( isize * sizeof(float) );

//! ========== building the MPI derived datatype ============
    if (icount_cm < icount_ic) {
      maximum = icount_ic;
    } else {
      maximum = icount_cm;
    }

    displacement = malloc( maximum * sizeof(int) );

    MPI_Get_address( &array_cm[0], &struct_displacement[0] );
    MPI_Get_address( &array_ic[0], &struct_displacement[1] );

    blocklength[0] = 1;
    blocklength[1] = 1;

    for( j = 0 ; j < icount_cm ; j++ ) {
      displacement[j] = temp_displacement_cm[idx2D(j,i,icount_cm)] * 3;
    }
    MPI_Type_create_indexed_block( icount_cm, 3, &displacement[0], MPI_FLOAT, &dtype_temp_t[0] );
    
    for( j = 0 ; j < icount_ic ; j++ ) {
      displacement[j] = temp_displacement_ic[idx2D(j,i,icount_ic)] * 3;
    }
    MPI_Type_create_indexed_block(icount_ic, 3, &displacement[0], MPI_FLOAT, &dtype_temp_t[1] );

    MPI_Type_create_struct( 2, &blocklength[0], &struct_displacement[0], &dtype_temp_t[0], &dtype_indexed_t );
    MPI_Type_commit( &dtype_indexed_t );
    
    MPI_Type_free( &dtype_temp_t[0] );
    MPI_Type_free( &dtype_temp_t[1] );

    free( displacement );
 
    if ( myrank == 0 ) {
      timing_record(1);
    }

//! =============== ping pong communication =================

    for( j=0 ; j<inner_loop ; j++ ) {
      if ( myrank == 0 ) {
//! pack the buffer
        pos = 0;
        MPI_Pack( MPI_BOTTOM, 1, dtype_indexed_t, &buffer[0], bytes, &pos, local_communicator );
        timing_record(2);
//! send the data from rank 0 to rank 1          
        MPI_Send( &buffer[0], pos, MPI_PACKED, 1, itag, local_communicator );
//! receive the data from rank 1 back
        MPI_Recv( &buffer[0], bytes, MPI_PACKED, 1, itag, local_communicator, MPI_STATUS_IGNORE );
        timing_record(3);
//! unpack the data
        pos = 0;
        MPI_Unpack( &buffer[0], bytes, &pos, MPI_BOTTOM, 1, dtype_indexed_t, local_communicator );
        timing_record(4);
//! now for rank 1
      } else {
//! receive the data from rank 0
        MPI_Recv( &buffer[0], bytes, MPI_PACKED, 0, itag, local_communicator, MPI_STATUS_IGNORE );
//! unpack the data
        pos = 0;
        MPI_Unpack( &buffer[0], bytes, &pos, MPI_BOTTOM, 1, dtype_indexed_t, local_communicator );
//! pack the data
        pos = 0;
        MPI_Pack( MPI_BOTTOM, 1, dtype_indexed_t, &buffer[0], bytes, &pos, local_communicator );
//! send the data back to rank 0
        MPI_Send( &buffer[0], pos, MPI_PACKED, 0, itag, local_communicator );
      }
    } //! inner loop

//! ======================= clean up ========================
    MPI_Type_free( &dtype_indexed_t );

    free( buffer );

    if ( myrank == 0 ) {
      timing_record(5);
    }

  } //! outer loop

  if ( myrank == 0 ) {
    timing_print( 1 );
  }

  free( array_cm );
  free( array_ic );

  free( temp_displacement_cm );
  free( temp_displacement_ic );

}

void timing_specfem3d_mt_ddt( int DIM1, int DIM2, int DIM3, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator ) {

  float* send_array;
  float* recv_array;

  int myrank;
  int i, j;
  int base, bytes, typesize;

  MPI_Datatype dtype_temp_t, dtype_send_t, dtype_recv_t;

  char method[50];

//! just some statements to prevent compiler warnings of unused variables
//! those parameter are included for future features
  *correct_flag = 0;
  *ptypesize = 0;
//  typesize = filehandle_debug

  send_array = malloc( DIM1 * DIM2 * DIM3 * sizeof(float) );
  recv_array = malloc( DIM1 * DIM3 * sizeof(float) );

  MPI_Comm_rank( local_communicator, &myrank );

  base = myrank * DIM1 * DIM2 * DIM3 + 1;
  utilities_fill_unique_array_3D_float( &send_array[0], DIM1, DIM2, DIM3, base );

  if ( myrank == 0 ) {
    snprintf(method, 50, "mpi_ddt" );

    MPI_Type_size( MPI_FLOAT, &typesize );
    bytes = DIM1 * DIM3 * typesize;

    timing_init( testname, &method[0], bytes );
  }

  for( i=0 ; i<outer_loop ; i++ ) {
        
    MPI_Type_contiguous( DIM1, MPI_FLOAT, &dtype_temp_t );
    MPI_Type_vector( DIM3, 1, DIM2, dtype_temp_t, &dtype_send_t );
    MPI_Type_commit( &dtype_send_t );
    MPI_Type_free( &dtype_temp_t );

    MPI_Type_contiguous( DIM3*DIM1, MPI_FLOAT, &dtype_recv_t );
    MPI_Type_commit( &dtype_recv_t );

    if ( myrank == 0 ) {
      timing_record(1);
    }

    for( j=0 ; j<inner_loop ; j++ ) {
          
      if ( myrank == 0 ) {
        MPI_Send( &send_array[0], 1, dtype_send_t, 1, itag, local_communicator );
        MPI_Recv( &recv_array[0], 1, dtype_recv_t, 1, itag, local_communicator, MPI_STATUS_IGNORE );
        timing_record(3);
      } else {
        MPI_Recv( &recv_array[0], 1, dtype_recv_t, 0, itag, local_communicator, MPI_STATUS_IGNORE );
        MPI_Send( &send_array[0], 1, dtype_send_t, 0, itag, local_communicator );
      }
    
    } //! inner loop

    MPI_Type_free( &dtype_send_t );
    MPI_Type_free( &dtype_recv_t );

    if ( myrank == 0 ) {
      timing_record(5);
    }

  } //! outer loop

  if ( myrank == 0 ) {
    timing_print( 1 );
  }

  free(send_array);
  free(recv_array);
}

void timing_specfem3d_mt_manual( int DIM1, int DIM2, int DIM3, int outer_loop, int inner_loop, int *correct_flag, int *ptypesize, char *testname, MPI_File filehandle_debug, MPI_Comm local_communicator ) {

  float* send_array;
  float* recv_array;
  
  float* buffer;

  int myrank;
  int i, j, k, base, bytes, typesize;

  char method[50];

//! just some statements to prevent compiler warnings of unused variables
//! those parameter are included for future features
  *correct_flag = 0;
  *ptypesize = 0;
// typesize = filehandle_debug

  send_array = malloc( DIM1 * DIM2 * DIM3 * sizeof(float) );
  recv_array = malloc( DIM1 * DIM3 * sizeof(float) );

  MPI_Comm_rank( local_communicator, &myrank );

  base = myrank * DIM1 * DIM2 * DIM3 + 1;
  utilities_fill_unique_array_3D_float( &send_array[0], DIM1, DIM2, DIM3, base );

  if ( myrank == 0 ) {
    snprintf( method, 50, "manual" );

    MPI_Type_size( MPI_FLOAT, &typesize );
    bytes = DIM1 * DIM3 * typesize ;

    timing_init( testname, &method[0], bytes );
  }

  for( i=0 ; i<outer_loop ; i++ ) {

    buffer = malloc( DIM1 * DIM3 * sizeof(float) );
    MPI_Type_size( MPI_FLOAT, &typesize );
    bytes = DIM1 * DIM3 * typesize;
        
    if ( myrank == 0 ) {
      timing_record(1);
    }

    for( j=0 ; j<inner_loop ; j++ ) {
          
      if ( myrank == 0 ) {
        for( k = 0; k<DIM3 ; k++ ) {
          memcpy( &buffer[k*DIM1], &send_array[k*DIM1*DIM2], DIM1*sizeof(float) );
        }
        timing_record(2);
        MPI_Send( &buffer[0], DIM1*DIM3, MPI_FLOAT, 1, itag, local_communicator );
        MPI_Recv( &recv_array[0], DIM1*DIM3, MPI_FLOAT, 1, itag, local_communicator, MPI_STATUS_IGNORE );
        timing_record(3);
      } else {
        MPI_Recv( &recv_array[0], DIM1*DIM3, MPI_FLOAT, 0, itag, local_communicator, MPI_STATUS_IGNORE );
        for( k = 0; k<DIM3 ; k++ ) {
          memcpy( &buffer[k*DIM1], &send_array[k*DIM1*DIM2], DIM1*sizeof(float) );
        }
        MPI_Send( &buffer[0], DIM1*DIM3, MPI_FLOAT, 0, itag, local_communicator );
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

  free( send_array );
  free( recv_array );
}

void timing_specfem3d_mt_mpi_pack_ddt( int DIM1, int DIM2, int DIM3, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator ) {

  float* send_array;
  float* recv_array;

  float* buffer;

  int myrank;
  int i, j, base, bytes, typesize, pos;

  MPI_Datatype dtype_temp_t, dtype_send_t, dtype_recv_t;

  char method[50];

//! just some statements to prevent compiler warnings of unused variables
//! those parameter are included for future features
  *correct_flag = 0;
  *ptypesize = 0;
//  typesize = filehandle_debug

  send_array = malloc( DIM1 * DIM2 * DIM3 * sizeof(float) );
  recv_array = malloc( DIM1 * DIM3 * sizeof(float) );

  MPI_Comm_rank( local_communicator, &myrank );

  base = myrank * DIM1 * DIM2 * DIM3 + 1;
  utilities_fill_unique_array_3D_float( &send_array[0], DIM1, DIM2, DIM3, base );

  if ( myrank == 0 ) {
    snprintf( method, 50, "mpi_pack_ddt" );

    MPI_Type_size( MPI_FLOAT, &typesize );
    bytes = DIM1 * DIM3 * typesize ;

    timing_init( testname, &method[0], bytes );
  }

  for( i=0 ; i<outer_loop ; i++ ) {

    buffer = malloc( DIM1 * DIM3 * sizeof(float) );
    MPI_Type_size( MPI_FLOAT, &typesize );
    bytes = DIM1 * DIM3 * typesize;
        
    MPI_Type_contiguous( DIM1, MPI_FLOAT, &dtype_temp_t );
    MPI_Type_vector( DIM3, 1, DIM2, dtype_temp_t, &dtype_send_t );
    MPI_Type_commit( &dtype_send_t );
    MPI_Type_free( &dtype_temp_t );

    MPI_Type_contiguous( DIM3*DIM1, MPI_FLOAT, &dtype_recv_t );
    MPI_Type_commit( &dtype_recv_t );

    if ( myrank == 0 ) {
      timing_record(1);
    }

    for( j=0 ; j<inner_loop ; j++ ) {
          
      if ( myrank == 0 ) {
        pos = 0;
        MPI_Pack( &send_array[0], 1, dtype_send_t, &buffer[0], bytes, &pos, local_communicator );
        timing_record(2);
        MPI_Send( &buffer[0], pos, MPI_PACKED, 1, itag, local_communicator );
        MPI_Recv( &buffer[0], bytes, MPI_PACKED, 1, itag, local_communicator, MPI_STATUS_IGNORE );
        timing_record(3);
        pos = 0;
        MPI_Unpack( &buffer[0], bytes, &pos, &recv_array[0], 1, dtype_recv_t, local_communicator );
        timing_record(4);
      } else {
        MPI_Recv( &buffer[0], bytes, MPI_PACKED, 0, itag, local_communicator, MPI_STATUS_IGNORE );
        pos = 0;
        MPI_Unpack( &buffer[0], bytes, &pos, &recv_array[0], 1, dtype_recv_t, local_communicator );
        pos = 0;
        MPI_Pack( &send_array[0], 1, dtype_send_t, &buffer[0], bytes, &pos, local_communicator );
        MPI_Send( &buffer[0], pos, MPI_PACKED, 0, itag, local_communicator );
      }
    
    } //! inner loop

    free( buffer );

    MPI_Type_free( &dtype_send_t );
    MPI_Type_free( &dtype_recv_t );

    if ( myrank == 0 ) {
      timing_record(5);
    }

  } //! outer loop

  if ( myrank == 0 ) {
    timing_print( 1 );
  }

  free( send_array );
  free( recv_array );
}
