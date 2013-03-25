//! contains the timing benchmarks for NAS tests (MG/LU)

#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>

#include "ddtbench.h"

#define itag 0

static inline int idx3D(int x, int y, int z, int DIM1, int DIM2) {
    return x+DIM1*(y+z*DIM2);
}

void timing_nas_lu_y_ddt( int DIM2, int DIM3, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator ) {

  int DIM1 = 5;

  double *array;

  int myrank;
  int i, j, base, bytes, typesize;

  MPI_Datatype dtype_y_t, dtype_temp_t;

  char method[50];

//! just some statements to prevent compiler warnings of unused variables      
//! those parameter are included for future features
  *correct_flag = 0;
  *ptypesize = 0;
  //typesize = filehandle_debug;
  
  array = malloc( DIM1*(DIM2+2)*(DIM3+2) * sizeof(double) );
      
  MPI_Comm_rank( local_communicator, &myrank );

  base = myrank * DIM1 * (DIM2+2) * (DIM3+2) + 1;
  utilities_fill_unique_array_3D_double( &array[0], DIM1, DIM2+2, DIM3+2, base );

  if ( myrank == 0 ) {
    snprintf(method,50,"mpi_ddt");

    MPI_Type_size( MPI_DOUBLE, &typesize );
    bytes = 5 * DIM3 * typesize;

    timing_init( testname, &method[0], bytes );
  }

  for( i=0 ; i < outer_loop ; i++ ) {

    MPI_Type_contiguous(5, MPI_DOUBLE, &dtype_temp_t );

    MPI_Type_vector( DIM3, 1, DIM2+2, dtype_temp_t, &dtype_y_t );
    MPI_Type_commit( &dtype_y_t );

    MPI_Type_free( &dtype_temp_t );

    if ( myrank == 0 ) {
      timing_record(1);
    }

    for (j=0 ; j<inner_loop ; j++ ) {

      if ( myrank == 0 ) {
        MPI_Send( &array[idx3D(0,DIM2,1,DIM1,DIM2+2)], 1, dtype_y_t, 1, itag, local_communicator );
        MPI_Recv( &array[idx3D(0,0,1,DIM1,DIM2+2)], 1, dtype_y_t, 1, itag, local_communicator, MPI_STATUS_IGNORE );
        timing_record(3);
      } else {
        MPI_Recv( &array[idx3D(0,0,1,DIM1,DIM2+2)], 1, dtype_y_t, 0, itag, local_communicator, MPI_STATUS_IGNORE );
        MPI_Send( &array[idx3D(0,DIM2,1,DIM1,DIM2+2)], 1, dtype_y_t, 0, itag, local_communicator );
      }

    }

    MPI_Type_free( &dtype_y_t );

    if ( myrank == 0 ) {
      timing_record(5);
    }

  } 

  if ( myrank == 0 ) {
    timing_print( 1 );
  }

  free(array);
}

void timing_nas_lu_y_manual( int DIM2, int DIM3, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator ) {

  int DIM1 = 5;

  double* array;
  double* buffer;

  int myrank;
  int i, j, k, l,  base, bytes, typesize;

  char method[50];

//! just some statements to prevent compiler warnings of unused variables      
//! those parameter are included for future features
  *correct_flag = 0;
  *ptypesize = 0;
//  typesize = filehandle_debug

  array = malloc( DIM1*(DIM2+2)*(DIM3+2) * sizeof(double) );

  MPI_Comm_rank( local_communicator, &myrank );

  base = myrank * DIM1 * (DIM2+2) * (DIM3+2) + 1;
  utilities_fill_unique_array_3D_double( &array[0], DIM1, DIM2+2, DIM3+2, base );
      
  if ( myrank == 0 ) {
    snprintf(method, 50, "manual");

    MPI_Type_size( MPI_DOUBLE, &typesize );
    bytes = 5 * DIM3 * typesize;

    timing_init( testname, &method[0], bytes );
  }

  for( i=0 ; i<outer_loop ; i++ ) {

    buffer = malloc( DIM1 * DIM3 * sizeof(double));
        
    if ( myrank == 0 ) {
       timing_record(1);
    }

    for( j=0 ; j<inner_loop ; j++ ) {
      if ( myrank == 0 ) {
//! pack the data
        base = 0;
        for( k=1 ; k<DIM3 ; k++ ) {
          for( l=0 ; l<DIM1 ; l++ ) {
            buffer[base++] = array[idx3D(l,DIM2,k,DIM1,DIM2+2)];
          }
        }
        timing_record(2);
        MPI_Send( &buffer[0], DIM1*DIM3, MPI_DOUBLE, 1, itag, local_communicator );
        MPI_Recv( &buffer[0], DIM1*DIM3, MPI_DOUBLE, 1, itag, local_communicator, MPI_STATUS_IGNORE );
        timing_record(3);
//! unpack the data
        base = 0;
        for( k=1 ; k<DIM3 ; k++ ) {
          for( l=0 ; l<DIM1 ; l++ ) {
            array[idx3D(l,0,k,DIM1,DIM2+2)] = buffer[base++];
          }
        }
        timing_record(4);
      } else {
        MPI_Recv( &buffer[0], DIM1*DIM3, MPI_DOUBLE, 0, itag, local_communicator, MPI_STATUS_IGNORE );
//! unpack the data
        base = 0;
        for( k=1 ; k<DIM3 ; k++ ) {
          for( l=0 ; l<DIM1 ; l++ ) {
            array[idx3D(l,0,k,DIM1,DIM2+2)] = buffer[base++];
          }
        }
//! pack the data
        base = 0;
        for( k=1 ; k<DIM3 ; k++ ) {
          for( l=0 ; l<DIM1 ; l++ ) {
            buffer[base++] = array[idx3D(l,DIM2,k,DIM1,DIM2+2)];
          }
        }
        MPI_Send( &buffer[0], DIM1*DIM3, MPI_DOUBLE, 0, itag, local_communicator );
      }
    } //! inner loop

    free( buffer );

    if ( myrank == 0 ) {
      timing_record(5);
    }
  }

  if ( myrank == 0 ) {
    timing_print( 1 );
  }
    
  free(array);
}

void timing_nas_lu_y_mpi_pack_ddt( int DIM2, int DIM3, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator ) {

  int DIM1  = 5;

  double* array;
  double* buffer;

  int myrank;
  int i, j, base, bytes, typesize, pos;

  MPI_Datatype dtype_y_t, dtype_temp_t;

  char method[50];

//! just some statements to prevent compiler warnings of unused variables      
//! those parameter are included for future features
  *correct_flag = 0;
  *ptypesize = 0;
//  typesize = filehandle_debug
  
  array = malloc( DIM1 * (DIM2+2) * (DIM3+2) * sizeof(double));
     
  MPI_Comm_rank( local_communicator, &myrank );

  base = myrank * DIM1 * (DIM2+2) * (DIM3+2) + 1;
  utilities_fill_unique_array_3D_double( &array[0], DIM1, DIM2+2, DIM3+2, base);

  if ( myrank == 0 ) {
    snprintf(method, 50, "mpi_pack_ddt");

    MPI_Type_size( MPI_DOUBLE, &typesize );
    bytes = DIM1 * DIM3 * typesize;
         
    timing_init( testname, &method[0], bytes );
  }

  for( i=0 ; i<outer_loop ; i++ ) {
    buffer = malloc( DIM1 * DIM3 * sizeof(double));
    MPI_Type_size( MPI_DOUBLE, &typesize );
    bytes = DIM1 * DIM3 * typesize;

    MPI_Type_contiguous(5, MPI_DOUBLE, &dtype_temp_t );

    MPI_Type_vector( DIM3, 1, DIM2+2, dtype_temp_t, &dtype_y_t );
    MPI_Type_commit( &dtype_y_t );

    MPI_Type_free( &dtype_temp_t );
       
    if ( myrank == 0 ) {
      timing_record(1);
    }

    for( j=0 ; j<inner_loop ; j++ ) {
      if ( myrank == 0 ) {
//! pack the data
        pos = 0;
        MPI_Pack( &array[idx3D(0,DIM2,1,DIM1,DIM2)], 1, dtype_y_t, &buffer[0], bytes, &pos, local_communicator );
        timing_record(2);
        MPI_Send( &buffer[0], pos, MPI_PACKED, 1, itag, local_communicator );
        MPI_Recv( &buffer[0], bytes, MPI_PACKED, 1, itag, local_communicator, MPI_STATUS_IGNORE );
        timing_record(3);
//! unpack the data
        pos = 0;
        MPI_Unpack( &buffer[0], bytes, &pos, &array[idx3D(0,0,1,DIM1,DIM2+2)], 1, dtype_y_t, local_communicator );
        timing_record(4);
      } else {
        MPI_Recv( &buffer[0], bytes, MPI_PACKED, 0, itag, local_communicator, MPI_STATUS_IGNORE );
//! unpack the data
        pos = 0;
        MPI_Unpack( &buffer[0], bytes, &pos, &array[idx3D(0,0,1,DIM1,DIM2+2)], 1, dtype_y_t, local_communicator );
//! pack the data
        pos = 0;
        MPI_Pack( &array[idx3D(0,DIM2,1,DIM1,DIM2+2)], 1, dtype_y_t, &buffer[0], bytes, &pos, local_communicator );
        MPI_Send( &buffer[0], pos, MPI_PACKED, 0, itag, local_communicator );
      }
    } //! inner loop

    free( buffer );

    MPI_Type_free( &dtype_y_t );

    if ( myrank == 0 ) {
      timing_record(5);
    }

  } //! outer loop 

  if ( myrank == 0 ) {
    timing_print( 1 );
  }

  free(array);

}

void timing_nas_lu_x_ddt( int DIM2, int DIM3, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator ) {

  int DIM1 = 5;
    
  double* array;

  int myrank;
  int i, j, base, bytes, typesize;

  MPI_Datatype dtype_x_t;

  char method[50];

//! just some statements to prevent compiler warnings of unused variables      
//! those parameter are included for future features
  *correct_flag = 0;
  *ptypesize = 0;
//  typesize = filehandle_debug

  array = malloc( DIM1 * (DIM2+2) * (DIM3+2) * sizeof(double));

  MPI_Comm_rank( local_communicator, &myrank );

  base = myrank * DIM1 * (DIM2+2) * (DIM3+2) + 1;
  utilities_fill_unique_array_3D_double( &array[0], DIM1, DIM2+2, DIM3+2, base );

  if ( myrank == 0 ) {
    snprintf(method, 50, "mpi_ddt");

    MPI_Type_size( MPI_DOUBLE, &typesize );
    bytes = DIM1 * DIM2 * typesize;

    timing_init( testname, &method[0], bytes );
  }

  for( i=0 ; i<outer_loop ; i++ ) {

    MPI_Type_contiguous(DIM1*DIM2, MPI_DOUBLE, &dtype_x_t );

    MPI_Type_commit( &dtype_x_t );

    if ( myrank == 0 ) {
      timing_record(1);
    }

    for( j=0 ; j<inner_loop ; j++ ) {

      if ( myrank == 0 ) {
        MPI_Send( &array[idx3D(0,1,DIM3,DIM1,DIM2+2)], 1, dtype_x_t, 1, itag, local_communicator );
        MPI_Recv( &array[idx3D(0,1,0,DIM1,DIM2+2)], 1, dtype_x_t, 1, itag, local_communicator, MPI_STATUS_IGNORE );
        timing_record(3);
      } else {
        MPI_Recv( &array[idx3D(0,1,0,DIM1,DIM2+2)], 1, dtype_x_t, 0, itag, local_communicator, MPI_STATUS_IGNORE );
        MPI_Send( &array[idx3D(0,1,DIM3,DIM1,DIM2+2)], 1, dtype_x_t, 0, itag, local_communicator );
      }
    }

    MPI_Type_free( &dtype_x_t );

    if ( myrank == 0 ) {
      timing_record(5);
    }
  }

  if ( myrank == 0 ) {
    timing_print( 1 );
  }

  free(array);
}

void timing_nas_lu_x_manual( int DIM2, int DIM3, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator ) {

  int DIM1 = 5;

  double* array;
  double* buffer;

  int myrank;
  int i, j, k, l,  base, bytes, typesize;

  char method[50];

//! just some statements to prevent compiler warnings of unused variables      
//! those parameter are included for future features
  *correct_flag = 0;
  *ptypesize = 0;
//  typesize = filehandle_debug;

  array = malloc( DIM1 * (DIM2+2) * (DIM3+2) * sizeof(double));

  MPI_Comm_rank( local_communicator, &myrank );

  base = myrank * DIM1 * (DIM2+2) * (DIM3+2) + 1;

  utilities_fill_unique_array_3D_double( &array[0], DIM1, DIM2+2, DIM3+2, base );

  if ( myrank == 0 ) {
    snprintf(method, 50, "manual");

    MPI_Type_size( MPI_DOUBLE, &typesize );
    bytes = DIM1 * DIM2 * typesize;

    timing_init( testname, &method[0], bytes );
  }

  for( i=0 ; i<outer_loop ; i++ ) {
    
    buffer = malloc(DIM1 * DIM2 * sizeof(double));
        
    if ( myrank == 0 ) {
      timing_record(1);
    }

    for( j=0 ; j<inner_loop ; j++ ) {
      if ( myrank == 0 ) {
//! pack the data
        base = 0;
        for( k=1 ; k<DIM2+1 ; k++ ) {
          for( l=0 ; l<DIM1 ; l++ ) {
            buffer[base++] = array[idx3D(l,k,DIM3,DIM1,DIM2+2)];
          }
        }
        timing_record(2);
        MPI_Send( &buffer[0], DIM1*DIM2, MPI_DOUBLE, 1, itag, local_communicator );
        MPI_Recv( &buffer[0], DIM1*DIM2, MPI_DOUBLE, 1, itag, local_communicator, MPI_STATUS_IGNORE );
        timing_record(3);
//! unpack the data
        base = 0;
        for( k=1 ; k<DIM2+1 ; k++ ) {
          for( l=0 ; l<DIM1 ; l++ ) {
            array[idx3D(l,k,0,DIM1,DIM2+2)] = buffer[base++];
          }
        }
        timing_record(4);
      } else {
        MPI_Recv( &buffer[0], DIM1*DIM2, MPI_DOUBLE, 0, itag, local_communicator, MPI_STATUS_IGNORE );
//! unpack the data
        base = 0;
        for( k=1 ; k<DIM2+1 ; k++ ) {
          for( l=0 ; l<DIM1 ; l++ ) {
            array[idx3D(l,k,0,DIM1,DIM2+2)] = buffer[base++];
          }              
        }
//! pack the data
        base = 0;
        for( k=1 ; k<DIM2+1 ; k++ ) {
          for( l=0 ; l<DIM1 ; l++ ) {
            buffer[base++] = array[idx3D(l,k,DIM3,DIM1,DIM2+2)];
          }
        }
        MPI_Send( &buffer[0], DIM1*DIM2, MPI_DOUBLE, 0, itag, local_communicator );
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

  free(array);
}

void timing_nas_lu_x_mpi_pack_ddt( int DIM2, int DIM3, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator ) {

  int DIM1 = 5;

  double* array;
  double* buffer;

  int myrank;
  int i, j, base, bytes, typesize, pos;

  MPI_Datatype dtype_x_t;

  char method[50];

//! just some statements to prevent compiler warnings of unused variables      
//! those parameter are included for future features
  *correct_flag = 0;
  *ptypesize = 0;
//  typesize = filehandle_debug

  array = malloc( DIM1 * (DIM2+2) * (DIM3+2) * sizeof(double));
       
  MPI_Comm_rank( local_communicator, &myrank );

  base = myrank * DIM1 * (DIM2+2) * (DIM3+2) + 1;
  utilities_fill_unique_array_3D_double( &array[0], DIM1, DIM2+2, DIM3+2, base );

  if ( myrank == 0 ) {
    snprintf(method, 50, "mpi_pack_ddt");
        
    MPI_Type_size( MPI_DOUBLE, &typesize );
    bytes = DIM1 * DIM2 * typesize;

    timing_init( testname, &method[0], bytes );
  }
  
  for( i=0 ; i<outer_loop ; i++ ) {

    buffer = malloc( DIM1 * DIM2 * sizeof(double));
    MPI_Type_size( MPI_DOUBLE, &typesize );
    bytes = DIM1 * DIM2 * typesize;

    MPI_Type_contiguous( DIM1*DIM2, MPI_DOUBLE, &dtype_x_t );
    MPI_Type_commit( &dtype_x_t );

    if ( myrank == 0 ) {
      timing_record(1);
    }

    for( j=0 ; j<inner_loop ; j++ ) {
      if ( myrank == 0 ) {
//! pack the data
        pos = 0;
        MPI_Pack( &array[idx3D(0,1,DIM3,DIM1,DIM2+2)], 1, dtype_x_t, &buffer[0], bytes, &pos, local_communicator );
        timing_record(2);
        MPI_Send( &buffer[0], pos, MPI_PACKED, 1, itag, local_communicator );
        MPI_Recv( &buffer[0], bytes, MPI_PACKED, 1, itag, local_communicator, MPI_STATUS_IGNORE );
        timing_record(3);
//! unpack the data
        pos = 0;
        MPI_Unpack( &buffer[0], bytes, &pos, &array[idx3D(0,1,0,DIM1,DIM2+2)], 1, dtype_x_t, local_communicator );
        timing_record(4);
      } else {
        MPI_Recv( &buffer[0], bytes, MPI_PACKED, 0, itag, local_communicator, MPI_STATUS_IGNORE );
//! unpack the data
        pos = 0;
        MPI_Unpack( &buffer[0], bytes, &pos, &array[idx3D(0,1,0,DIM1,DIM2+2)], 1, dtype_x_t, local_communicator );
//! pack the data
        pos = 0;
        MPI_Pack( &array[idx3D(0,1,DIM3,DIM1,DIM2+2)], 1, dtype_x_t, &buffer[0], bytes, &pos, local_communicator );
        MPI_Send( &buffer[0], pos, MPI_PACKED, 0, itag, local_communicator );
      }
    } //! inner loop

    free( buffer );

    MPI_Type_free( &dtype_x_t );

    if ( myrank == 0 ) {
      timing_record(5);
    }
  } //! outer loop

  if ( myrank == 0 ) {
    timing_print( 1 );
  }

  free(array);
}

void timing_nas_mg_x_ddt( int DIM1, int DIM2, int DIM3, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator ) {

  double* array;
      
  int myrank;
  int base, i, j, typesize, bytes;

  MPI_Aint stride;
  MPI_Datatype dtype_face_x_t, dtype_temp_t;

  char method[50];

//! just some statements to prevent compiler warnings of unused variables      
//! those parameter are included for future features
  *correct_flag = 0;
  *ptypesize = 0;
//  typesize = filehandle_debug
  
  array = malloc( DIM1 * DIM2 * DIM3 * sizeof(double));

  MPI_Comm_rank( local_communicator, &myrank );

  base = myrank * DIM1 * DIM2 * DIM3 + 1;
  utilities_fill_unique_array_3D_double( &array[0], DIM1, DIM2, DIM3, base);

  if ( myrank == 0 ) {
    snprintf(method, 50, "mpi_ddt");
  
    MPI_Type_size( MPI_DOUBLE, &typesize );
    bytes = (DIM2-2)*(DIM3-2) * typesize;

    timing_init( testname, &method[0], bytes );
  }

  for( i=0 ; i<outer_loop ; i++ ) {
  
    MPI_Type_vector( DIM2-2, 1, DIM1, MPI_DOUBLE, &dtype_temp_t );

    MPI_Type_size( MPI_DOUBLE, &typesize );
    stride = DIM1 * DIM2 * typesize;

    MPI_Type_create_hvector( DIM3-2, 1, stride, dtype_temp_t, &dtype_face_x_t );
    MPI_Type_commit( &dtype_face_x_t );

    MPI_Type_free( &dtype_temp_t );

    if ( myrank == 0 ) {
      timing_record(1);
    }

    for( j=0 ; j<inner_loop ; j++ ) {

      if ( myrank == 0 ) {
        MPI_Send( &array[idx3D(DIM1-2,1,1,DIM1,DIM2)], 1, dtype_face_x_t, 1, itag, local_communicator );
        MPI_Recv( &array[idx3D(DIM1-1,1,1,DIM1,DIM2)], 1, dtype_face_x_t, 1, itag, local_communicator, MPI_STATUS_IGNORE );
        timing_record(3);
      } else {
        MPI_Recv( &array[idx3D(DIM1-1,1,1,DIM1,DIM2)], 1, dtype_face_x_t, 0, itag, local_communicator, MPI_STATUS_IGNORE );
        MPI_Send( &array[idx3D(DIM1-2,1,1,DIM1,DIM2)], 1, dtype_face_x_t, 0, itag, local_communicator );
      }

    } //! inner loop

    MPI_Type_free( &dtype_face_x_t );

    if ( myrank == 0 ) {
      timing_record(5);
    }

  }

  if ( myrank == 0 ) {
    timing_print( 1 );
  }

  free(array);
}

void timing_nas_mg_x_manual( int DIM1, int DIM2, int DIM3, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator ) {
      
  double* array;
  double* buffer;

  int myrank;
  int base, i, j, k, l, typesize, bytes, psize;

  char method[50];

//! just some statements to prevent compiler warnings of unused variables      
//! those parameter are included for future features
  *correct_flag = 0;
  *ptypesize = 0;
//  typesize = filehandle_debug

  array = malloc( DIM1 * DIM2 * DIM3 * sizeof(double));

  MPI_Comm_rank( local_communicator, &myrank );

  base = myrank * DIM1 * DIM2 * DIM3 + 1;
  utilities_fill_unique_array_3D_double( &array[0], DIM1, DIM2, DIM3, base);

  if ( myrank == 0 ) {
    snprintf(method, 50, "manual");

    MPI_Type_size( MPI_DOUBLE, &typesize );
    bytes = (DIM2-2)*(DIM3-2) * typesize;

    timing_init( testname, &method[0], bytes );
  }

  for( i=0 ; i<outer_loop ; i++ ) {

    psize = (DIM2-2)*(DIM3-2);
    buffer = malloc( psize * sizeof(double) );

    if ( myrank == 0 ) {
      timing_record(1);
    }

    for( j=0 ; j<inner_loop ; j++ ) {

      if ( myrank == 0 ) {
        base = 0;
        for( k=1 ; k<DIM3-1 ; k++ ) {
          for( l=1 ; l<DIM2-1 ; l++ ) {
            buffer[base++] = array[idx3D(DIM1-2,l,k,DIM1,DIM2)];
          }
        }
        timing_record(2);
        MPI_Send( &buffer[0], psize, MPI_DOUBLE, 1, itag, local_communicator );
        MPI_Recv( &buffer[0], psize, MPI_DOUBLE, 1, itag, local_communicator, MPI_STATUS_IGNORE );
        timing_record(3);
        base = 0;
        for( k=1 ; k<DIM3-1 ; k++ ) {
          for( l=1 ; l<DIM2-1 ; l++ ) {
            array[idx3D(DIM1-1,l,k,DIM1,DIM2)] = buffer[base++];
          }
        }
        timing_record(4);
      } else {
        MPI_Recv( &buffer[0], psize, MPI_DOUBLE, 0, itag, local_communicator, MPI_STATUS_IGNORE );
        base = 0;
        for( k=1 ; k<DIM3-1 ; k++ ) {
          for( l=1 ; l<DIM2-1 ; l++ ) {
            array[idx3D(DIM1-1,l,k,DIM1,DIM2)] = buffer[base++];
          }
        }
        base = 0;
        for( k=1 ; k<DIM3-1 ; k++ ) {
          for( l=1 ; l<DIM2-1 ; l++ ) {
            buffer[base++] = array[idx3D(DIM1-2,l,k,DIM1,DIM2)];
          }
        }
        MPI_Send( &buffer[0], psize, MPI_DOUBLE, 0, itag, local_communicator );
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

  free(array);
}
      
void timing_nas_mg_x_mpi_pack_ddt( int DIM1, int DIM2, int DIM3, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator ) {
      
  double* array;
  double* buffer;
      
  int myrank;
  int base, i, j, typesize, bytes, pos;

  MPI_Aint stride;
  MPI_Datatype dtype_face_x_t, dtype_temp_t;

  char method[50];

//! just some statements to prevent compiler warnings of unused variables      
//! those parameter are included for future features
  *correct_flag = 0;
  *ptypesize = 0;
//  typesize = filehandle_debug

  array = malloc( DIM1 * DIM2 * DIM3 * sizeof(double) );

  MPI_Comm_rank( local_communicator, &myrank );

  base = myrank * DIM1 * DIM2 * DIM3 + 1;
  utilities_fill_unique_array_3D_double( &array[0], DIM1, DIM2, DIM3, base);

  if ( myrank == 0 ) {
    snprintf(method, 50, "mpi_pack_ddt");

    MPI_Type_size( MPI_DOUBLE, &typesize );
    bytes = (DIM2-2)*(DIM3-2) * typesize;

    timing_init( testname, &method[0], bytes );
  }

  for( i=0 ; i<outer_loop ; i++ ) {

    buffer = malloc( (DIM2-2) * (DIM3-2) * sizeof(double) );
    MPI_Type_size( MPI_DOUBLE, &typesize );
    bytes = (DIM2-2)*(DIM3-2) * typesize;
  
    MPI_Type_vector( DIM2-2, 1, DIM1, MPI_DOUBLE, &dtype_temp_t );

    MPI_Type_size( MPI_DOUBLE, &typesize );
    stride = DIM1 * DIM2 * typesize;

    MPI_Type_create_hvector( DIM3-2, 1, stride, dtype_temp_t, &dtype_face_x_t );
    MPI_Type_commit( &dtype_face_x_t );

    MPI_Type_free( &dtype_temp_t );

    if ( myrank == 0 ) {
      timing_record(1);
    }

    for( j=0 ; j<inner_loop ; j++ ) {

      if ( myrank == 0 ) {
        pos = 0;
        MPI_Pack( &array[idx3D(DIM1-2,1,1,DIM1,DIM2)], 1, dtype_face_x_t, &buffer[0], bytes, &pos, local_communicator );
        timing_record(2);
        MPI_Send( &buffer[0], pos, MPI_PACKED, 1, itag, local_communicator );
        MPI_Recv( &buffer[0], bytes, MPI_PACKED, 1, itag, local_communicator, MPI_STATUS_IGNORE );
        timing_record(3);
        pos = 0;
        MPI_Unpack( &buffer[0], bytes, &pos, &array[idx3D(DIM1-1,1,1,DIM1,DIM2)], 1, dtype_face_x_t, local_communicator );
        timing_record(4);
      } else {
        MPI_Recv( &buffer[0], bytes, MPI_PACKED, 0, itag, local_communicator, MPI_STATUS_IGNORE );
        pos = 0;
        MPI_Unpack( &buffer[0], bytes, &pos, &array[idx3D(DIM1-1,1,1,DIM1,DIM2)], 1, dtype_face_x_t, local_communicator );
        pos = 0;
        MPI_Pack( &array[idx3D(DIM1-2,1,1,DIM1,DIM2)], 1, dtype_face_x_t, &buffer[0], bytes, &pos, local_communicator );
        MPI_Send( &buffer[0], pos, MPI_PACKED, 0, itag, local_communicator );
      }

    } //! inner loop

    free( buffer );

    MPI_Type_free( &dtype_face_x_t );

    if ( myrank == 0 ) {
      timing_record(5);
    }

  } //! outer loop

  if ( myrank == 0 ) {
    timing_print( 1 );
  }

  free(array);
}

void timing_nas_mg_y_ddt( int DIM1, int DIM2, int DIM3, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator ) {
      
  double* array;

  int myrank;
  int base, i, j, typesize, bytes;

  MPI_Datatype dtype_face_y_t;

  char method[50];

//! just some statements to prevent compiler warnings of unused variables      
//! those parameter are included for future features
  *correct_flag = 0;
  *ptypesize = 0;
//  typesize = filehandle_debug

  array = malloc( DIM1 * DIM2 * DIM3 * sizeof(double) );
  
  MPI_Comm_rank( local_communicator, &myrank );

  base = myrank * DIM1 * DIM2 * DIM3 + 1;
  utilities_fill_unique_array_3D_double( &array[0], DIM1, DIM2, DIM3, base );

  if ( myrank == 0 ) {
    snprintf(method, 50, "mpi_ddt");

    MPI_Type_size( MPI_DOUBLE, &typesize );
    bytes = (DIM1-2) * (DIM3-2) * typesize;

    timing_init( testname, &method[0], bytes );
  }

  for( i=0 ; i<outer_loop ; i++ ) {
  
    MPI_Type_vector( DIM3-2, DIM1-2, DIM1*DIM2, MPI_DOUBLE, &dtype_face_y_t );

    MPI_Type_commit( &dtype_face_y_t );

    if ( myrank == 0 ) {
      timing_record(1);
    }

    for( j=0 ; j<inner_loop ; j++ ) {

      if ( myrank == 0 ) {
        MPI_Send( &array[idx3D(1,DIM2-2,1,DIM1,DIM2)], 1, dtype_face_y_t, 1, itag, local_communicator );
        MPI_Recv( &array[idx3D(1,DIM2-1,1,DIM1,DIM2)], 1, dtype_face_y_t, 1, itag, local_communicator, MPI_STATUS_IGNORE );
        timing_record(3);
      } else {
        MPI_Recv( &array[idx3D(1,DIM2-1,1,DIM1,DIM2)], 1, dtype_face_y_t, 0, itag, local_communicator, MPI_STATUS_IGNORE );
        MPI_Send( &array[idx3D(1,DIM2-2,1,DIM1,DIM2)], 1, dtype_face_y_t, 0, itag, local_communicator );
      }

    } //! inner loop

    MPI_Type_free( &dtype_face_y_t );

    if ( myrank == 0 ) {
      timing_record(5);
    }

  } // outer loop

  if ( myrank == 0 ) {
    timing_print( 1 );
  }

  free(array);
}

void timing_nas_mg_y_manual( int DIM1, int DIM2, int DIM3, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator ) {

 double* array;
 double* buffer;

 int myrank;
 int base, i, j, k, l, typesize, bytes, psize;

 char method[50];

//! just some statements to prevent compiler warnings of unused variables      
//! those parameter are included for future features
 *correct_flag = 0;
 *ptypesize = 0;
// typesize = filehandle_debug
  
  array = malloc( DIM1 * DIM2 * DIM3 * sizeof(double) );

  MPI_Comm_rank( local_communicator, &myrank );

  base = myrank * DIM1 * DIM2 * DIM3 + 1;
  utilities_fill_unique_array_3D_double( &array[0], DIM1, DIM2, DIM3, base);

  if ( myrank == 0 ) {
    snprintf(method, 50, "manual");

    MPI_Type_size( MPI_DOUBLE, &typesize );
    bytes = (DIM1-2)*(DIM3-2) * typesize;

    timing_init( testname, method, bytes );
  }

  for( i=0 ; i<outer_loop ; i++ ) {

    psize = (DIM1-2)*(DIM3-2);
    buffer = malloc( psize * sizeof(double) );

    if ( myrank == 0 ) {
      timing_record(1);
    }

    for( j=0 ; j<inner_loop ; j++ ) {

      if ( myrank == 0 ) {
        base = 0;
        for( k=1 ; k<DIM3-1 ; k++ ) {
          for( l=1 ; l<DIM1-1 ; l++ ) {
            buffer[base++] = array[idx3D(l,DIM2-2,k,DIM1,DIM2)];
          }
        }
        timing_record(2);
        MPI_Send( &buffer[0], psize, MPI_DOUBLE, 1, itag, local_communicator );
        MPI_Recv( &buffer[0], psize, MPI_DOUBLE, 1, itag, local_communicator, MPI_STATUS_IGNORE );
        timing_record(3);
        base = 0;
        for( k=1 ; k<DIM3-1 ; k++ ) {
          for( l=1 ; l<DIM1-1 ; l++ ) {
            array[idx3D(l,DIM2-1,k,DIM1,DIM2)] = buffer[base++];
          }
        }
        timing_record(4);
      } else {
        MPI_Recv( &buffer[0], psize, MPI_DOUBLE, 0, itag, local_communicator, MPI_STATUS_IGNORE );
        base = 0;
        for( k=1 ; k<DIM3-1 ; k++ ) {
          for( l=1 ; l<DIM1-1 ; l++ ) {
            array[idx3D(l,DIM2-1,k,DIM1,DIM2)] = buffer[base++];
          }
        }
        base = 0;
        for( k=1 ; k<DIM3-1 ; k++ ) {
          for( l=1 ; l<DIM1-1 ; l++ ) {
            buffer[base++] = array[idx3D(l,DIM2-2,k,DIM1,DIM2)];
          }
        }
        MPI_Send( &buffer[0], psize, MPI_DOUBLE, 0, itag, local_communicator );
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

  free(array);
}
      
void timing_nas_mg_y_mpi_pack_ddt( int DIM1, int DIM2, int DIM3, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator ) {

  double* array;
  double* buffer;

  int myrank;
  int base, i, j, typesize, bytes, pos;

  MPI_Datatype dtype_face_y_t;

  char method[50];

//! just some statements to prevent compiler warnings of unused variables      
//! those parameter are included for future features
  *correct_flag = 0;
  *ptypesize = 0;
//  typesize = filehandle_debug

  array = malloc( DIM1 * DIM2 * DIM3 * sizeof(double) );

  MPI_Comm_rank( local_communicator, &myrank );

  base = myrank * DIM1 * DIM2 * DIM3 + 1;
  utilities_fill_unique_array_3D_double( &array[0], DIM1, DIM2, DIM3, base);

  if ( myrank == 0 ) {
    snprintf( method, 50, "mpi_pack_ddt" );

    MPI_Type_size( MPI_DOUBLE, &typesize );
    bytes = (DIM1-2) * (DIM3-2) * typesize;

    timing_init( testname, &method[0], bytes );
  }

  for( i=0 ; i<outer_loop ; i++ ) {

    buffer = malloc( (DIM1-2) * (DIM3-2) * sizeof(double) );
    MPI_Type_size( MPI_DOUBLE, &typesize );
    bytes = (DIM1-2) * (DIM3-2) * typesize;
  
    MPI_Type_vector( DIM3-2, DIM1-2, DIM1*DIM2, MPI_DOUBLE, &dtype_face_y_t );
    MPI_Type_commit( &dtype_face_y_t );

    if ( myrank == 0 ) {
      timing_record(1);
    }

    for( j=0 ; j<inner_loop ; j++ ) {

      if ( myrank == 0 ) {
        pos = 0;
        MPI_Pack( &array[idx3D(1,DIM2-2,1,DIM1,DIM2)], 1, dtype_face_y_t, &buffer[0], bytes, &pos, local_communicator );
        timing_record(2);
        MPI_Send( &buffer[0], pos, MPI_PACKED, 1, itag, local_communicator );
        MPI_Recv( &buffer[0], bytes, MPI_PACKED, 1, itag, local_communicator, MPI_STATUS_IGNORE );
        timing_record(3);
        pos = 0;
        MPI_Unpack( &buffer[0], bytes, &pos, &array[idx3D(1,DIM2-1,1,DIM1,DIM2)], 1, dtype_face_y_t, local_communicator );
        timing_record(4);
      } else {
        MPI_Recv( &buffer[0], bytes, MPI_PACKED, 0, itag, local_communicator, MPI_STATUS_IGNORE );
        pos = 0;
        MPI_Unpack( &buffer[0], bytes, &pos, &array[idx3D(1,DIM2-1,1,DIM1,DIM2)], 1, dtype_face_y_t, local_communicator );
        pos = 0;
        MPI_Pack( &array[idx3D(1,DIM2-2,1,DIM1,DIM2)], 1, dtype_face_y_t, &buffer[0], bytes, &pos, local_communicator );
        MPI_Send( &buffer[0], pos, MPI_PACKED, 0, itag, local_communicator );
      }

    } //! inner loop

    free( buffer );

    MPI_Type_free( &dtype_face_y_t );

    if ( myrank == 0 ) {
      timing_record(5);
    }

  } //! outer_loop

  if ( myrank == 0 ) {
    timing_print( 1 );
  }

  free( array );
}

void timing_nas_mg_z_ddt( int DIM1, int DIM2, int DIM3, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator ) {

  double* array;

  int myrank;
  int base, i, j, typesize, bytes;

  MPI_Datatype dtype_face_z_t;

  char method[50];

//! just some statements to prevent compiler warnings of unused variables      
//! those parameter are included for future features
  *correct_flag = 0;
  *ptypesize = 0;
//  typesize = filehandle_debug

  array = malloc( DIM1 * DIM2 * DIM3 * sizeof(double));
  
  MPI_Comm_rank( local_communicator, &myrank );

  base = myrank * DIM1 * DIM2 * DIM3 + 1;
  utilities_fill_unique_array_3D_double( &array[0], DIM1, DIM2, DIM3, base);

  if ( myrank == 0 ) {
    snprintf( method, 50, "mpi_ddt" );

    MPI_Type_size( MPI_DOUBLE, &typesize );
    bytes = (DIM1-2) * (DIM2-2) * typesize;

    timing_init( testname, &method[0], bytes );
  }

  for( i=0 ; i<outer_loop ; i++ ) {
  
    MPI_Type_vector( DIM2-2, DIM1-2, DIM1, MPI_DOUBLE, &dtype_face_z_t );
    MPI_Type_commit( &dtype_face_z_t );

    if ( myrank == 0 ) {
      timing_record(1);
    }

    for( j=0 ; j<inner_loop ; j++ ) {

      if ( myrank == 0 ) {
        MPI_Send( &array[idx3D(1,1,1,DIM1,DIM2)], 1, dtype_face_z_t, 1, itag, local_communicator );
        MPI_Recv( &array[idx3D(1,1,0,DIM1,DIM2)], 1, dtype_face_z_t, 1, itag, local_communicator, MPI_STATUS_IGNORE );
        timing_record(3);
      } else {
        MPI_Recv( &array[idx3D(1,1,0,DIM1,DIM2)], 1, dtype_face_z_t, 0, itag, local_communicator, MPI_STATUS_IGNORE );
        MPI_Send( &array[idx3D(1,1,1,DIM1,DIM2)], 1, dtype_face_z_t, 0, itag, local_communicator );
      }

    } //! inner loop

    MPI_Type_free( &dtype_face_z_t );

    if ( myrank == 0 ) {
      timing_record(5);
    }

  }

  if ( myrank == 0 ) {
    timing_print( 1 );
  }
  
  free(array);
}

void timing_nas_mg_z_manual( int DIM1, int DIM2, int DIM3, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator ) {

  double* array;
  double* buffer;

  int myrank;
  int base, i, j, k, l, typesize, bytes, psize;

  char method[50];

//! just some statements to prevent compiler warnings of unused variables      
//! those parameter are included for future features
  *correct_flag = 0;
  *ptypesize = 0;
//typesize = filehandle_debug

  array = malloc( DIM1 * DIM2 * DIM3 * sizeof(double) );

  MPI_Comm_rank( local_communicator, &myrank );

  base = myrank * DIM1 * DIM2 * DIM3 + 1;
  utilities_fill_unique_array_3D_double( &array[0], DIM1, DIM2, DIM3, base );

  if ( myrank == 0 ) {
    snprintf( method, 50, "manual" );

    MPI_Type_size( MPI_DOUBLE, &typesize );
    bytes = (DIM1-2) * (DIM2-2) * typesize;
         
    timing_init( testname, &method[0], bytes );
  }

  for( i=0 ; i<outer_loop ; i++ ) {

    psize = (DIM1-2) * (DIM2-2);
    buffer = malloc( psize * sizeof(double) );

    if ( myrank == 0 ) {
      timing_record(1);
    }

    for( j=0 ; j<inner_loop ; j++ ) {

      if ( myrank == 0 ) {
        base = 0;
        for( k=1 ; k<DIM2-1 ; k++ ) {
          for( l=1 ; l<DIM1-1 ; l++ ) {
            buffer[base++] = array[idx3D(l,k,1,DIM1,DIM2)];
          }
        }
        timing_record(2);
        MPI_Send( &buffer[0], psize, MPI_DOUBLE, 1, itag, local_communicator );
        MPI_Recv( &buffer[0], psize, MPI_DOUBLE, 1, itag, local_communicator, MPI_STATUS_IGNORE );
        timing_record(3);
        base = 0;
        for( k=1 ; k<DIM2-1 ; k++ ) {
          for( l=1 ; l<DIM1-1 ; l++ ) {
            array[idx3D(l,k,0,DIM1,DIM2)] = buffer[base++];
          }
        }
        timing_record(4);
      } else {
        MPI_Recv( &buffer[0], psize, MPI_DOUBLE, 0, itag, local_communicator, MPI_STATUS_IGNORE );
        base = 0;
        for( k=1 ; k<DIM2-1 ; k++ ) {
          for( l=1 ; l<DIM1-1 ; l++ ) {
            array[idx3D(l,k,0,DIM1,DIM2)] = buffer[base++];
          }
        }
        base = 0;
        for( k=1 ; k<DIM2-1 ; k++ ) {
          for( l=1 ; l<DIM1-1 ; l++ ) {
            buffer[base++] = array[idx3D(l,k,1,DIM1,DIM2)];
          }
        }
        MPI_Send( &buffer[0], psize, MPI_DOUBLE, 0, itag, local_communicator );
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

  free(array);
}
      
void timing_nas_mg_z_mpi_pack_ddt( int DIM1, int DIM2, int DIM3, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator ) {

  double* array;
  double* buffer;

  int myrank;
  int base, i, j, typesize, bytes, pos;

  MPI_Datatype dtype_face_z_t;

  char method[50];

//! just some statements to prevent compiler warnings of unused variables      
//! those parameter are included for future features
  *correct_flag = 0;
  *ptypesize = 0;
//  typesize = filehandle_debug

  array = malloc( DIM1 * DIM2 * DIM3 * sizeof(double) );  

  MPI_Comm_rank( local_communicator, &myrank );

  base = myrank * DIM1 * DIM2 * DIM3 + 1;
  utilities_fill_unique_array_3D_double( array, DIM1, DIM2, DIM3, base );

  if ( myrank == 0 ) {
    snprintf(method, 50, "mpi_pack_ddt");

    MPI_Type_size( MPI_DOUBLE, &typesize );
    bytes = (DIM1-2) * (DIM2-2) * typesize;

    timing_init( testname, &method[0], bytes );
  }

  for( i=0 ; i<outer_loop ; i++ ) {

    buffer = malloc( (DIM1-2) * (DIM2-2) * sizeof(double) );
    MPI_Type_size( MPI_DOUBLE, &typesize );
    bytes = (DIM1-2) * (DIM2-2) * typesize;
  
    MPI_Type_vector( DIM2-2, DIM1-2, DIM1, MPI_DOUBLE, &dtype_face_z_t );
    MPI_Type_commit( &dtype_face_z_t );

    if ( myrank == 0 ) {
      timing_record(1);
    }

    for( j=0 ; j<inner_loop ; j++ ) {

      if ( myrank == 0 ) {
        pos = 0;
        MPI_Pack( &array[idx3D(1,1,1,DIM1,DIM2)], 1, dtype_face_z_t, &buffer[0], bytes, &pos, local_communicator );
        timing_record(2);
        MPI_Send( &buffer[0], pos, MPI_PACKED, 1, itag, local_communicator );
        MPI_Recv( &buffer[0], bytes, MPI_PACKED, 1, itag, local_communicator, MPI_STATUS_IGNORE );
        timing_record(3);
        pos = 0;
        MPI_Unpack( &buffer[0], bytes, &pos, &array[idx3D(1,1,0,DIM1,DIM2)], 1, dtype_face_z_t, local_communicator );
        timing_record(4);
      } else {
        MPI_Recv( &buffer[0], bytes, MPI_PACKED, 0, itag, local_communicator, MPI_STATUS_IGNORE );
        pos = 0;
        MPI_Unpack( &buffer[0], bytes, &pos, &array[idx3D(1,1,0,DIM1,DIM2)], 1, dtype_face_z_t, local_communicator );
        pos = 0;
        MPI_Pack( &array[idx3D(1,1,1,DIM1,DIM2)], 1, dtype_face_z_t, &buffer[0], bytes, &pos, local_communicator );
        MPI_Send( &buffer[0], pos, MPI_PACKED, 0, itag, local_communicator );
      }

    } //! inner loop

    free( buffer );

    MPI_Type_free( &dtype_face_z_t );

    if ( myrank == 0 ) {
      timing_record(5);
    }

  }

  if ( myrank == 0 ) {
    timing_print( 1 );
  }

  free(array);
}
