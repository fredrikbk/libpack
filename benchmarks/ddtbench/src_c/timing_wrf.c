#include <stdio.h>
#include <stdlib.h>

#include "mpi.h"

#include "ddtbench.h"

#define itag 0

static inline int idx2D(int x, int y, int DIM1) {
  return x+DIM1*y;
}
static inline int idx3D(int x, int y, int z, int DIM1, int DIM2) {
  return x+DIM1*(y+z*DIM2);
}
static inline int idx4D(int x, int y, int z, int t, int DIM1, int DIM2, int DIM3) {
  return x+DIM1*(y+DIM2*(z+DIM3*t));
}

void timing_wrf_manual ( int number_2D, int number_3D, int number_4D, int ims, int ime, int jms, int jme, int kms, int kme, int* limit_4D_arrays, int is, int ie, int js, int je, 
  int ks, int ke, int param_first_scalar, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator ) {

  float** array2Ds;
  float** array3Ds;
  float** array4Ds;

  float* buffer;

  int counter, i, j, bytes, typesize, base;
  int k, l, m, n, o;
  int element_number;
  int myrank;
  int dim1, dim2, dim3;
  int sub_dim1, sub_dim2, sub_dim3;

  char method[50];

//! just some statements to prevent compiler warnings of unused variables
//! those parameter are included for future features
  *correct_flag = 0;
  *ptypesize = 0;
//  typesize = filehandle_debug

  MPI_Comm_rank( local_communicator, &myrank );

//! some conversion from fortran to c
  dim1 = ime-ims+1;
  dim2 = kme-kms+1;
  dim3 = jme-jms+1;

  sub_dim1 = ie-is+1;
  sub_dim2 = ke-ks+1;
  sub_dim3 = je-js+1;

  is = is - ims;
  ks = ks - kms;
  js = js - jms;

  ie = ie - ims;
  ke = ke - kms;
  je = je - jms;

  param_first_scalar--;

//! ================= initialize the arrays =================

//! allocate all needed arrays first

  array2Ds = malloc( number_2D * sizeof(float*) );
  array3Ds = malloc( number_3D * sizeof(float*) );
  array4Ds = malloc( number_4D * sizeof(float*) );

//! allocate and initialize the arrays
//! compute the number of elements in the arrays
  counter = ( number_2D + number_3D * dim2 ) * dim1 * dim3 ;
  for( m=0 ; m<number_4D ; m++ ) {
    counter = counter + limit_4D_arrays[m] * dim1 * dim2 * dim3;
  }
  base = myrank * counter + 1;

  for( m=0 ; m<number_2D ; m++ ) {
    array2Ds[m] = malloc( dim1 * dim3 * sizeof(float) );
    utilities_fill_unique_array_2D_float( array2Ds[m], dim1, dim3, base );
    base = base + dim1 * dim3;
  }

  for( m=0 ; m<number_3D ; m++ ) {
    array3Ds[m] = malloc( dim1 * dim2 * dim3 * sizeof(float) );
    utilities_fill_unique_array_3D_float( array3Ds[m], dim1, dim2, dim3, base );
    base = base + dim1 * dim2 * dim3;
  }

  for( m=0 ; m<number_4D ; m++ ) {
    array4Ds[m] = malloc( dim1 * dim2 * dim3 * limit_4D_arrays[m] * sizeof(float) );
    utilities_fill_unique_array_4D_float( array4Ds[m], dim1, dim2, dim3, limit_4D_arrays[m], base );
    base = base + limit_4D_arrays[m] * dim1 * dim2  * dim3;
  }

  if ( myrank == 0 ) {
    snprintf( &method[0], 50, "manual" );

//! compute the number of bytes to be communicated
//! first compute the number of elements in the subarrays
    counter = number_2D * sub_dim1 * sub_dim3 + number_3D * sub_dim1 * sub_dim2 * sub_dim3;
    for( m=0 ; m<number_4D ; m++ ) {
      if (limit_4D_arrays[m] > param_first_scalar) {
        counter = counter + (limit_4D_arrays[m]-param_first_scalar) * sub_dim1 * sub_dim2 * sub_dim3;
      }
    }
    MPI_Type_size( MPI_FLOAT, &typesize );
    bytes = counter * typesize;

    timing_init( testname, method, bytes );
  }

  for( i=0 ; i<outer_loop ; i++ ) {

//! compute the number of elements in the subarray
    element_number = number_2D * sub_dim1 * sub_dim3 + number_3D * sub_dim1 * sub_dim2 * sub_dim3;
    for( m=0 ; m<number_4D ; m++ ) {
      if (limit_4D_arrays[m] > param_first_scalar) {
        element_number = element_number + (limit_4D_arrays[m]-param_first_scalar) * sub_dim1 * sub_dim2 * sub_dim3;
      }
    }

    buffer = malloc( element_number * sizeof(float) );

    if ( myrank == 0 ) {
      timing_record(1);
    }

    for( j=0 ; j<inner_loop ; j++ ) {

//! =============== ping pong communication =================

//! send the data from rank 0 to rank 1 
      if ( myrank == 0 ) {
//! ==================== pack the data ======================
        counter = 0;
        for( m=0 ; m<number_2D ; m++ ) {
          for( k=js ; k<=je ; k++ ) {
            for( l=is ; l<=ie ; l++ ) {
              buffer[counter++] = *(array2Ds[m]+idx2D(l,k,dim1));
            }
          }
        }
        for( m=0 ; m<number_3D ; m++ ) {
          for( k=js ; k<=je ; k++ ) {
            for( l=ks ; l<=ke ; l++ ) {
              for( n=is ; n<=ie ; n++ ) {
                buffer[counter++] = *(array3Ds[m]+idx3D(n,l,k,dim1,dim2));
              }
            }
          }
        }
        for( m=0 ; m<number_4D ; m++ ) {
          for( k=param_first_scalar ; k<limit_4D_arrays[m] ; k++) {
            for( l=js ; l<=je ; l++ ) {
              for( n=ks ; n<=ke ; n++ ) {
                for( o=is ; o<=ie ; o++ ) {
                  buffer[counter++] = *(array4Ds[m]+idx4D(o,n,l,k,dim1,dim2,dim3));
                }
              }
            }
          }
        }
        timing_record(2);
        MPI_Send( &buffer[0], element_number, MPI_FLOAT, 1, itag, local_communicator );
//! receive the data back from rank 1
        MPI_Recv( &buffer[0], element_number, MPI_FLOAT, 1, itag, local_communicator, MPI_STATUS_IGNORE );
        timing_record(3);
//! =================== unpack the data =====================
        counter = 0;
        for( m=0 ; m<number_2D ; m++ ) {
          for( k=js ; k<=je ; k++ ) {
            for( l=is ; l<=ie ; l++ ) {
              *(array2Ds[m]+idx2D(l,k,dim1)) = buffer[counter++];
            }
          }
        }
        for( m=0 ; m<number_3D ; m++ ) {
          for( k=js ; k<=je ; k++ ) {
            for( l=ks ; l<=ke ; l++ ) {
              for( n=is ; n<=ie ; n++ ) {
                *(array3Ds[m]+idx3D(n,l,k,dim1,dim2)) = buffer[counter++];
              }
            }
          }
        }
        for( m=0 ; m<number_4D ; m++ ) {
          for( k=param_first_scalar ; k<limit_4D_arrays[m] ; k++) {
            for( l=js ; l<=je ; l++ ) {
              for( n=ks ; n<=ke ; n++ ) {
                for( o=is ; o<=ie ; o++ ) {
                  *(array4Ds[m]+idx4D(o,n,l,k,dim1,dim2,dim3)) = buffer[counter++];
                }
              }
            }
          }
        }
        timing_record(4);
//! now for rank 1
      } else {
//! receive from rank 0      
        MPI_Recv( &buffer[0], element_number, MPI_FLOAT, 0, itag, local_communicator, MPI_STATUS_IGNORE );
//! unpack the data 
        counter = 0;
        for( m=0 ; m<number_2D ; m++ ) {
          for( k=js ; k<=je ; k++ ) {
            for( l=is ; l<=ie ; l++ ) {
              *(array2Ds[m]+idx2D(l,k,dim1)) = buffer[counter++];
            }
          }
        }
        for( m=0 ; m<number_3D ; m++ ) {
          for( k=js ; k<=je ; k++ ) {
            for( l=ks ; l<=ke ; l++ ) {
              for( n=is ; n<=ie ; n++ ) {
                *(array3Ds[m]+idx3D(n,l,k,dim1,dim2)) = buffer[counter++];
              }
            }
          }
        }
        for( m=0 ; m<number_4D ; m++ ) {
          for( k=param_first_scalar ; k<limit_4D_arrays[m] ; k++) {
            for( l=js ; l<=je ; l++ ) {
              for( n=ks ; n<=ke ; n++ ) {
                for( o=is ; o<=ie ; o++ ) {
                  *(array4Ds[m]+idx4D(o,n,l,k,dim1,dim2,dim3)) = buffer[counter++];
                }
              }
            }
          }
        }
//! pack the data
        counter = 0;
        for( m=0 ; m<number_2D ; m++ ) {
          for( k=js ; k<=je ; k++ ) {
            for( l=is ; l<=ie ; l++ ) {
              buffer[counter++] = *(array2Ds[m]+idx2D(l,k,dim1));
            }
          }
        }
        for( m=0 ; m<number_3D ; m++ ) {
          for( k=js ; k<=je ; k++ ) {
            for( l=ks ; l<=ke ; l++ ) {
              for( n=is ; n<=ie ; n++ ) {
                buffer[counter++] = *(array3Ds[m]+idx3D(n,l,k,dim1,dim2));
              }
            }
          }
        }
        for( m=0 ; m<number_4D ; m++ ) {
          for( k=param_first_scalar ; k<limit_4D_arrays[m] ; k++) {
            for( l=js ; l<=je ; l++ ) {
              for( n=ks ; n<=ke ; n++ ) {
                for( o=is ; o<=ie ; o++ ) {
                  buffer[counter++] = *(array4Ds[m]+idx4D(o,n,l,k,dim1,dim2,dim3));
                }
              }
            }
          }
        }

//!> send to rank 0
        MPI_Send( buffer, element_number, MPI_FLOAT, 0, itag, local_communicator );
      } //! of myrank .EQ. 0?

    } //! inner_loop

//! ======================= clean up ========================

    free( buffer );

    if ( myrank == 0 ) {
      timing_record(5);
    }

  } //! outer_loop

  if ( myrank == 0 ) {
    timing_print( 1 );
  }

//! ======================= clean up ========================

  for( m=0 ; m<number_2D ; m++ ) {
    free( array2Ds[m] );
  }

  for( m=0 ; m<number_3D ; m++ ) {
    free( array3Ds[m] );
  }

  for( m=0 ; m<number_4D ; m++ ) {
    free( array4Ds[m] );
  }
     
  free( array2Ds );
  free( array3Ds );
  free( array4Ds );
}

void timing_wrf_vec_ddt ( int number_2D, int number_3D, int number_4D, int ims, int ime, int jms, int jme, int kms, int kme, int* limit_4D_arrays, int is, int ie, int js, int je, 
  int ks, int ke, int param_first_scalar, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator ) {

  float** array2Ds;
  float** array3Ds;
  float** array4Ds;

  int m, counter, i, j, bytes, base;
  int myrank;
  int dim1, dim2, dim3;
  int sub_dim1, sub_dim2, sub_dim3;

  char method[50];

//! variables for the MPI derived datatypes
  MPI_Datatype dtype_subarray_t, dtype_temp_2D_t, dtype_temp_t, dtype_temp_3D_t;
  MPI_Datatype* oldtype;
  int* blocklength;
  MPI_Aint* displacement;
  MPI_Aint stride;
  int typesize;

//! just some statements to prevent compiler warnings of unused variables
//! those parameter are included for future features
  *correct_flag = 0;
  *ptypesize = 0;
//  typesize = filehandle_debug

//! some conversion from fortran to c
  dim1 = ime-ims+1;
  dim2 = kme-kms+1;
  dim3 = jme-jms+1;

  sub_dim1 = ie-is+1;
  sub_dim2 = ke-ks+1;
  sub_dim3 = je-js+1;

  is = is - ims;
  ks = ks - kms;
  js = js - jms;

  ie = ie - ims;
  ke = ke - kms;
  je = je - jms;

  param_first_scalar--;

  MPI_Comm_rank( local_communicator, &myrank );

//! ================= initialize the arrays =================

//! allocate all needed arrays first
  array2Ds = malloc( number_2D * sizeof(float*) );
  array3Ds = malloc( number_3D * sizeof(float*) );
  array4Ds = malloc( number_4D * sizeof(float*) );

//! allocate and initialize the arrays
//! compute the number of elements in the arrays
  counter = ( number_2D + number_3D * dim2 ) * dim1 * dim3;
  for( m=0 ; m<number_4D ; m++ ) {
    counter = counter + limit_4D_arrays[m] * dim1 * dim2 * dim3;
  }
  base = myrank * counter + 1;
 
  for( m=0 ; m<number_2D ; m++ ) {
    array2Ds[m] = malloc( dim1 * dim3 * sizeof(float) );
    utilities_fill_unique_array_2D_float( array2Ds[m], dim1, dim3, base );
    base = base + dim1 * dim3;
  }

  for( m=0 ; m<number_3D ; m++ ) {
    array3Ds[m] = malloc( dim1 * dim2 * dim3 * sizeof(float) );
    utilities_fill_unique_array_3D_float( array3Ds[m], dim1, dim2, dim3, base );
    base = base + dim1 * dim2 * dim3;
  }

  for( m=0 ; m<number_4D ; m++ ) {
    array4Ds[m] = malloc( dim1 * dim2 * dim3 * limit_4D_arrays[m] * sizeof(float) );
    utilities_fill_unique_array_4D_float( array4Ds[m], dim1, dim2, dim3, limit_4D_arrays[m], base );
    base = base + limit_4D_arrays[m] * dim1 * dim2 * dim3;
  }

  if ( myrank == 0 ) {
    snprintf( &method[0], 50, "mpi_ddt" );

//! compute the number of bytes to be communicated
//! first compute the number of elements in the subarrays

    counter = number_2D * sub_dim1 * sub_dim3 + number_3D * sub_dim1 * sub_dim2 * sub_dim3;
    for( m=0 ; m<number_4D ; m++ ) {
      if (limit_4D_arrays[m] > param_first_scalar) {
        counter = counter + (limit_4D_arrays[m]-param_first_scalar) * sub_dim1 * sub_dim2 * sub_dim3;
      }
    }
    MPI_Type_size( MPI_FLOAT, &typesize );
    bytes = counter * typesize;

    timing_init( testname, &method[0], bytes );
  }

  for( i=0 ; i<outer_loop ; i++ ) {

//! ========== building the MPI derived datatype ============

//! compute the number of subarrays, so that the buffer for the datatype
//! creation are big enough

    counter = number_2D + number_3D;
    for( m=0 ; m<number_4D ; m++ ) {
      if ( limit_4D_arrays[m] >= param_first_scalar ) {
        counter++;
      }
    }
 
    oldtype = malloc( counter * sizeof(MPI_Datatype) );
    blocklength = malloc( counter * sizeof(int) );
    displacement = malloc( counter * sizeof(MPI_Aint) );

    for( j=0 ; j<counter ; j++) {
      blocklength[j] = 1;
    }
  
    counter = 0;
    
//! building the 2D datatype    
    MPI_Type_vector( sub_dim3, sub_dim1, dim1, MPI_FLOAT, &dtype_temp_2D_t );

//! set the parameter for the struct datatype for the 2D parts
    for( m=0 ; m<number_2D ; m++ ) {
      MPI_Get_address( array2Ds[m]+idx2D(is,js,dim1), &displacement[counter] );
      oldtype[counter++] = dtype_temp_2D_t;
    }

//! building the 3D datatype
    MPI_Type_size (MPI_FLOAT, &typesize);
    MPI_Type_vector( sub_dim2, sub_dim1, dim1, MPI_FLOAT, &dtype_temp_t );
    stride = dim1 * dim2 * typesize;
    MPI_Type_create_hvector( sub_dim3, 1, stride, dtype_temp_t, &dtype_temp_3D_t );
    MPI_Type_free( &dtype_temp_t );

//! set the parameter for the struct datatype for the 3D parts
    for( m=0 ; m<number_3D ; m++ ) {
      MPI_Get_address( array3Ds[m]+idx3D(is,ks,js,dim1,dim2), &displacement[counter] );
      oldtype[counter++] = dtype_temp_3D_t;
    }

//! building the 4D datatypes for each 4D array and setting the parameters for the struct datatype
    stride = stride * dim3;
    for( m=0 ; m<number_4D ; m++ ) {
      if (limit_4D_arrays[m] > param_first_scalar) {
        MPI_Get_address( array4Ds[m]+idx4D(is,ks,js,param_first_scalar,dim1,dim2,dim3), &displacement[counter] );
        MPI_Type_create_hvector( limit_4D_arrays[m]-param_first_scalar, 1, stride, dtype_temp_3D_t, &oldtype[counter++] );
      }
    }

//! create the datatype which is really used
    MPI_Type_create_struct( counter, &blocklength[0], &displacement[0], &oldtype[0], &dtype_subarray_t );
    MPI_Type_commit( &dtype_subarray_t );

//! some cleanup
    MPI_Type_free( &dtype_temp_2D_t );
    MPI_Type_free( &dtype_temp_3D_t );
    counter = number_2D + number_3D;
    for( m=0 ; m<number_4D ; m++ ) {
      if (limit_4D_arrays[m] > param_first_scalar) {
        MPI_Type_free( &oldtype[counter++] );
      }
    }

    free( oldtype );
    free( blocklength );
    free( displacement );

    if ( myrank == 0 ) {
      timing_record(1);
    }

    for( j=0 ; j<inner_loop ; j++ ) {

//! =============== ping pong communication =================

//! send the data from rank 0 to rank 1 
      if ( myrank == 0 ) {
        MPI_Send( MPI_BOTTOM, 1, dtype_subarray_t, 1, itag, local_communicator );
//! receive the data back from rank 1
        MPI_Recv( MPI_BOTTOM, 1, dtype_subarray_t, 1, itag, local_communicator, MPI_STATUS_IGNORE );
        timing_record(3);
//! now for rank 1
      } else {
//! receive from rank 0      
        MPI_Recv( MPI_BOTTOM, 1, dtype_subarray_t, 0, itag, local_communicator, MPI_STATUS_IGNORE );
//! send to rank 0
        MPI_Send( MPI_BOTTOM, 1, dtype_subarray_t, 0, itag, local_communicator );
      }

    } //! inner_loop

//! ======================= clean up ========================

    MPI_Type_free( &dtype_subarray_t );

    if ( myrank == 0 ) {
      timing_record(5);
    }

  } //! outer_loop

  if ( myrank == 0 ) {
    timing_print( 1 );
  }

//! ======================= clean up ========================

  for( m=0 ; m<number_2D ; m++ ) {
    free( array2Ds[m] );
  }

  for( m=0 ; m<number_3D ; m++ ) {
    free( array3Ds[m] );
  }

  for( m=0 ; m<number_4D ; m++ ) {
    free( array4Ds[m] );
  }
     
  free( array2Ds );
  free( array3Ds );
  free( array4Ds );
}

void timing_wrf_vec_mpi_pack_ddt ( int number_2D, int number_3D, int number_4D, int ims, int ime, int jms, int jme, int kms, int kme, int* limit_4D_arrays, int is, int ie, int js, 
  int je, int ks, int ke, int param_first_scalar, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator ) {

  float** array2Ds;
  float** array3Ds;
  float** array4Ds;

  float* buffer;

  int m, counter, i, j, bytes, typesize, pos, base;
  int element_number;
  int myrank;
  int dim1, dim2, dim3;
  int sub_dim1, sub_dim2, sub_dim3;

//! some variables for writing output
  char method[50];

//! variables for the MPI derived datatypes
  MPI_Datatype dtype_subarray_t, dtype_temp_2D_t, dtype_temp_t, dtype_temp_3D_t;
  MPI_Datatype* oldtype;
  int* blocklength;
  MPI_Aint* displacement;
  MPI_Aint stride;

//! just some statements to prevent compiler warnings of unused variables
//! those parameter are included for future features
  *correct_flag = 0;
  *ptypesize = 0;
//  typesize = filehandle_debug;

  MPI_Comm_rank( local_communicator, &myrank );

//! some conversion from fortran to c
  dim1 = ime-ims+1;
  dim2 = kme-kms+1;
  dim3 = jme-jms+1;

  sub_dim1 = ie-is+1;
  sub_dim2 = ke-ks+1;
  sub_dim3 = je-js+1;

  is = is - ims;
  ks = ks - kms;
  js = js - jms;

  ie = ie - ims;
  ke = ke - kms;
  je = je - jms;

  param_first_scalar--;

//! ================= initialize the arrays =================

//! allocate all needed arrays first
  array2Ds = malloc( number_2D * sizeof(float*) );
  array3Ds = malloc( number_3D * sizeof(float*) );
  array4Ds = malloc( number_4D * sizeof(float*) );

//! allocate and initialize the arrays
//! compute the number of elements in the arrays
  counter = ( number_2D + number_3D * dim2 ) * dim1 * dim3;
  for( m=0 ; m<number_4D ; m++ ) {
    counter = counter + limit_4D_arrays[m] * dim1 * dim2 * dim3;
  }
  base = myrank * counter + 1;

  for( m=0 ; m<number_2D ; m++ ) {
    array2Ds[m] = malloc( dim1 * dim3 * sizeof(float) );
    utilities_fill_unique_array_2D_float( array2Ds[m], dim1, dim3, base );
    base = base + dim1 * dim3;
  }

  for( m=0 ; m<number_3D ; m++ ) {
    array3Ds[m] = malloc( dim1 * dim2 * dim3 * sizeof(float) );
    utilities_fill_unique_array_3D_float( array3Ds[m], dim1, dim2, dim3, base );
    base = base + dim1 * dim2 * dim3 ;
  }

  for( m=0 ; m<number_4D ; m++ ) {
    array4Ds[m] = malloc( dim1 * dim2 * dim3 * limit_4D_arrays[m] * sizeof(float) );
    utilities_fill_unique_array_4D_float( array4Ds[m], dim1, dim2, dim3, limit_4D_arrays[m], base );
    base = base + limit_4D_arrays[m] * dim1 * dim2 * dim3;
  }

  if ( myrank == 0 ) {
    snprintf( &method[0], 50, "mpi_pack_ddt" );

//! compute the number of bytes to be communicated
//! first compute the number of elements in the subarrays

    counter = number_2D * sub_dim1 * sub_dim3 + number_3D * sub_dim1 * sub_dim2 * sub_dim3;
    for( m=0 ; m<number_4D ; m++ ) {
      if (limit_4D_arrays[m] > param_first_scalar) {
        counter = counter + (limit_4D_arrays[m]-param_first_scalar) * sub_dim1 * sub_dim2 * sub_dim3;
      }
    }
        
    MPI_Type_size( MPI_FLOAT, &typesize );
    bytes = counter * typesize;

    timing_init( testname, &method[0], bytes );
  }

  for( i=0 ; i<outer_loop ; i++ ) {

//! compute the number of elements in the subarray
    element_number = number_2D * sub_dim1 * sub_dim3 + number_3D * sub_dim1 * sub_dim2 * sub_dim3;
    for( m=0 ; m<number_4D ; m++ ) {
      if (limit_4D_arrays[m] > param_first_scalar) {
        element_number = element_number + (limit_4D_arrays[m]-param_first_scalar) * sub_dim1 * sub_dim2 * sub_dim3;
      }
    }

    buffer = malloc( element_number * sizeof(float) );
    MPI_Type_size( MPI_FLOAT, &typesize );
    bytes = element_number * typesize;

//! compute the number of arrays, so that the buffer for the datatype creation are big enough
    counter = number_2D + number_3D;
    for( m=0 ; m<number_4D ; m++ ) {
      if ( limit_4D_arrays[m] > param_first_scalar ) {
        counter++;
      }
    }
 
    oldtype = malloc( counter * sizeof(MPI_Datatype) );
    blocklength = malloc( counter * sizeof(int) );
    displacement = malloc( counter * sizeof(MPI_Aint) );

    for( j=0 ; j<counter ; j++ ) {
      blocklength[j] = 1;
    }
  
    MPI_Type_vector( sub_dim3, sub_dim1, dim1, MPI_FLOAT, &dtype_temp_2D_t );
 
    counter = 0;
    for( m=0 ; m<number_2D ; m++ ) {
      MPI_Get_address( array2Ds[m]+idx2D(is,js,dim1), &displacement[counter] );
      oldtype[counter++] = dtype_temp_2D_t;
    }

    MPI_Type_size( MPI_FLOAT, &typesize );
    MPI_Type_vector( sub_dim2, sub_dim1, dim1, MPI_FLOAT, &dtype_temp_t );
    stride = dim1 * dim2 * typesize;
    MPI_Type_create_hvector( sub_dim3, 1, stride, dtype_temp_t, &dtype_temp_3D_t );
    MPI_Type_free( &dtype_temp_t );

    for( m=0 ; m<number_3D ; m++ ) {
      MPI_Get_address( array3Ds[m]+idx3D(is,ks,js,dim1,dim2), &displacement[counter] );
      oldtype[counter++] = dtype_temp_3D_t;
    }

    stride = stride * dim3;
    for( m=0 ; m<number_4D ; m++ ) {
      if (limit_4D_arrays[m] > param_first_scalar) {
        MPI_Get_address( array4Ds[m]+idx4D(is,ks,js,param_first_scalar,dim1,dim2,dim3), &displacement[counter] );
        MPI_Type_create_hvector( limit_4D_arrays[m]-param_first_scalar, 1, stride, dtype_temp_3D_t, &oldtype[counter++] );
      }
    }

    MPI_Type_create_struct( counter, &blocklength[0], &displacement[0], &oldtype[0], &dtype_subarray_t );
    MPI_Type_commit( &dtype_subarray_t );

    MPI_Type_free( &dtype_temp_2D_t );
    MPI_Type_free( &dtype_temp_3D_t );
    counter = number_2D + number_3D;
    for( m=0 ; m<number_4D ; m++ ) {
      if (limit_4D_arrays[m] > param_first_scalar) {
        MPI_Type_free( &oldtype[counter++] );
      }
    }

    free( oldtype );
    free( blocklength );
    free( displacement );

    if ( myrank == 0 ) {
      timing_record(1);
    }

    for( j=0 ; j<inner_loop ; j++ ) {

//! =============== ping pong communication =================

//! send the data from rank 0 to rank 1 
      if ( myrank == 0 ) {
//! ==================== pack the data ======================
        pos = 0;
        MPI_Pack( MPI_BOTTOM, 1, dtype_subarray_t, &buffer[0], bytes, &pos, local_communicator );
        timing_record(2);
        MPI_Send( &buffer[0], pos, MPI_PACKED, 1, itag, local_communicator );
//! receive the data back from rank 1
        MPI_Recv( &buffer[0], bytes, MPI_PACKED, 1, itag, local_communicator, MPI_STATUS_IGNORE );
        timing_record(3);
//! =================== unpack the data =====================
        pos = 0;
        MPI_Unpack( &buffer[0], bytes, &pos, MPI_BOTTOM, 1, dtype_subarray_t, local_communicator );
        timing_record(4);
//! now for rank 1
      } else {
//! receive from rank 0      
        MPI_Recv( &buffer[0], element_number, MPI_FLOAT, 0, itag, local_communicator, MPI_STATUS_IGNORE );
//! unpack the data 
        pos = 0;
        MPI_Unpack( &buffer[0], bytes, &pos, MPI_BOTTOM, 1, dtype_subarray_t, local_communicator );
//! pack the data
        pos = 0;
        MPI_Pack( MPI_BOTTOM, 1, dtype_subarray_t, &buffer[0], bytes, &pos, local_communicator );
//! send to rank 0
        MPI_Send( &buffer[0], element_number, MPI_FLOAT, 0, itag, local_communicator );
      } //! of myrank .EQ. 0?

    } //! inner_loop

//! ======================= clean up ========================

    MPI_Type_free( &dtype_subarray_t );

    free( buffer );

    if ( myrank == 0 ) {
      timing_record(5);
    }

  } //! outer_loop

  if ( myrank == 0 ) {
    timing_print( 1 );
  }

//!> ======================= clean up ========================

  for( m=0 ; m<number_2D ; m++ ) {
    free( array2Ds[m] );
  }

  for( m=0 ; m<number_3D ; m++ ) {
    free( array3Ds[m] );
  }

  for( m=0 ; m<number_4D ; m++ ) {
    free( array4Ds[m] );
  }
     
  free( array2Ds );
  free( array3Ds );
  free( array4Ds );
}

#if 0
void timing_wrf_sa_ddt ( int number_2D, int number_3D, int number_4D, int ims, int ime, int jms, int jme, int kms, int kme, int* limit_4D_arrays, int is, int ie, int js, int je, 
  int ks, int ke, int param_first_scalar, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator ) {

  float** array2Ds;
  float** array3Ds;
  float** array4Ds;

  int m, counter, i, j, bytes, base;
  int myrank;
  int dim1, dim2, dim3;
  int sub_dim1, sub_dim2, sub_dim3;

//! some variables for writing output
  char method[50];

//! variables for the MPI derived datatype
  MPI_Datatype dtype_subarray_t, dtype_temp_2D_t, dtype_temp_3D_t;
  MPI_Datatype* oldtype;
  int* blocklength;
  MPI_Aint* displacement;
  int typesize;
  int arraysize[4], subarraysize[4], subarraystart[4];

//! just some statements to prevent compiler warnings of unused variables
//! those parameter are included for future features
  correct_flag = 0;
  ptypesize = 0;
//  typesize = filehandle_debug

  MPI_Comm_rank( local_communicator, &myrank );

//! some conversion from fortran to c
  dim1 = ime-ims+1;
  dim2 = kme-kms+1;
  dim3 = jme-jms+1;

  sub_dim1 = ie-is+1;
  sub_dim2 = ke-ks+1;
  sub_dim3 = je-js+1;

  is = is - ims;
  ks = ks - kms;
  js = js - jms;

  ie = ie - ims;
  ke = ke - kms;
  je = je - jms;

  param_first_scalar--;

//! ================= initialize the arrays =================

//! allocate all needed arrays first
  array2Ds = malloc( number_2D * sizeof(float*) );
  array3Ds = malloc( number_3D * sizeof(float*) );
  array4Ds = malloc( number_4D * sizeof(float*) );

//! allocate and initialize the arrays
//! compute the number of elements in the arrays
  counter = ( number_2D + number_3D * dim2 ) * dim1 * dim3;
  for( m=0 ; m<number_4D ; m++ ) {
    counter = counter + limit_4D_arrays[m] * dim1 * dim2 * dim3;
  }
  base = myrank * counter + 1;

  for( m=0 ; m<number_2D ; m++ ) {
    array2Ds[m] = malloc( dim1 * dim3 * sizeof(float) );
    utilities_fill_unique_array_2D_float( array2Ds[m], dim1, dim3, base );
    base = base + dim1 * dim3;
  }

  for( m=0 ; m<number_3D ; m++ ) {
    array3Ds[m] = malloc( dim1 * dim2 * dim3 * sizeof(float) );
    utilities_fill_unique_array_3D_float( array3Ds[m], dim1, dim2, dim3, base );
    base = base + dim1 * dim2 * dim3;
  }

  for( m=0 ; m<number_4D ; m++ ) {
    array4Ds[m] = malloc( dim1 * dim2 * dim3 * limit_4D_arrays[m] * sizeof(float) );
    utilities_fill_unique_array_4D_float( array4Ds[m], dim1, dim2, dim3, limit_4D_arrays[m], base );
    base = base + limit_4D_arrays[m] * dim1 * dim2 * dim3;
  }

  if ( myrank == 0 ) {
    snprintf( &method[0], 50, "mpi_ddt" );

//! compute the number of bytes to be communicated
//! first compute the number of elements in the subarrays

    counter = number_2D * sub_dim1 * sub_dim3 + number_3D * sub_dim1 * sub_dim2 * sub_dim3;
    for( m=0 ; m<number_4D ; m++ ) {
      if (limit_4D_arrays[m] > param_first_scalar) {
        counter = counter + (limit_4D_arrays[m]-param_first_scalar) * sub_dim1 * sub_dim2 * sub_dim3;
      }
    }
    
    MPI_Type_size( MPI_FLOAT, &typesize );
    bytes = counter * typesize;

    timing_init( testname, &method[0], bytes );
  }

  for( i=0 ; i<outer_loop ; i++ ) {

//! ========== building the MPI derived datatype ============

//! compute the number of subarrays, so that the buffer for the datatype
//! creation are big enough
    counter = number_2D + number_3D;
    for( m=0 ; m<number_4D ; m++ ) {
      if ( limit_4D_arrays[m] > param_first_scalar ) {
        counter++;
      }
    }
 
    oldtype = malloc( counter * sizeof(MPI_Datatype) );
    blocklength = malloc( counter * sizeof(int) );
    displacement = malloc( counter * sizeof(MPI_Aint) );

    for( j=0 ; j<counter ; j++ ) {
      blocklength[j] = 1;
    }
  
    counter = 0;
       
//! create the 2D subarray type 
    arraysize[2] = dim3;
    arraysize[3] = dim1;
    subarraysize[2] = sub_dim3;
    subarraysize[3] = sub_dim1;
    subarraystart[2] = js;
    subarraystart[3] = is;
    MPI_Type_create_subarray( 2, &arraysize[2], &subarraysize[2], &subarraystart[2], MPI_ORDER_C, MPI_FLOAT, &dtype_temp_2D_t );

//! create the parameters of the 2D arrays for the struct type
    for( m=0 ; m<number_2D ; m++ ) {
      MPI_Get_address( array2Ds[m], &displacement[counter] );
      oldtype[counter++] = dtype_temp_2D_t;
    }

//! create the 3D subarray type
    arraysize[1] = dim3;
    arraysize[2] = dim2;
//    arraysize[3] = dim1; already set
    subarraysize[1] = sub_dim3;
    subarraysize[2] = sub_dim2;
//    subarraysize[3] = sub_dim1; already set
    subarraystart[1] = js;
    subarraystart[2] = ks;
//    subarraystart[3] = is; already set
    MPI_Type_create_subarray( 3, &arraysize[1], &subarraysize[1], &subarraystart[1], MPI_ORDER_C, MPI_FLOAT, &dtype_temp_3D_t );

    for( m=0 ; m<number_3D ; m++ ) {
      MPI_Get_address( array3Ds[m], &displacement[counter] );
      oldtype[counter++] = dtype_temp_3D_t;
    }

//! create the 4D subarray types
    subarraystart[0] = param_first_scalar;
    for( m=0 ; m<number_4D ; m++ ) {
      if (limit_4D_arrays[m] > param_first_scalar) {
        arraysize[0] = limit_4D_arrays[m];
        subarraysize[0] = limit_4D_arrays[m]-param_first_scalar;
        MPI_Type_create_subarray( 4, &arraysize[0], &subarraysize[0], &subarraystart[0], MPI_ORDER_C, MPI_FLOAT, &oldtype[counter] );
        MPI_Get_address( array4Ds[m], &displacement[counter++] );
      }
    }

//! create the all-embracing struct type
    MPI_Type_create_struct(counter, &blocklength[0], &displacement[0], &oldtype[0], &dtype_subarray_t);
    MPI_Type_commit( &dtype_subarray_t );

    MPI_Type_free( &dtype_temp_2D_t );
    MPI_Type_free( &dtype_temp_3D_t );
    counter = number_2D + number_3D;
    for( m=0 ; m<number_4D ; m++ ) {
      if (limit_4D_arrays[m] > param_first_scalar) {
        MPI_Type_free( &oldtype[counter++] );
      }
    }

    free( oldtype );
    free( blocklength );
    free( displacement );

    if ( myrank == 0 ) {
      timing_record(1);
    }

    for( j=0 ; j<inner_loop ; j++ ) {

//! =============== ping pong communication =================

//! send the data from rank 0 to rank 1 
      if ( myrank == 0 ) {
        MPI_Send( MPI_BOTTOM, 1, dtype_subarray_t, 1, itag, local_communicator );
//! receive the data back from rank 1
        MPI_Recv( MPI_BOTTOM, 1, dtype_subarray_t, 1, itag, local_communicator, MPI_STATUS_IGNORE );
        timing_record(3);
//! now for rank 1
      } else {
//! receive from rank 0      
        MPI_Recv( MPI_BOTTOM, 1, dtype_subarray_t, 0, itag, local_communicator, MPI_STATUS_IGNORE );
//! send to rank 0
        MPI_Send( MPI_BOTTOM, 1, dtype_subarray_t, 0, itag, local_communicator );
      }

    } //! inner_loop

//! ======================= clean up ========================

    MPI_Type_free( &dtype_subarray_t );

    if ( myrank == 0 ) {
      timing_record(5);
    }

  } //! outer_loop

  if ( myrank == 0 ) {
    timing_print( 1 );
  }

//! ======================= clean up ========================

  for( m=0 ; m<number_2D ; m++ ) {
    free( array2Ds[m] );
  }

  for( m=0 ; m<number_3D ; m++ ) {
    free( array3Ds[m] );
  }

  for( m=0 ; m<number_4D ; m++ ) {
    free( array4Ds[m] );
  }
     
  free( array2Ds );
  free( array3Ds );
  free( array4Ds );
}

void timing_wrf_sa_mpi_pack_ddt ( int number_2D, int number_3D, int number_4D, int ims, int ime, int jms, int jme, int kms, int kme, int* limit_4D_arrays, int is, int ie, int js, int je, 
  int ks, int ke, int param_first_scalar, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator ) {

  float** array2Ds;
  float** array3Ds;
  float** array4Ds;
      
  int m, counter, i, j, bytes, typesize, pos, base;
  int element_number;
  int myrank;
  int dim1, dim2, dim3;
  int sub_dim1, sub_dim2, sub_dim3;

  float* buffer;

  char method[50];

//! variables for the MPI derived datatypes
  MPI_Datatype dtype_subarray_t, dtype_temp_2D_t, dtype_temp_3D_t;
  MPI_Datatype* oldtype;
  int* blocklength;
  MPI_Aint* displacement;
  int arraysize[4], subarraysize[4], subarraystart[4];

//! just some statements to prevent compiler warnings of unused variables
//! those parameter are included for future features
  *correct_flag = 0;
  *ptypesize = 0;
//   typesize = filehandle_debug

   MPI_Comm_rank( local_communicator, &myrank );

//! some conversion from fortran to c
  dim1 = ime-ims+1;
  dim2 = kme-kms+1;
  dim3 = jme-jms+1;

  sub_dim1 = ie-is+1;
  sub_dim2 = ke-ks+1;
  sub_dim3 = je-js+1;

  is = is - ims;
  ks = ks - kms;
  js = js - jms;

  ie = ie - ims;
  ke = ke - kms;
  je = je - jms;

  param_first_scalar--;

//! ================= initialize the arrays =================

//! allocate all needed arrays first
  array2Ds = malloc( number_2D * sizeof(float*) );
  array3Ds = malloc( number_3D * sizeof(float*) );
  array4Ds = malloc( number_4D * sizeof(float*) );
 
//! allocate and initialize the arrays
//! compute the number of elements in the arrays
  counter = ( number_2D + number_3D * dim2 ) * dim1 * dim3;
  for( m=0 ; m<number_4D ; m++ ) {
    counter = counter + limit_4D_arrays[m] * dim1 * dim2 * dim3;
  }
  base = myrank * counter + 1;
 
  for( m=0 ; m<number_2D ; m++ ) {
    array2Ds[m] = malloc( dim1 * dim3 * sizeof(float) );
    utilities_fill_unique_array_2D_float( array2Ds[m], dim1, dim3, base );
    base = base + dim1 * dim3;
  }

  for( m=0 ; m<number_3D ; m++ ) {
    array3Ds[m] = malloc( dim1 * dim2 * dim3 * sizeof(float) );
    utilities_fill_unique_array_3D_float( array3Ds[m], dim1, dim2, dim3, base );
    base = base + dim1 * dim2 * dim3;
  }

  for( m=0 ; m<number_4D ; m++ ) {
    array4Ds[m] = malloc( dim1 * dim2 * dim3 * limit_4D_arrays[m] * sizeof(float) );
    utilities_fill_unique_array_4D_float( array4Ds[m], dim1, dim2, dim3, limit_4D_arrays[m], base );
    base = base + limit_4D_arrays[m] * dim1 * dim2 * dim3;
  }

  if ( myrank == 0 ) {
    snprintf( &method[0], 50, "mpi_pack_ddt" );

//! compute the number of bytes to be communicated
//! first compute the number of elements in the subarrays

    counter = number_2D * sub_dim1 * sub_dim3 + number_3D * sub_dim1 * sub_dim2 * sub_dim3;
    for( m=0 ; m<number_4D ; m++ ) {
      if (limit_4D_arrays[m] > param_first_scalar) {
        counter = counter + (limit_4D_arrays[m]-param_first_scalar) * sub_dim1 * sub_dim2 * sub_dim3;
      }
    }
    MPI_Type_size( MPI_FLOAT, &typesize );
    bytes = counter * typesize;

    timing_init( testname, &method[0], bytes );
  }

  for( i=0 ; i<outer_loop ; i++ ) {

//! compute the number of elements in the subarray
    element_number = number_2D * sub_dim1 * sub_dim3 + number_3D * sub_dim1 * sub_dim2 * sub_dim3;
    for( m=0 ; m<number_4D ; m++ ) {
      if (limit_4D_arrays[m] > param_first_scalar) {
        element_number = element_number + (limit_4D_arrays[m]-param_first_scalar) * sub_dim1 * sub_dim2 * sub_dim3;
      }
    }

    buffer = malloc( element_number * sizeof(float) );
    MPI_Type_size( MPI_FLOAT, &typesize );
    bytes = element_number * typesize;

//! compute the number of subarrays, so that buffer for the datatype
//! creation are big enough
    counter = number_2D + number_3D;
    for( m=0 ; m<number_4D ; m++ ) {
      if ( limit_4D_arrays[m] > param_first_scalar ) {
        counter++;
      }
    }
 
    oldtype = malloc( counter * sizeof(MPI_Datatype) );
    blocklength = malloc( counter * sizeof(int) );
    displacement = malloc( counter * sizeof(MPI_Aint) );

    for( j=0 ; j<counter ; j++ ) {
      blocklength[j] = 1;
    }

    counter = 0;
 
//! create the subtype for the 2D case
    arraysize[2] = dim3;
    arraysize[3] = dim1;
    subarraysize[2] = sub_dim3;
    subarraysize[3] = sub_dim1;
    subarraystart[2] = js;
    subarraystart[3] = is;
    MPI_Type_create_subarray( 2, &arraysize[2], &subarraysize[2], &subarraystart[2], MPI_ORDER_C, MPI_FLOAT, &dtype_temp_2D_t );

    for( m=0 ; m<number_2D ; m++ ) {
      MPI_Get_address( array2Ds[m], &displacement[counter] );
      oldtype[counter++] = dtype_temp_2D_t;
    }

//! create the subtype for the 3D case
    arraysize[1] = dim3;
    arraysize[2] = dim2;
//    arraysize[3] = dim1, already set
    subarraysize[1] = sub_dim3;
    subarraysize[2] = sub_dim2;
//    subarraysize[3] = sub_dim1; already set
    subarraystart[1] = js;
    subarraystart[2] = ks;
//    subarraystart[3] = is;
    MPI_Type_create_subarray( 3, &arraysize[1], &subarraysize[1], &subarraystart[1], MPI_ORDER_C, MPI_FLOAT, &dtype_temp_3D_t );

    for( m=0 ; m<number_3D ; m++ ) {
      MPI_Get_address( array3Ds[m], &displacement[counter] );
      oldtype[counter++] = dtype_temp_3D_t;
    }

//! create the subtypes for the 4D cases
    subarraystart[0] = param_first_scalar;
    for( m=0 ; m<number_4D ; m++ ) {
      if (limit_4D_arrays[m] > param_first_scalar) {
        arraysize[0] = limit_4D_arrays[m];
        subarraysize[0] = limit_4D_arrays[m]-param_first_scalar;
        MPI_Type_create_subarray( 4, &arraysize[0], &subarraysize[0], &subarraystart[0], MPI_ORDER_C, MPI_FLOAT, &oldtype[counter] );
        MPI_Get_address( array4Ds[m], &displacement[counter++] );
      }
    }
 
//! create the all-embracing struct type
    MPI_Type_create_struct( counter, blocklength, displacement, oldtype, &dtype_subarray_t );
    MPI_Type_commit( &dtype_subarray_t );

    MPI_Type_free( &dtype_temp_2D_t );
    MPI_Type_free( &dtype_temp_3D_t );
    counter = number_2D + number_3D;
    for( m=0 ; m<number_4D ; m++ ) {
      if (limit_4D_arrays[m] > param_first_scalar) {
        MPI_Type_free( &oldtype[counter++] );
      }
    }

    free( oldtype );
    free( blocklength );
    free( displacement );

    if ( myrank == 0 ) {
      timing_record(1);
    }

    for( j=0 ; j<inner_loop ; j++ ) {

//! =============== ping pong communication =================

//! send the data from rank 0 to rank 1 
      if ( myrank == 0 ) {
//! ==================== pack the data ======================
        pos = 0;
        MPI_Pack( MPI_BOTTOM, 1, dtype_subarray_t, &buffer[0], bytes, &pos, local_communicator );
        timing_record(2);
        MPI_Send( &buffer[0], pos, MPI_PACKED, 1, itag, local_communicator );
//! receive the data back from rank 1
        MPI_Recv( &buffer[0], bytes, MPI_PACKED, 1, itag, local_communicator, MPI_STATUS_IGNORE );
        timing_record(3);
//! =================== unpack the data =====================
        pos = 0;
        MPI_Unpack( &buffer[0], bytes, &pos, MPI_BOTTOM, 1, dtype_subarray_t, local_communicator );
        timing_record(4);
//! now for rank 1
      } else {
//! receive from rank 0      
        MPI_Recv( &buffer[0], bytes, MPI_PACKED, 0, itag, local_communicator, MPI_STATUS_IGNORE );
//! unpack the data 
        pos = 0;
        MPI_Unpack( &buffer[0], bytes, &pos, MPI_BOTTOM, 1, dtype_subarray_t, local_communicator );
//! pack the data
        pos = 0;
        MPI_Pack( MPI_BOTTOM, 1, dtype_subarray_t, &buffer[0], bytes, &pos, local_communicator );
//! send to rank 0
        MPI_Send( &buffer[0], pos, MPI_PACKED, 0, itag, local_communicator );
      } //! of myrank .EQ. 0?

    } //! inner_loop

//! ======================= clean up ========================

    MPI_Type_free( &dtype_subarray_t );

    free( buffer );

    if ( myrank == 0 ) {
      timing_record(5);
    }

  } //! outer_loop

  if ( myrank == 0 ) {
    timing_print( 1 );
  }

//! ======================= clean up ========================

  for( m=0 ; m<number_2D ; m++ ) {
    free( array2Ds[m] );
  }

  for( m=0 ; m<number_3D ; m++ ) {
    free( array3Ds[m] );
  }

  for( m=0 ; m<number_4D ; m++ ) {
    free(array4Ds[m] );
  }
     
  free( array2Ds );
  free( array3Ds );
  free( array4Ds );
}

#endif

