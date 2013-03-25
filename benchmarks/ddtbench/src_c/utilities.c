#include <stdlib.h>
#include <string.h>

//! gives list_dim unique numbers in index list back (from the range of
//! 1:global_dim

void utilities_random_array_shuffle( int* index_list, int list_dim, int global_dim ) {

  int i, temp, irandom;
  int* shuffle_array;

  shuffle_array = malloc( global_dim * sizeof(int) );

  for( i=0 ; i<global_dim ; i++ ) {
//! it is i+1 to be compliant with the fortran code
    shuffle_array[i] = i+1;
  }

  for( i=0 ; i<global_dim ; i++ ) {
    irandom = (rand() % global_dim);

    temp = shuffle_array[i];
    shuffle_array[i] = shuffle_array[irandom];
    shuffle_array[irandom] = temp;
  }

  memcpy( &index_list[0], &shuffle_array[0], list_dim*sizeof(int) );

  free( shuffle_array );

}

//! the following subroutines initialize the elements of an (1D/2D/3D/4D) array
//! with a unique number beginning from base
//! for each needed datatype there is a extra subroutine

void utilities_fill_unique_array_1D_float( float* array, int DIM1, int base ) {

  int i;

  for( i=1 ; i<DIM1 ; i++ ) {
    array[i] = base + i;
  }
}

void utilities_fill_unique_array_2D_float( float* array, int DIM1, int DIM2, int base ) {
      
  int i, j;

  for( j=0 ; j<DIM2 ; j++ ) {
    for( i=0 ; i<DIM1 ; i++ ) {
      array[i+j*DIM1] = base++;
    }
  }
}

void utilities_fill_unique_array_3D_float( float* array, int DIM1, int DIM2, int DIM3, int base ) {

  int i, j, k;

  for( k=0 ; k<DIM3 ; k++ ) {
    for( j=0 ; j<DIM2 ; j++ ) {
      for( i=0 ; i<DIM1 ; i++ ) {
        array[i+DIM1*(j+DIM2*k)] = base++;
      }
    }
  }
}

void utilities_fill_unique_array_4D_float( float* array, int DIM1, int DIM2, int DIM3, int DIM4, int base ) {

  int i, j, k, l;

  for( l=0 ; l<DIM4 ; l++ ) {
    for( k=0 ; k<DIM3 ; k++ ) {
      for( j=0 ; j<DIM2 ; j++ ) {
        for( i=0 ; i<DIM1 ; i++ ) {
          array[i+DIM1*(j+DIM2*(k+DIM3*l))] = base++;
        }
      }
    }
  }
}

void utilities_fill_unique_array_5D_float( float* array, int DIM1, int DIM2, int DIM3, int DIM4, int DIM5, int base ) { 

  int i, j, k, l, m;

  for( m=0 ; m<DIM5 ; m++ ) {
    for( l=0 ; l<DIM4 ; l++ ) {
      for( k=0 ; k<DIM3 ; k++ ) {
        for( j=0 ; j<DIM2 ; j++ ) {
          for( i=0 ; i<DIM1 ; i++ ) {
            array[i+DIM1*(j+DIM2*(k+DIM3*(l+DIM4*m)))] = base++;
          }
        }
      }
    }
  }
}

void utilities_fill_unique_array_1D_double( double* array, int DIM1, int base ) {

  int i;

  for( i=0 ; i<DIM1 ; i++) {
    array[i] = base++;
  }
}

void utilities_fill_unique_array_2D_double( double* array, int DIM1, int DIM2, int base ) {

  int i, j;

  for( j=0 ; j<DIM2 ; j++ ) {
    for( i=0 ; i<DIM1 ; i++ ) {
      array[i+j*DIM1] = base++;
    }
  }
}

void utilities_fill_unique_array_3D_double( double* array, int DIM1, int DIM2, int DIM3, int base ) {

  int i, j, k;

  for( k=0 ; k<DIM3 ; k++ ) {
    for( j=0 ; j<DIM2 ; j++ ) {
      for( i=0 ; i<DIM1 ; i++ ) {
        array[i+DIM1*(j+DIM2*k)] = base++;
      }
    }
  }
}
