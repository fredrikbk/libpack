#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mpi.h"

#include "ddtbench.h"

//! contains a wrapper for each test, so that the main program don't need
//! to call every method by himself
//! also this is the part, where the parameter checking should be put
//! also the output for the correctness checking should be placed here

void wrapper_timing_wrf_vec( int number_2D, int number_3D, int number_4D, int ims, int ime, int jms, int jme, int kms, int kme, int* limit_4D_arrays, int is, int ie, int js, int je, 
  int ks, int ke, int param_first_scalar, int outer_loop, int inner_loop, MPI_File filehandle_correctness, MPI_File filehandle_debug, char* ptestname, MPI_Comm local_communicator ) {

  int correct_flag;
  int typesize;
  char testname[50];
  int nelements, loops;
  int m;

//! just some statements to prevent compiler warnings of unused variables
//! those parameter are included for future features
//   typesize = filehandle_correctness

  if ( strlen(ptestname) == 0 ) {
    snprintf(&testname[0], 50, "WRF_vec" );
  } else {
    strncpy(&testname[0], &ptestname[0], 50 );
  }

  timing_wrf_vec_ddt( number_2D, number_3D, number_4D, ims, ime, jms, jme, kms, kme, limit_4D_arrays, is, ie, js, je, ks, ke, param_first_scalar, outer_loop, inner_loop, 
    &correct_flag, &typesize, &testname[0], filehandle_debug, local_communicator );

  timing_wrf_manual( number_2D, number_3D, number_4D, ims, ime, jms, jme, kms, kme, limit_4D_arrays, is, ie, js, je, ks, ke, param_first_scalar, outer_loop, inner_loop, 
    &correct_flag, &typesize, &testname[0], filehandle_debug, local_communicator );

  timing_wrf_vec_mpi_pack_ddt( number_2D, number_3D, number_4D, ims, ime, jms, jme, kms, kme, limit_4D_arrays, is, ie, js, je, ks, ke, param_first_scalar, outer_loop, inner_loop, 
    &correct_flag, &typesize, &testname[0], filehandle_debug, local_communicator );

  nelements = number_2D * (ie-is+1) * (je-js+1) + number_3D * (ie-is+1) * (je-js+1) * (ke-ks+1);
  for( m=0 ; m<number_4D ; m++ ) {
    if (limit_4D_arrays[m] >= param_first_scalar) {
      nelements = nelements + (limit_4D_arrays[m]-param_first_scalar+1) * (ie-is+1) * (je-js+1) * (ke-ks+1);
    }
  }
  loops = outer_loop * inner_loop;
  timing_basic_ping_pong_nelements( nelements, loops, &testname[0], local_communicator );

}

#if 0
void wrapper_timing_wrf_sa( int number_2D, int number_3D, int number_4D, int ims, int ime, int jms, int jme, int kms, int kme, int* limit_4D_arrays, int is, int ie, int js, int je, 
  int ks, int ke, int param_first_scalar, int outer_loop, int inner_loop, MPI_File filehandle_correctness, MPI_File filehandle_debug, char* ptestname, MPI_Comm local_communicator ) {
      
  int correct_flag;
  int typesize;
  char testname[50];
  int nelements, loops;
  int m;

//! just some statements to prevent compiler warnings of unused variables
//! those parameter are included for future features
//  typesize = filehandle_correctness

  if ( strlen(ptestname) == 0 ) {
    snprintf( &testname[0], 50, "WRF_sa" );
  } else {
    strncpy( &testname[0], &ptestname[0], 50 );
  }

  timing_wrf_sa_ddt( number_2D, number_3D, number_4D, ims, ime, jms, jme, kms, kme, limit_4D_arrays, is, ie, js, je, ks, ke, param_first_scalar, outer_loop, inner_loop, 
    &correct_flag, &typesize, &testname[0], filehandle_debug, local_communicator );

  timing_wrf_manual( number_2D, number_3D, number_4D, ims, ime, jms, jme, kms, kme, limit_4D_arrays, is, ie, js, je, ks, ke, param_first_scalar, outer_loop, inner_loop, 
    &correct_flag, &typesize, &testname[0], filehandle_debug, local_communicator );

  timing_wrf_sa_mpi_pack_ddt( number_2D, number_3D, number_4D, ims, ime, jms, jme, kms, kme, limit_4D_arrays, is, ie, js, je, ks, ke, param_first_scalar, outer_loop, inner_loop, 
    &correct_flag, &typesize, &testname[0], filehandle_debug, local_communicator );

  nelements = number_2D * (ie-is+1) * (je-js+1) + number_3D * (ie-is+1) * (je-js+1) * (ke-ks+1);
  for( m=0 ; m<number_4D ; m++ ) {
    if (limit_4D_arrays[m] >= param_first_scalar) {
      nelements = nelements + (limit_4D_arrays[m]-param_first_scalar+1) * (ie-is+1) * (je-js+1) * (ke-ks+1);
    }
  }
  loops = outer_loop * inner_loop;
  timing_basic_ping_pong_nelements( nelements, loops, &testname[0], local_communicator );
}
#endif

void wrapper_timing_wrf( int number_2D, int number_3D, int number_4D, int ims, int ime, int jms, int jme, int kms, int kme, int* limit_4D_arrays, int is, int ie, int js, int je, 
  int ks, int ke, int param_first_scalar, int outer_loop, int inner_loop, MPI_File filehandle_correctness, MPI_File filehandle_debug, char* ptestname, MPI_Comm local_communicator ){

  wrapper_timing_wrf_vec( number_2D, number_3D, number_4D, ims, ime, jms, jme, kms, kme, limit_4D_arrays, is, ie, js, je, ks, ke, param_first_scalar, outer_loop, inner_loop, 
    filehandle_correctness, filehandle_debug, &ptestname[0], local_communicator );
#if 0
  wrapper_timing_wrf_sa( number_2D, number_3D, number_4D, ims, ime, jms, jme, kms, kme, limit_4D_arrays, is, ie, js, je, ks, ke, param_first_scalar, outer_loop, inner_loop, 
    filehandle_correctness, filehandle_debug, &ptestname[50], local_communicator );
#endif
}

void wrapper_timing_milc_su3_zdown( int DIM2, int DIM3, int DIM4, int DIM5, int outer_loop, int inner_loop, MPI_File filehandle_correctness, MPI_File filehandle_debug, char* ptestname, MPI_Comm local_communicator ) {

  int correct_flag;
  int typesize;
  char testname[50];
  int nelements, loops;

//! just some statements to prevent compiler warnings of unused variables
//! those parameter are included for future features
//  typesize = filehandle_correctness

  if ( strlen(ptestname) == 0 ) {
    snprintf( &testname[0], 50, "MILC_su3_zd" );
  } else {
    strncpy( &testname[0], &ptestname[0], 50 );
  }

  timing_milc_su3_zdown_ddt( DIM2, DIM3, DIM4, DIM5, outer_loop, inner_loop, &correct_flag, &typesize, &testname[0], filehandle_debug, local_communicator );

  timing_milc_su3_zdown_manual( DIM2, DIM3, DIM4, DIM5, outer_loop, inner_loop, &correct_flag, &typesize, &testname[0], filehandle_debug, local_communicator );

  timing_milc_su3_zdown_mpi_pack_ddt( DIM2, DIM3, DIM4, DIM5, outer_loop, inner_loop, &correct_flag, &typesize, &testname[0], filehandle_debug, local_communicator );

//! not necessarily correct, since it assumes that a complex uses twice
//! the bytes a real does
   nelements = DIM2*DIM3/2*DIM5*2*3*2;
   loops = outer_loop * inner_loop;
   timing_basic_ping_pong_nelements( nelements, loops, &testname[0], local_communicator );
}

void wrapper_timing_nas_lu_x( int DIM2, int DIM3, int outer_loop, int inner_loop, MPI_File filehandle_correctness, MPI_File filehandle_debug, char* ptestname, MPI_Comm local_communicator ) {
      
  int correct_flag;
  int typesize;
  char testname[50];
  int nelements, loops;

//! just some statements to prevent compiler warnings of unused variables
//! those parameter are included for future features
//      typesize = filehandle_correctness

  if ( strlen(ptestname) == 0 ) {
    snprintf( testname, 50, "NAS_LU_x" );
  } else {
    strncpy( &testname[0], &ptestname[0], 50 );
  }

  timing_nas_lu_x_ddt( DIM2, DIM3, outer_loop, inner_loop, &correct_flag, &typesize, &testname[0], filehandle_debug, local_communicator );

  timing_nas_lu_x_manual( DIM2, DIM3, outer_loop, inner_loop, &correct_flag, &typesize, &testname[0], filehandle_debug, local_communicator );

  timing_nas_lu_x_mpi_pack_ddt( DIM2, DIM3, outer_loop, inner_loop, &correct_flag, &typesize, &testname[0], filehandle_debug, local_communicator );

//! not necessarily correct, since it assumes that a double uses twice
//! the bytes a real does
  nelements = 5 * DIM2 * 2;
  loops = outer_loop * inner_loop;
  timing_basic_ping_pong_nelements( nelements, loops, &testname[0], local_communicator );
}

void wrapper_timing_nas_lu_y( int DIM2, int DIM3, int outer_loop, int inner_loop, MPI_File filehandle_correctness, MPI_File filehandle_debug, char* ptestname, MPI_Comm local_communicator ) {

  int correct_flag;
  int typesize;
  char testname[50];
  int nelements, loops;

//! just some statements to prevent compiler warnings of unused variables
//! those parameter are included for future features
//      typesize = filehandle_correctness

  if ( strlen(ptestname) == 0 ) {
    snprintf(&testname[0],50,"NAS_LU_y");
  } else {
    strncpy(&testname[0], &ptestname[0], 50 );
  }

  timing_nas_lu_y_ddt( DIM2, DIM3, outer_loop, inner_loop, &correct_flag, &typesize, &testname[0], filehandle_debug, local_communicator );

  timing_nas_lu_y_manual( DIM2, DIM3, outer_loop, inner_loop, &correct_flag, &typesize, &testname[0], filehandle_debug, local_communicator );

  timing_nas_lu_y_mpi_pack_ddt( DIM2, DIM3, outer_loop, inner_loop, &correct_flag, &typesize, &testname[0], filehandle_debug, local_communicator );

//! not necessarily correct, since it assumes that a double uses twice
//! the bytes a real does
  nelements = 5 * DIM3 * 2;
  loops = outer_loop * inner_loop;
  timing_basic_ping_pong_nelements( nelements, loops, &testname[0], local_communicator );
}

void wrapper_timing_nas_lu( int DIM2, int DIM3, int outer_loop, int inner_loop, MPI_File filehandle_correctness, MPI_File filehandle_debug, char* ptestname, MPI_Comm local_communicator ) {

  wrapper_timing_nas_lu_x( DIM2, DIM3, outer_loop, inner_loop, filehandle_correctness, filehandle_debug, &ptestname[0], local_communicator );

  wrapper_timing_nas_lu_y( DIM2, DIM3, outer_loop, inner_loop, filehandle_correctness, filehandle_debug, &ptestname[50], local_communicator );

}

void wrapper_timing_nas_mg_x( int DIM1, int DIM2, int DIM3, int outer_loop, int inner_loop, MPI_File filehandle_correctness, MPI_File filehandle_debug, char* ptestname, MPI_Comm local_communicator ) {
      
  int correct_flag;
  int typesize;
  char testname[50];
  int nelements, loops;

//! just some statements to prevent compiler warnings of unused variables
//! those parameter are included for future features
//  typesize = filehandle_correctness

  if ( strlen(ptestname) == 0 ) {
    snprintf(&testname[0], 50, "NAS_MG_x");
  } else {
    strncpy( &testname[0], &ptestname[0], 50 );
  }

  timing_nas_mg_x_ddt( DIM1, DIM2, DIM3, outer_loop, inner_loop, &correct_flag, &typesize, &testname[0], filehandle_debug, local_communicator );

  timing_nas_mg_x_manual( DIM1, DIM2, DIM3, outer_loop, inner_loop, &correct_flag, &typesize, &testname[0], filehandle_debug, local_communicator );

  timing_nas_mg_x_mpi_pack_ddt( DIM1, DIM2, DIM3, outer_loop, inner_loop, &correct_flag, &typesize, &testname[0], filehandle_debug, local_communicator );

//! not necessarily correct, since it assumes that a double uses twice
//! the bytes a real does
  nelements = (DIM2-2)*(DIM3-2)*2;
  loops = outer_loop * inner_loop;
  timing_basic_ping_pong_nelements( nelements, loops, &testname[0], local_communicator );

}

void wrapper_timing_nas_mg_y( int DIM1, int DIM2, int DIM3, int outer_loop, int inner_loop, MPI_File filehandle_correctness, MPI_File filehandle_debug, char* ptestname, MPI_Comm local_communicator ) {
      
  int correct_flag;
  int typesize;
  char testname[50];
  int nelements, loops;

//! just some statements to prevent compiler warnings of unused variables
//! those parameter are included for future features
//  typesize = filehandle_correctness

  if ( strlen(ptestname) == 0 ) {
    snprintf(&testname[0], 50, "NAS_MG_y");
   } else {
    strncpy( &testname[0], &ptestname[0], 50 );
   }

   timing_nas_mg_y_ddt( DIM1, DIM2, DIM3, outer_loop, inner_loop, &correct_flag, &typesize, &testname[0], filehandle_debug, local_communicator );

   timing_nas_mg_y_manual( DIM1, DIM2, DIM3, outer_loop, inner_loop, &correct_flag, &typesize, &testname[0], filehandle_debug, local_communicator );

   timing_nas_mg_y_mpi_pack_ddt( DIM1, DIM2, DIM3, outer_loop, inner_loop, &correct_flag, &typesize, &testname[0], filehandle_debug, local_communicator );

//! not necessarily correct, since it assumes that a double uses twice
//! the bytes a real does
    nelements = (DIM1-2) * (DIM3-2) * 2;
    loops = outer_loop * inner_loop;
    timing_basic_ping_pong_nelements( nelements, loops, &testname[0], local_communicator );

}

void wrapper_timing_nas_mg_z( int DIM1, int DIM2, int DIM3, int outer_loop, int inner_loop, MPI_File filehandle_correctness, MPI_File filehandle_debug, char* ptestname, MPI_Comm local_communicator ) {

 int correct_flag;
 int typesize;
 char testname[50];
 int nelements, loops;

//! just some statements to prevent compiler warnings of unused variables
//! those parameter are included for future features
//  typesize = filehandle_correctness

  if ( strlen(ptestname) == 0 ) {
    snprintf( &testname[0], 50, "NAS_MG_z" );
  } else {
    strncpy( &testname[0], &ptestname[0], 50 );
  }

  timing_nas_mg_z_ddt( DIM1, DIM2, DIM3, outer_loop, inner_loop, &correct_flag, &typesize, &testname[0], filehandle_debug, local_communicator );

  timing_nas_mg_z_manual( DIM1, DIM2, DIM3, outer_loop, inner_loop, &correct_flag, &typesize, &testname[0], filehandle_debug, local_communicator );

  timing_nas_mg_z_mpi_pack_ddt( DIM1, DIM2, DIM3, outer_loop, inner_loop, &correct_flag, &typesize, &testname[0], filehandle_debug, local_communicator );

//! not necessarily correct, since it assumes that a double uses twice
//! the bytes a real does
  
  nelements = (DIM1-2) * (DIM2-2) * 2;
  loops = outer_loop * inner_loop;
  timing_basic_ping_pong_nelements( nelements, loops, &testname[0], local_communicator );

}

void wrapper_timing_nas_mg( int DIM1, int DIM2, int DIM3, int outer_loop, int inner_loop, MPI_File filehandle_correctness, MPI_File filehandle_debug, char* ptestname, MPI_Comm local_communicator ) {

  wrapper_timing_nas_mg_x( DIM1, DIM2, DIM3, outer_loop, inner_loop, filehandle_correctness, filehandle_debug, &ptestname[0], local_communicator );

  wrapper_timing_nas_mg_y( DIM1, DIM2, DIM3, outer_loop, inner_loop, filehandle_correctness, filehandle_debug, &ptestname[50], local_communicator );

  wrapper_timing_nas_mg_z( DIM1, DIM2, DIM3, outer_loop, inner_loop, filehandle_correctness, filehandle_debug, &ptestname[100], local_communicator );

}

void wrapper_timing_fft( int DIM1, int outer_loop, int inner_loop, MPI_File filehandle_correctness, MPI_File filehandle_debug, char* ptestname, MPI_Comm local_communicator ) {

  int correct_flag;
  int typesize;
  char testname[50];
  int nelements, loops;
  int myrank, procs;

//! just some statements to prevent compiler warnings of unused variables
//! those parameter are included for future features
//  typesize = filehandle_correctness

  MPI_Comm_rank( local_communicator, &myrank );
  MPI_Comm_size( local_communicator, &procs );

  DIM1 = DIM1 + DIM1 % procs;

  if ( strlen(ptestname) == 0 ) {
    snprintf( &testname[0], 50, "FFT%i", procs);
  } else {
    strncpy( &testname[0], &ptestname[0], 50 );
  }

  timing_fft2d_ddt( DIM1, procs, outer_loop, inner_loop, &correct_flag, &typesize, &testname[0], filehandle_debug, local_communicator );

  timing_fft2d_manual( DIM1, procs, outer_loop, inner_loop, &correct_flag, &typesize, &testname[0], filehandle_debug, local_communicator );

  timing_fft2d_mpi_pack_ddt( DIM1, procs, outer_loop, inner_loop, &correct_flag, &typesize, &testname[0], filehandle_debug, local_communicator );

//! not necessarily correct, since it assumes that a double uses twice
//! the bytes a real does
  nelements = DIM1/procs*DIM1/procs * 2 * 2;
  loops = outer_loop * inner_loop;
  timing_basic_alltoall_nelements( nelements, procs, loops, &testname[0], local_communicator );
}
      
void wrapper_timing_specfem3d_mt( int DIM1, int DIM2, int DIM3, int outer_loop, int inner_loop, MPI_File filehandle_correctness, MPI_File filehandle_debug, char* ptestname, MPI_Comm local_communicator ) {

  int correct_flag;
  int typesize;
  char testname[50];
  int nelements, loops;

//! just some statements to prevent compiler warnings of unused variables
//! those parameter are included for future features
//  typesize = filehandle_correctness;

  if ( strlen(ptestname) == 0 ) {
    snprintf( &testname[0], 50, "SPECFEM3D_mt" );
  } else {
    strncpy( &testname[0], &ptestname[0], 50 );
  }

  timing_specfem3d_mt_ddt( DIM1, DIM2, DIM3, outer_loop, inner_loop, &correct_flag, &typesize, &testname[0], filehandle_debug, local_communicator );

  timing_specfem3d_mt_manual( DIM1, DIM2, DIM3, outer_loop, inner_loop, &correct_flag, &typesize, &testname[0], filehandle_debug, local_communicator );

  timing_specfem3d_mt_mpi_pack_ddt( DIM1, DIM2, DIM3, outer_loop, inner_loop, &correct_flag, &typesize, &testname[0], filehandle_debug, local_communicator );

  nelements = DIM1 * DIM3;
  loops = outer_loop * inner_loop;

  timing_basic_ping_pong_nelements( nelements, loops, &testname[0], local_communicator );

}

void wrapper_timing_lammps_full( int DIM1, int icount, int outer_loop, int inner_loop, MPI_File filehandle_correctness, MPI_File filehandle_debug, char* ptestname, MPI_Comm local_communicator ) {

  int correct_flag;
  int typesize;
  char testname[50];
  int nelements, loops, i;
  int* list;

//! just some statements to prevent compiler warnings of unused variables
//! those parameter are included for future features
//  typesize = filehandle_correctness
  
  list = malloc( icount * outer_loop * sizeof(int) );

  if ( strlen(ptestname) == 0 ) {
    snprintf( &testname[0], 50, "LAMMPS_full" );
  } else {
    strncpy( &testname[0], &ptestname[0], 50 );
  }

  for( i=0 ; i<outer_loop ; i++ ) {
    utilities_random_array_shuffle( &list[i*icount], icount, DIM1 );
  }

  timing_lammps_full_ddt( DIM1, icount, list, outer_loop, inner_loop, &correct_flag, &typesize, &testname[0], filehandle_debug, local_communicator );

  timing_lammps_full_manual( DIM1, icount, list, outer_loop, inner_loop, &correct_flag, &typesize, &testname[0], filehandle_debug, local_communicator );

  timing_lammps_full_mpi_pack_ddt( DIM1, icount, list, outer_loop, inner_loop, &correct_flag, &typesize, &testname[0], filehandle_debug, local_communicator );

//! not necessarily correct, since it assumes that a double uses twice
//! the bytes a real does
  nelements = icount * 8 * 2;
  loops = outer_loop * inner_loop;
  timing_basic_ping_pong_nelements( nelements, loops, testname, local_communicator );

  free(list);
}

void wrapper_timing_lammps_atomic( int DIM1, int icount, int outer_loop, int inner_loop, MPI_File filehandle_correctness, MPI_File filehandle_debug, char* ptestname, MPI_Comm local_communicator ) {

  int correct_flag;
  int typesize;
  char testname[50];
  int nelements, loops, i;
  int* list;

//! just some statements to prevent compiler warnings of unused variables
//! those parameter are included for future features
//  typesize = filehandle_correctness

  list = malloc( icount * outer_loop * sizeof(int) );

  if ( strlen(ptestname) == 0 ) {
    snprintf( &testname[0], 50, "LAMMPS_atomic" );
  } else {
    strncpy( &testname[0], &ptestname[0], 50 );
  }

  for( i=0 ; i<outer_loop ; i++ ) {
    utilities_random_array_shuffle( &list[i*icount], icount, DIM1 );
  }

  timing_lammps_atomic_ddt( DIM1, icount, list, outer_loop, inner_loop, &correct_flag, &typesize, &testname[0], filehandle_debug, local_communicator );

  timing_lammps_atomic_manual( DIM1, icount, list, outer_loop, inner_loop, &correct_flag, &typesize, &testname[0], filehandle_debug, local_communicator );

  timing_lammps_atomic_mpi_pack_ddt( DIM1, icount, list, outer_loop, inner_loop, &correct_flag, &typesize, &testname[0], filehandle_debug, local_communicator );

//! not necessarily correct, since it assumes that a double uses twice
//! the bytes a real does
  nelements = icount * 6 * 2;
  loops = outer_loop * inner_loop;
  timing_basic_ping_pong_nelements( nelements, loops, &testname[0], local_communicator );

  free(list);
}

void wrapper_timing_specfem3D_oc( int DIM1, int icount, int outer_loop, int inner_loop, MPI_File filehandle_correctness, MPI_File filehandle_debug, char* ptestname, MPI_Comm local_communicator ) {

  int correct_flag;
  int typesize;
  char testname[50];
  int nelements, loops, i;
  int* list;

  list = malloc( icount * outer_loop * sizeof(int) );

//! just some statements to prevent compiler warnings of unused variables
//! those parameter are included for future features
//  typesize = filehandle_correctness

  if ( strlen(ptestname) == 0 ) {
    snprintf( &testname[0], 50, "SPECFEM3D_oc" );
  } else {
    strncpy( &testname[0], &ptestname[0], 50 );
  }

  for( i=0 ; i<outer_loop ; i++ ) {
    utilities_random_array_shuffle( &list[i*icount], icount, DIM1 );
  }

  timing_specfem3D_oc_ddt( DIM1, icount, &list[0], outer_loop, inner_loop, &correct_flag, &typesize, &testname[0], filehandle_debug, local_communicator );

  timing_specfem3D_oc_manual( DIM1, icount, &list[0], outer_loop, inner_loop, &correct_flag, &typesize, &testname[0], filehandle_debug, local_communicator );

  timing_specfem3D_oc_mpi_pack_ddt( DIM1, icount, &list[0], outer_loop, inner_loop, &correct_flag, &typesize, &testname[0], filehandle_debug, local_communicator );

  nelements = icount;
  loops = outer_loop * inner_loop;
  timing_basic_ping_pong_nelements( nelements, loops, &testname[0], local_communicator );

  free(list);
}

void wrapper_timing_specfem3D_cm( int DIM2_cm, int DIM2_ic, int icount_cm, int icount_ic, int outer_loop, int inner_loop, MPI_File filehandle_correctness, MPI_File filehandle_debug, char* ptestname, MPI_Comm local_communicator ) {

  int correct_flag;
  int typesize;
  char testname[50];
  int nelements, loops, i;

  int* list_cm;
  int* list_ic;

  list_cm = malloc( icount_cm * outer_loop * sizeof(int) );
  list_ic = malloc( icount_ic * outer_loop * sizeof(int) );

//! just some statements to prevent compiler warnings of unused variables
//! those parameter are included for future features
//   typesize = filehandle_correctness;

  if ( strlen(ptestname) == 0 ) {
    snprintf( &testname[0], 50, "SPECFEM3D_cm" );
  } else {
    strncpy( &testname[0], &ptestname[0], 50 );
  }

  for( i=0 ; i<outer_loop ; i++ ) {
    utilities_random_array_shuffle( &list_cm[i*icount_cm], icount_cm, DIM2_cm );
    utilities_random_array_shuffle( &list_ic[i*icount_ic], icount_ic, DIM2_ic );
  }

  timing_specfem3D_cm_ddt( DIM2_cm, DIM2_ic, icount_cm, icount_ic, &list_cm[0], &list_ic[0], outer_loop, inner_loop, &correct_flag, &typesize, &testname[0], filehandle_debug, local_communicator );

  timing_specfem3D_cm_manual( DIM2_cm, DIM2_ic, icount_cm, icount_ic, &list_cm[0], &list_ic[0], outer_loop, inner_loop, &correct_flag, &typesize, &testname[0], filehandle_debug, local_communicator );

  timing_specfem3D_cm_mpi_pack_ddt( DIM2_cm, DIM2_ic, icount_cm, icount_ic, &list_cm[0], &list_ic[0], outer_loop, inner_loop, &correct_flag, &typesize, &testname[0], filehandle_debug, local_communicator );

  nelements = (icount_cm+icount_ic) * 3;
  loops = outer_loop * inner_loop;
  timing_basic_ping_pong_nelements( nelements, loops, &testname[0], local_communicator );

  free(list_cm);
  free(list_ic);
}
