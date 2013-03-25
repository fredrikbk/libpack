#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "mpi.h"
#include "ddtbench.h"

#include "../hrtimer/hrtimer.h"


int main( int argc, char **argv, char *envp[]) {

  int outer_loop = 10;
  int inner_loop = 20;

  int myrank;
  int max_epochs;

  MPI_Comm local_comm_pp, local_comm_all2all;
      
//! variables for file handling
  char filename[50];
  char testname[3][50];   
  MPI_File filehandle_correctness, filehandle_debug;

  int DIM1, DIM2, DIM3, DIM4;

  int number_2D, number_3D, number_4D;
  int ims, ime, jms, jme, kms, kme;
  int is, ie, js, je, ks, ke;
  int param_first_scalar;
  int *limit_4D_arrays;

  int icount1, icount2;

  int i;

//! ====================================================================
//! =============================== intro ==============================
//! ====================================================================


//! init some variables for file handling
  testname[0][0] = '\0';
  testname[1][0] = '\0';
  testname[2][0] = '\0';
  filehandle_correctness = MPI_FILE_NULL;
  filehandle_debug = MPI_FILE_NULL;

  MPI_Init( &argc, &argv );

#if TEST_TYPE != 2
  timing_hrt_init();
#endif
 
  MPI_Comm_rank( MPI_COMM_WORLD, &myrank );

  MPI_Comm_dup( MPI_COMM_WORLD, &local_comm_all2all );

  if ( myrank < 2 ) {
    MPI_Comm_split( MPI_COMM_WORLD, 0, myrank, &local_comm_pp );
  } else {
    MPI_Comm_split( MPI_COMM_WORLD, MPI_UNDEFINED, myrank, &local_comm_pp );
  }

//  init the random number generator
  unsigned int iseed = (unsigned int)time(NULL);
  srand (iseed);

//! reading the arguments
  if ( argc != 1 ) {
    if ( argc != 3 ) {
      if ( myrank  == 0 ) {
        printf("usage ddtbench <outer_loop> <inner_loop>\n");
        printf("proceed with standard parameter for the inner and outer loop\n\n");
      }
    } else {
      outer_loop = atoi(argv[1]);
      inner_loop = atoi(argv[2]);
    }
  } 

//! some intro
  if ( myrank == 0 ) {
    printf(" Welcome to our DDT benchmark suite\n\n");
    printf(" outer loop parameter = %i\n", outer_loop);
    printf(" inner loop parameter = %i\n\n", inner_loop);
//! open the filehandle for the output of the timing
    snprintf(filename,50,"ddtbench.out");
    timing_open_file( &filename[0] );
#if TEST_TYPE > 1
    init_papi();
#endif
//! set the maximum test number for the progression bar
//! ddt has 1 epoch in inner loop, 2 epochs in outer loop, and 2 epochs 
//! for each test
//! manual has 3 epochs in inner loop, 2 epochs in outer loop
//! mpi_pack_ddt has 3 epochs in inner loop, 2 epochs in outer looper
//! reference has 1 epoch in loops (=outer_loops*inner_loops) and 4
//! epochs for each test
//! there are 65 tests in sum
    max_epochs = 65 * (outer_loop * (inner_loop * (1+3+3+1) + (2+2+2+0)) + (2+0+0+4)) + 4 * (outer_loop * (inner_loop * (1+2+3+1) + (2+2+2+0)) + (2+1+0+4));
    timing_set_max_tests( max_epochs );
  }

  if ( myrank < 2 ) {

//! ====================================================================
//! ========================== WRF y direction =========================
//! ====================================================================

    number_2D = 4;
    number_3D = 3;
    number_4D = 2;

    limit_4D_arrays = malloc( number_4D * sizeof(int));

    snprintf(&testname[0][0],50,"WRF_y_vec");
    snprintf(&testname[1][0],50,"WRF_y_sa");

//! 2x2 ym send, em_b_wave case
    ims = -4;
    ime = 27;
    kms = 1;
    kme = 65;
    jms = 34;
    jme = 85;
    for (i = 0 ; i < number_4D ; i++) {
      limit_4D_arrays[i] = 2;
    }

    is = 1;
    ie = 23;
    ks = 1;
    ke = 65;
    js = 41;
    je = 43;
    param_first_scalar = 2;

    wrapper_timing_wrf( number_2D, number_3D, number_4D, ims, ime, jms, jme, kms, kme, limit_4D_arrays, is, ie, js, je, ks, ke, param_first_scalar, outer_loop, 
      inner_loop, filehandle_correctness, filehandle_debug, &testname[0][0], local_comm_pp );

//! 3x3 ym send, em_b_wave case
    ims = 8;
    ime = 34;
    kms = 1;
    kme = 65;
    jms = 21;
    jme = 61;
    for (i = 0 ; i < number_4D ; i++) {
      limit_4D_arrays[i] = 2;
    }

    is = 12;
    ie = 30;
    ks = 1;
    ke = 65;
    js = 28;
    je = 30;
    param_first_scalar = 2;

    wrapper_timing_wrf( number_2D, number_3D, number_4D, ims, ime, jms, jme, kms, kme, limit_4D_arrays, is, ie, js, je, ks, ke, param_first_scalar, outer_loop, 
      inner_loop, filehandle_correctness, filehandle_debug, &testname[0][0], local_comm_pp );

//! 4x4 ym send, em_b_wave case
    ims = 4;
    ime = 27;
    kms = 1;
    kme = 65;
    jms = 14;
    jme = 47;
    for (i = 0 ; i < number_4D ; i++) {
      limit_4D_arrays[i] = 2;
    }

    is = 8;
    ie = 23;
    ks = 1;
    ke = 65;
    js = 21;
    je = 23;
    param_first_scalar = 2;

    wrapper_timing_wrf( number_2D, number_3D, number_4D, ims, ime, jms, jme, kms, kme, limit_4D_arrays, is, ie, js, je, ks, ke, param_first_scalar, outer_loop, 
      inner_loop, filehandle_correctness, filehandle_debug, &testname[0][0], local_comm_pp );

//! 5x5 ym send, em_b_wave case
    ims = 10;
    ime = 31;
    kms = 1;
    kme = 65;
    jms = 10;
    jme = 39;
    for (i = 0 ; i < number_4D ; i++) {
      limit_4D_arrays[i] = 2;
    }

    is = 14;
    ie = 27;
    ks = 1;
    ke = 65;
    js = 17;
    je = 19;
    param_first_scalar = 2;

    wrapper_timing_wrf( number_2D, number_3D, number_4D, ims, ime, jms, jme, kms, kme, limit_4D_arrays, is, ie, js, je, ks, ke, param_first_scalar, outer_loop, 
      inner_loop, filehandle_correctness, filehandle_debug, &testname[0][0], local_comm_pp );

//! 8x8 ym send, em_b_wave case
    ims = 4;
    ime = 22;
    kms = 1;
    kme = 65;
    jms = 4;
    jme = 27;
    for (i = 0 ; i < number_4D ; i++) {
      limit_4D_arrays[i] = 2;
    }
    
    is = 8;
    ie = 18;
    ks = 1;
    ke = 65;
    js = 11;
    je = 13;
    param_first_scalar = 2;

    wrapper_timing_wrf( number_2D, number_3D, number_4D, ims, ime, jms, jme, kms, kme, limit_4D_arrays, is, ie, js, je, ks, ke, param_first_scalar, outer_loop, 
      inner_loop, filehandle_correctness, filehandle_debug, &testname[0][0], local_comm_pp );

//! ====================================================================
//! ========================== WRF x direction =========================
//! ====================================================================

    snprintf(&testname[0][0],50,"WRF_x_vec");
    snprintf(&testname[1][0],50,"WRF_x_sa");

//! 2x2 xp send, em_b_wave
    ims = -4;
    ime = 27;
    kms = 1;
    kme = 65;
    jms = 34;
    jme = 85;
    for (i = 0 ; i < number_4D ; i++) {
      limit_4D_arrays[i] = 2;
    }

    is = 18;
    ie = 20;
    ks = 1;
    ke = 65;
    js = 38;
    je = 81;
    param_first_scalar = 2;

    wrapper_timing_wrf( number_2D, number_3D, number_4D, ims, ime, jms, jme, kms, kme, limit_4D_arrays, is, ie, js, je, ks, ke, param_first_scalar, outer_loop, 
      inner_loop, filehandle_correctness, filehandle_debug, &testname[0][0], local_comm_pp );

//! 4x4 xp send, em_b_wave
    ims = 4;
    ime = 27;
    kms = 1;
    kme = 65;
    jms = 14;
    jme = 47;
    for (i = 0 ; i < number_4D ; i++) {
      limit_4D_arrays[i] = 2;
    }

    is = 18;
    ie = 20;
    ks = 1;
    ke = 65;
    js = 18;
    je = 43;
    param_first_scalar = 2;

    wrapper_timing_wrf( number_2D, number_3D, number_4D, ims, ime, jms, jme, kms, kme, limit_4D_arrays, is, ie, js, je, ks, ke, param_first_scalar, outer_loop, 
      inner_loop, filehandle_correctness, filehandle_debug, &testname[0][0], local_comm_pp );

//! 8x8 xp send, em_b_wave
    ims = 4;
    ime = 22;
    kms = 1;
    kme = 65;
    jms = 4;
    jme = 27;
    for (i = 0 ; i < number_4D ; i++) {
      limit_4D_arrays[i] = 2;
    }

    is = 13;
    ie = 15;
    ks = 1;
    ke = 65;
    js = 8;
    je = 23;
    param_first_scalar = 2;

    wrapper_timing_wrf( number_2D, number_3D, number_4D, ims, ime, jms, jme, kms, kme, limit_4D_arrays, is, ie, js, je, ks, ke, param_first_scalar, outer_loop, 
      inner_loop, filehandle_correctness, filehandle_debug, &testname[0][0], local_comm_pp );

    free( limit_4D_arrays );

    testname[0][0] = '\0';
    testname[1][0] = '\0';

//! ====================================================================
//! ==================== MILC su3 zdown direction ======================
//! ====================================================================

    DIM1 = 16;
    DIM2 = 16;
    DIM3 = 16;
    DIM4 = 16;

    wrapper_timing_milc_su3_zdown( DIM1, DIM2, DIM3, DIM4, outer_loop, inner_loop, filehandle_correctness, filehandle_debug, &testname[0][0], local_comm_pp );
 
    DIM1 = 8;
    DIM2 = 8;
    DIM3 = 16;
    DIM4 = 16;

    wrapper_timing_milc_su3_zdown( DIM1, DIM2, DIM3, DIM4, outer_loop, inner_loop, filehandle_correctness, filehandle_debug, &testname[0][0], local_comm_pp );
  
    DIM1 = 8;
    DIM2 = 8;
    DIM3 = 8;
    DIM4 = 16;

    wrapper_timing_milc_su3_zdown( DIM1, DIM2, DIM3, DIM4, outer_loop, inner_loop, filehandle_correctness, filehandle_debug, &testname[0][0], local_comm_pp );

    DIM1 = 8;
    DIM2 = 8;
    DIM3 = 8;
    DIM4 = 8;

    wrapper_timing_milc_su3_zdown( DIM1, DIM2, DIM3, DIM4, outer_loop, inner_loop, filehandle_correctness, filehandle_debug, &testname[0][0], local_comm_pp );

    printf("MILC finished without error\n");

//! ====================================================================
//! =============================== NAS_LU =============================
//! ====================================================================

    DIM2 = 12;
    DIM3 = 12;

    wrapper_timing_nas_lu( DIM2, DIM3, outer_loop, inner_loop, filehandle_correctness, filehandle_debug, &testname[0][0], MPI_COMM_WORLD );

    DIM2 = 33;
    DIM3 = 33;

    wrapper_timing_nas_lu( DIM2, DIM3, outer_loop, inner_loop, filehandle_correctness, filehandle_debug, &testname[0][0], MPI_COMM_WORLD );

    DIM2 = 64;
    DIM3 = 64;

    wrapper_timing_nas_lu( DIM2, DIM3, outer_loop, inner_loop, filehandle_correctness, filehandle_debug, &testname[0][0], MPI_COMM_WORLD );

    DIM2 = 102;
    DIM3 = 102;

    wrapper_timing_nas_lu( DIM2, DIM3, outer_loop, inner_loop, filehandle_correctness, filehandle_debug, &testname[0][0], MPI_COMM_WORLD );

    DIM2 = 162;
    DIM3 = 162;

    wrapper_timing_nas_lu( DIM2, DIM3, outer_loop, inner_loop, filehandle_correctness, filehandle_debug, &testname[0][0], MPI_COMM_WORLD );

    DIM2 = 408;
    DIM3 = 408;

    wrapper_timing_nas_lu( DIM2, DIM3, outer_loop, inner_loop, filehandle_correctness, filehandle_debug, &testname[0][0], MPI_COMM_WORLD );

    DIM2 = 1020;
    DIM3 = 1020;

    wrapper_timing_nas_lu( DIM2, DIM3, outer_loop, inner_loop, filehandle_correctness, filehandle_debug, &testname[0][0], MPI_COMM_WORLD );


    printf("NAS LU finished without error\n");

//! ====================================================================
//! =============================== NAS_MG =============================
//! ====================================================================

    DIM1 = 34;
    DIM2 = 18;
    DIM3 = 18;

    wrapper_timing_nas_mg( DIM1, DIM2, DIM3, outer_loop, inner_loop, filehandle_correctness, filehandle_debug, &testname[0][0], local_comm_pp );

    DIM1 = 130;
    DIM2 = 66;
    DIM3 = 66;

    wrapper_timing_nas_mg( DIM1, DIM2, DIM3, outer_loop, inner_loop, filehandle_correctness, filehandle_debug, &testname[0][0], local_comm_pp );

    DIM1 = 258;
    DIM2 = 130;
    DIM3 = 130;

    wrapper_timing_nas_mg( DIM1, DIM2, DIM3, outer_loop, inner_loop, filehandle_correctness, filehandle_debug, &testname[0][0], local_comm_pp );

    DIM1 = 514;
    DIM2 = 258;
    DIM3 = 258;

    wrapper_timing_nas_mg( DIM1, DIM2, DIM3, outer_loop, inner_loop, filehandle_correctness, filehandle_debug, &testname[0][0], local_comm_pp );

    printf("NAS MG finished without error\n");

//! ====================================================================
//! ============================ LAMMPS_full ===========================
//! ====================================================================

// It is too huge...
#if 0

//! peptide example with 2 process (maximum)
    icount1 = 3062;
    DIM1 = 3534;

    wrapper_timing_lammps_full( DIM1, icount1, outer_loop, inner_loop, filehandle_correctness, filehandle_debug, &testname[0][0], local_comm_pp );

    printf("a\n");

//! peptide example with 4 process (maximum)
    icount1 = 2243;
    DIM1 = 2597;

    wrapper_timing_lammps_full( DIM1, icount1, outer_loop, inner_loop, filehandle_correctness, filehandle_debug, &testname[0][0], local_comm_pp );

    printf("b\n");
//! peptide example with 8 process (maximum)
    icount1 = 1662;
    DIM1 = 1907;

    wrapper_timing_lammps_full( DIM1, icount1, outer_loop, inner_loop, filehandle_correctness, filehandle_debug, &testname[0][0], local_comm_pp );

    printf("c\n");

    printf("LAMMPS_full finished without error\n");
#endif 

//! ====================================================================
//! =========================== LAMMPS_atomic ==========================
//! ====================================================================

//! crack example with 2 process (maximum)
    icount1 = 243;
    DIM1 = 4084;

    wrapper_timing_lammps_atomic( DIM1, icount1, outer_loop, inner_loop, filehandle_correctness, filehandle_debug, &testname[0][0], local_comm_pp );

//! crack example with 4 process (maximum)
    icount1 = 145;
    DIM1 = 2157;

    wrapper_timing_lammps_atomic( DIM1, icount1, outer_loop, inner_loop, filehandle_correctness, filehandle_debug, &testname[0][0], local_comm_pp );

//! crack example with 8 process (maximum)
    icount1 = 114;
    DIM1 = 1217;

    wrapper_timing_lammps_atomic( DIM1, icount1, outer_loop, inner_loop, filehandle_correctness, filehandle_debug, &testname[0][0], local_comm_pp );

    printf("LAMMPS_atomic finished without error\n");

//! ====================================================================
//! ============================ SPECFEM3D_oc ==========================
//! ====================================================================

//! 10x10x6, c=4
    icount1 = 3225;
    DIM1 = 88881;

    wrapper_timing_specfem3D_oc( DIM1, icount1, outer_loop, inner_loop, filehandle_correctness, filehandle_debug, &testname[0][0], local_comm_pp );

//! 10x10x6, c=3
    icount1 = 1897;
    DIM1 = 38585;

    wrapper_timing_specfem3D_oc( DIM1, icount1, outer_loop, inner_loop, filehandle_correctness, filehandle_debug, &testname[0][0], local_comm_pp );

//! 10x10x6, c=2
    icount1 = 877;
    DIM1 = 12857;

    wrapper_timing_specfem3D_oc( DIM1, icount1, outer_loop, inner_loop, filehandle_correctness, filehandle_debug, &testname[0][0], local_comm_pp );

//! 10x10x6, c=1
    icount1 = 493;
    DIM1 = 3697;

    wrapper_timing_specfem3D_oc( DIM1, icount1, outer_loop, inner_loop, filehandle_correctness, filehandle_debug, &testname[0][0], local_comm_pp );

    printf("SPECFEM3D_oc finished without error\n");

//! ====================================================================
//! ============================ SPECFEM3D_cm ==========================
//! ====================================================================
// too huge
#if 0

//! 10x10x6, c=4
    icount1 = 11797;
    icount2 = 3009;
    DIM1 = 834917;
    DIM2 = 51153;

    wrapper_timing_specfem3D_cm( DIM1, DIM2, icount1, icount2, outer_loop, inner_loop, filehandle_correctness, filehandle_debug, &testname[0][0], local_comm_pp );

//! 10x10x6, c=3
    icount1 = 7125;
    icount2 = 1729;
    DIM1 = 396849;
    DIM2 = 22477;

    wrapper_timing_specfem3D_cm( DIM1, DIM2, icount1, icount2, outer_loop, inner_loop, filehandle_correctness, filehandle_debug, &testname[0][0], local_comm_pp );

//! 10x10x6, c=2
    icount1 = 3801;
    icount2 = 801;
    DIM1 = 152001;
    DIM2 = 7209;

    wrapper_timing_specfem3D_cm( DIM1, DIM2, icount1, icount2, outer_loop, inner_loop, filehandle_correctness, filehandle_debug, &testname[0][0], local_comm_pp );

//! 10x10x6, c=1
    icount1 = 1957;
    icount2 = 245;
    DIM1 = 39929;
    DIM2 = 1225;

    wrapper_timing_specfem3D_cm( DIM1, DIM2, icount1, icount2, outer_loop, inner_loop, filehandle_correctness, filehandle_debug, &testname[0][0], local_comm_pp );


    printf("SPECFEM3D_cm finished without error\n");
#endif
//! ====================================================================
//! ============================ SPECFEM3D_mt ==========================
//! ====================================================================

//! 10x10x6, c=1
    DIM1 = 3;
    DIM2 = 2;
    DIM3 = 7600;

    wrapper_timing_specfem3d_mt( DIM1, DIM2, DIM3, outer_loop, inner_loop, filehandle_correctness, filehandle_debug, &testname[0][0], local_comm_pp );

//! 10x10x6, c=2
    DIM1 = 3;
    DIM2 = 2;
    DIM3 = 6400;

    wrapper_timing_specfem3d_mt( DIM1, DIM2, DIM3, outer_loop, inner_loop, filehandle_correctness, filehandle_debug, &testname[0][0], local_comm_pp );

//! 10x10x6, c=3
    DIM1 = 3;
    DIM2 = 2;
    DIM3 = 5600;

    wrapper_timing_specfem3d_mt( DIM1, DIM2, DIM3, outer_loop, inner_loop, filehandle_correctness, filehandle_debug, &testname[0][0], local_comm_pp );

//! 10x10x6, c=4
    DIM1 = 3;
    DIM2 = 2;
    DIM3 = 6200;

    wrapper_timing_specfem3d_mt( DIM1, DIM2, DIM3, outer_loop, inner_loop, filehandle_correctness, filehandle_debug, &testname[0][0], local_comm_pp );


  } //! of myrank < 2

//! ====================================================================
//! ================================ FFT ===============================
//! ====================================================================

  MPI_Barrier( local_comm_all2all );

#if 0

  DIM1 = 256;
  wrapper_timing_fft( DIM1, outer_loop, inner_loop, filehandle_correctness, filehandle_debug, &testname[0][0], local_comm_all2all );

  DIM1 = 512;
  wrapper_timing_fft( DIM1, outer_loop, inner_loop, filehandle_correctness, filehandle_debug, &testname[0][0], local_comm_all2all );

  DIM1 = 768;
  wrapper_timing_fft( DIM1, outer_loop, inner_loop, filehandle_correctness, filehandle_debug, &testname[0][0], local_comm_all2all );

  DIM1 = 1024;
  wrapper_timing_fft( DIM1, outer_loop, inner_loop, filehandle_correctness, filehandle_debug, &testname[0][0], local_comm_all2all );

  DIM1 = 1536;
  wrapper_timing_fft( DIM1, outer_loop, inner_loop, filehandle_correctness, filehandle_debug, &testname[0][0], local_comm_all2all );

#endif

//! ====================================================================
//! =============================== outro ==============================
//! ====================================================================

//! close the filehandle for the output of the timing
  if ( myrank == 0 ) {
    timing_close_file();
#if TEST_TYPE > 1
    cleanup_papi();
#endif
  }

  MPI_Comm_free( &local_comm_all2all );
  if ( myrank < 2 ) {
    MPI_Comm_free( &local_comm_pp );
  }

  MPI_Barrier( MPI_COMM_WORLD );

  MPI_Finalize();

  return 0;
}
