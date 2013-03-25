void timing_close_file();
void timing_init( char* ptestname, char* pmethod, int pbytes );
void timing_open_file( char* filename );
void timing_print( int last );
void timing_record( int id );
void timing_set_max_tests(int value);
#if TEST_TYPE > 1
  void init_papi();
  void cleanup_papi();
#endif
void timing_hrt_init();

void timing_basic_ping_pong_nelements( int DIM1, int loop, char* testname, MPI_Comm local_communicator);
void timing_basic_alltoall_nelements( int DIM1, int procs, int loop, char* testname, MPI_Comm local_communicator);

void timing_fft2d_ddt( int DIM1, int procs, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator );
void timing_fft2d_manual( int DIM1, int procs, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator );
void timing_fft2d_mpi_pack_ddt( int DIM1, int procs, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator );

void timing_lammps_atomic_ddt( int DIM1, int icount, int* list, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator);
void timing_lammps_atomic_manual( int DIM1, int icount, int* list, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator);
void timing_lammps_atomic_mpi_pack_ddt( int DIM1, int icount, int* list, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator );

void timing_lammps_full_ddt( int DIM1, int icount, int* list, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator);
void timing_lammps_full_manual( int DIM1, int icount, int* list, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator);
void timing_lammps_full_mpi_pack_ddt( int DIM1, int icount, int* list, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator);

void timing_milc_su3_zdown_ddt( int DIM2, int DIM3, int DIM4, int DIM5, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator );
void timing_milc_su3_zdown_manual( int DIM2, int DIM3, int DIM4, int DIM5, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator );
void timing_milc_su3_zdown_mpi_pack_ddt( int DIM2, int DIM3, int DIM4, int DIM5, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator );

void timing_nas_lu_x_ddt( int DIM2, int DIM3, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator );
void timing_nas_lu_x_manual( int DIM2, int DIM3, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator );
void timing_nas_lu_x_mpi_pack_ddt( int DIM2, int DIM3, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator );
void timing_nas_lu_y_ddt( int DIM2, int DIM3, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator );
void timing_nas_lu_y_manual( int DIM2, int DIM3, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator );
void timing_nas_lu_y_mpi_pack_ddt( int DIM2, int DIM3, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator );

void timing_nas_mg_x_ddt( int DIM1, int DIM2, int DIM3, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator );
void timing_nas_mg_x_manual( int DIM1, int DIM2, int DIM3, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator );
void timing_nas_mg_x_mpi_pack_ddt( int DIM1, int DIM2, int DIM3, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator );
void timing_nas_mg_y_ddt( int DIM1, int DIM2, int DIM3, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator );
void timing_nas_mg_y_manual( int DIM1, int DIM2, int DIM3, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator );
void timing_nas_mg_y_mpi_pack_ddt( int DIM1, int DIM2, int DIM3, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator );
void timing_nas_mg_z_ddt( int DIM1, int DIM2, int DIM3, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator );
void timing_nas_mg_z_manual( int DIM1, int DIM2, int DIM3, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator );
void timing_nas_mg_z_mpi_pack_ddt( int DIM1, int DIM2, int DIM3, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator );

void timing_specfem3D_cm_ddt( int DIM2_cm, int DIM2_ic, int icount_cm, int icount_ic, int* list_cm, int* list_ic, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator );
void timing_specfem3D_cm_manual( int DIM2_cm, int DIM2_ic, int icount_cm, int icount_ic, int* list_cm, int* list_ic, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator );
void timing_specfem3D_cm_mpi_pack_ddt( int DIM2_cm, int DIM2_ic, int icount_cm, int icount_ic, int* list_cm, int* list_ic, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator );
void timing_specfem3d_mt_ddt( int DIM1, int DIM2, int DIM3, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator );
void timing_specfem3d_mt_manual( int DIM1, int DIM2, int DIM3, int outer_loop, int inner_loop, int *correct_flag, int *ptypesize, char *testname, MPI_File filehandle_debug, MPI_Comm local_communicator );
void timing_specfem3d_mt_mpi_pack_ddt( int DIM1, int DIM2, int DIM3, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator );
void timing_specfem3D_oc_ddt( int DIM1, int icount, int* list, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator );
void timing_specfem3D_oc_manual( int DIM1, int icount, int* list, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator );
void timing_specfem3D_oc_mpi_pack_ddt( int DIM1, int icount, int* list, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator );

void timing_wrf_manual ( int number_2D, int number_3D, int number_4D, int ims, int ime, int jms, int jme, int kms, int kme, int* limit_4D_arrays, int is, int ie, int js, int je,
  int ks, int ke, int param_first_scalar, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator );
void timing_wrf_sa_ddt ( int number_2D, int number_3D, int number_4D, int ims, int ime, int jms, int jme, int kms, int kme, int* limit_4D_arrays, int is, int ie, int js, int je,
  int ks, int ke, int param_first_scalar, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator );
void timing_wrf_sa_mpi_pack_ddt ( int number_2D, int number_3D, int number_4D, int ims, int ime, int jms, int jme, int kms, int kme, int* limit_4D_arrays, int is, int ie, int js, int je,
  int ks, int ke, int param_first_scalar, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator );
void timing_wrf_vec_ddt ( int number_2D, int number_3D, int number_4D, int ims, int ime, int jms, int jme, int kms, int kme, int* limit_4D_arrays, int is, int ie, int js, int je,
  int ks, int ke, int param_first_scalar, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator );
void timing_wrf_vec_mpi_pack_ddt ( int number_2D, int number_3D, int number_4D, int ims, int ime, int jms, int jme, int kms, int kme, int* limit_4D_arrays, int is, int ie, int js,
  int je, int ks, int ke, int param_first_scalar, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator );

void utilities_fill_unique_array_1D_float( float* array, int DIM1, int base );
void utilities_fill_unique_array_2D_float( float* array, int DIM1, int DIM2, int base );
void utilities_fill_unique_array_3D_float( float* array, int DIM1, int DIM2, int DIM3, int base );
void utilities_fill_unique_array_4D_float( float* array, int DIM1, int DIM2, int DIM3, int DIM4, int base );
void utilities_fill_unique_array_5D_float( float* array, int DIM1, int DIM2, int DIM3, int DIM4, int DIM5, int base );
void utilities_fill_unique_array_1D_double( double* array, int DIM1, int base );
void utilities_fill_unique_array_2D_double( double* array, int DIM1, int DIM2, int base );
void utilities_fill_unique_array_3D_double( double* array, int DIM1, int DIM2, int DIM3, int base );
void utilities_random_array_shuffle( int* index_list, int list_dim, int global_dim );

void wrapper_timing_fft( int DIM1, int outer_loop, int inner_loop, MPI_File filehandle_correctness, MPI_File filehandle_debug, char* ptestname, MPI_Comm local_communicator );
void wrapper_timing_lammps_atomic( int DIM1, int icount, int outer_loop, int inner_loop, MPI_File filehandle_correctness, MPI_File filehandle_debug, char* ptestname, MPI_Comm local_communicator);
void wrapper_timing_lammps_full( int DIM1, int icount, int outer_loop, int inner_loop, MPI_File filehandle_correctness, MPI_File filehandle_debug, char* ptestname, MPI_Comm local_communicator );
void wrapper_timing_milc_su3_zdown( int DIM2, int DIM3, int DIM4, int DIM5, int outer_loop, int inner_loop, MPI_File filehandle_correctness, MPI_File filehandle_debug, char* ptestname, MPI_Comm local_communicator );
void wrapper_timing_nas_lu( int DIM2, int DIM3, int outer_loop, int inner_loop, MPI_File filehandle_correctness, MPI_File filehandle_debug, char* ptestname, MPI_Comm local_communicator );
void wrapper_timing_nas_mg( int DIM1, int DIM2, int DIM3, int outer_loop, int inner_loop, MPI_File filehandle_correctness, MPI_File filehandle_debug, char* ptestname, MPI_Comm local_communicator );
void wrapper_timing_specfem3D_cm( int DIM2_cm, int DIM2_ic, int icount_cm, int icount_ic, int outer_loop, int inner_loop, MPI_File filehandle_correctness, MPI_File filehandle_debug, char* ptestname, MPI_Comm local_communicator );
void wrapper_timing_specfem3d_mt( int DIM1, int DIM2, int DIM3, int outer_loop, int inner_loop, MPI_File filehandle_correctness, MPI_File filehandle_debug, char* ptestname, MPI_Comm local_communicator );
void wrapper_timing_specfem3D_oc( int DIM1, int icount, int outer_loop, int inner_loop, MPI_File filehandle_correctness, MPI_File filehandle_debug, char* ptestname, MPI_Comm local_communicator);
void wrapper_timing_wrf( int number_2D, int number_3D, int number_4D, int ims, int ime, int jms, int jme, int kms, int kme, int* limit_4D_arrays, int is, int ie, int js, int je,
  int ks, int ke, int param_first_scalar, int outer_loop, int inner_loop, MPI_File filehandle_correctness, MPI_File filehandle_debug, char* ptestname, MPI_Comm local_communicator );
