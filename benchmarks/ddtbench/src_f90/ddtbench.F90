      program ddtbench

!      use timing, only: timing_open_file, timing_close_file, &
!        timing_set_max_tests
      use utilities, only: init_random_seed
      use wrapper, only: wrapper_timing_nas_mg, wrapper_timing_nas_lu, &
        wrapper_timing_wrf, wrapper_timing_milc_su3_zdown, &
        wrapper_timing_specfem3D_oc, wrapper_timing_specfem3D_cm, &
        wrapper_timing_lammps_full, wrapper_timing_lammps_atomic, &
        wrapper_timing_fft, wrapper_timing_specfem3d_mt

      implicit none

      include 'mpif.h'

      integer :: outer_loop = 10
      integer :: inner_loop = 20

      integer :: myrank, ier
      integer :: max_epochs
      integer :: local_comm_pp, local_comm_all2all
      
!> variables for file handling
      character(50) :: filename
      character(50) :: testname(3) 
      character :: newline
      integer :: filehandle_correctness, filehandle_debug

      integer :: argc, arg_status, arg_length, arg_number
      character(32) :: arg

      integer :: DIM1, DIM2, DIM3, DIM4

      integer :: number_2D, number_3D, number_4D
      integer :: ims, ime, jms, jme, kms, kme
      integer :: is, ie, js, je, ks, ke
      integer :: param_first_scalar
      integer, dimension(:), allocatable :: limit_4D_arrays

      integer :: icount1, icount2

!> ====================================================================
!> =============================== intro ==============================
!> ====================================================================

!> init some variables for file handling
      newline = char(10)
      testname = ""
      filehandle_correctness = MPI_FILE_NULL
      filehandle_debug = MPI_FILE_NULL

      call MPI_Init( ier )

#if TEST_TYPE != 2
      call timing_hrt_init()
#endif

      call MPI_Comm_rank( MPI_COMM_WORLD, myrank, ier )

      call MPI_Comm_dup( MPI_COMM_WORLD, local_comm_all2all, ier )

      if ( myrank .LT. 2 ) then
        call MPI_Comm_split( MPI_COMM_WORLD, 0, myrank, local_comm_pp, &
          ier )
      else
        call MPI_Comm_split( MPI_COMM_WORLD, MPI_UNDEFINED, myrank, &
          local_comm_pp, ier )
      endif

      call init_random_seed()

!> reading the arguments
!> only in compliance with the fortran 2003 standard
      argc = command_argument_count()

      if ( argc .NE. 0 ) then
        if ( argc .NE. 2 ) then
          if ( myrank  .EQ. 0 ) then
            write (*,*) "usage ddtbench <outer_loop> <inner_loop>"
            write (*,*) "proceed with standard parameter for the ", &
              "inner and outer loop", newline
          endif
        else
          arg_number = 1
          arg_length = 32
          call get_command_argument( arg_number, arg, arg_length, &
            arg_status );
          if ( arg_status .EQ. 0 ) then
            read (arg,*) outer_loop
            arg_number = 2
            call get_command_argument( arg_number, arg, arg_length, &
              arg_status );
            if ( arg_status .EQ. 0 ) then
              read (arg,*) inner_loop
            else
              if ( myrank .EQ. 0 ) then
                write (*,*) "could not retrieve the second argument"
                write (*,*) "proceed with standard parameter for the", &
                  " inner and outer loop", newline
              endif
            endif
          else
            if ( myrank .EQ. 0 ) then
              write (*,*) "could not retrieve the first argument"
              write (*,*) "proceed with standard parameter for the ", &
                "inner and outer loop", newline
            endif
          endif
        endif
      endif

!> some intro
      if ( myrank .EQ. 0 ) then
        write (*,*) "Welcome to our DDT benchmark suite", newline
        write (*,'(A,I0,A)') " outer loop parameter = ", outer_loop
        write (*,'(A,I0,A,A)') " inner loop parameter = ", inner_loop, &
          newline
!> open the filehandle for the output of the timing
        write (filename,'(A)') "ddtbench.out"
        call timing_open_file( filename )
#if TEST_TYPE > 1
        call init_papi()
#endif
!> set the maximum test number for the progression bar
!> ddt has 1 epoch in inner loop, 2 epochs in outer loop, and 2 epochs 
!> for each test
!> manual has 3 epochs in inner loop, 2 epochs in outer loop
!> mpi_pack_ddt has 3 epochs in inner loop, 2 epochs in outer looper
!> reference has 1 epoch in loops (=outer_loops*inner_loops) and 4
!> epochs for each test
!> there are 65 tests in sum
!> the SPECFEM3D_mt tests need to be treaded differently, because in the
!> manual pack routine there is no unpack
        max_epochs = 65 * (outer_loop * (inner_loop * (1+3+3+1) + &
         (2+2+2+0)) + (2+0+0+4)) + 4 * (outer_loop * (inner_loop * &
         (1+2+3+1) + (2+2+2+0)) + (2+1+0+4))
        call timing_set_max_tests( max_epochs )
      endif

      if ( myrank .LT. 2 ) then

!> ====================================================================
!> ========================== WRF y direction =========================
!> ====================================================================

      number_2D = 4
      number_3D = 3
      number_4D = 2

      allocate( limit_4D_arrays(number_4D), stat=ier )

      write (testname(1),'(A)') "WRF_y_vec"
      write (testname(2),'(A)') "WRF_y_sa"

! 2x2 ym send, em_b_wave case
      ims = -4
      ime = 27
      kms = 1
      kme = 65
      jms = 34
      jme = 85
      limit_4D_arrays(1:number_4D) = 2

      is = 1
      ie = 23
      ks = 1
      ke = 65
      js = 41
      je = 43
      param_first_scalar = 2

      call wrapper_timing_wrf( number_2D, number_3D, number_4D, &
        ims, ime, jms, jme, kms, kme, limit_4D_arrays, &
        is, ie, js, je, ks, ke, param_first_scalar, outer_loop, &
        inner_loop, filehandle_correctness, filehandle_debug, &
        testname, local_comm_pp )

! 3x3 ym send, em_b_wave case
      ims = 8
      ime = 34
      kms = 1
      kme = 65
      jms = 21
      jme = 61
      limit_4D_arrays(1:number_4D) = 2

      is = 12
      ie = 30
      ks = 1
      ke = 65
      js = 28
      je = 30
      param_first_scalar = 2

      call wrapper_timing_wrf( number_2D, number_3D, number_4D, &
        ims, ime, jms, jme, kms, kme, limit_4D_arrays, &
        is, ie, js, je, ks, ke, param_first_scalar, outer_loop, &
        inner_loop, filehandle_correctness, filehandle_debug, &
        testname, local_comm_pp )

! 4x4 ym send, em_b_wave case
      ims = 4
      ime = 27
      kms = 1
      kme = 65
      jms = 14
      jme = 47
      limit_4D_arrays(1:number_4D) = 2

      is = 8
      ie = 23
      ks = 1
      ke = 65
      js = 21
      je = 23
      param_first_scalar = 2

      call wrapper_timing_wrf( number_2D, number_3D, number_4D, &
        ims, ime, jms, jme, kms, kme, limit_4D_arrays, &
        is, ie, js, je, ks, ke, param_first_scalar, outer_loop, &
        inner_loop, filehandle_correctness, filehandle_debug, &
        testname, local_comm_pp )

! 5x5 ym send, em_b_wave case
      ims = 10
      ime = 31
      kms = 1
      kme = 65
      jms = 10
      jme = 39
      limit_4D_arrays(1:number_4D) = 2

      is = 14
      ie = 27
      ks = 1
      ke = 65
      js = 17
      je = 19
      param_first_scalar = 2

      call wrapper_timing_wrf( number_2D, number_3D, number_4D, &
        ims, ime, jms, jme, kms, kme, limit_4D_arrays, &
        is, ie, js, je, ks, ke, param_first_scalar, outer_loop, &
        inner_loop, filehandle_correctness, filehandle_debug, &
        testname, local_comm_pp )

! 8x8 ym send, em_b_wave case
      ims = 4
      ime = 22
      kms = 1
      kme = 65
      jms = 4
      jme = 27
      limit_4D_arrays(1:number_4D) = 2

      is = 8
      ie = 18
      ks = 1
      ke = 65
      js = 11
      je = 13
      param_first_scalar = 2

      call wrapper_timing_wrf( number_2D, number_3D, number_4D, &
        ims, ime, jms, jme, kms, kme, limit_4D_arrays, &
        is, ie, js, je, ks, ke, param_first_scalar, outer_loop, &
        inner_loop, filehandle_correctness, filehandle_debug, &
        testname, local_comm_pp )

!> ====================================================================
!> ========================== WRF x direction =========================
!> ====================================================================

      write (testname(1),'(A)') "WRF_x_vec"
      write (testname(2),'(A)') "WRF_x_sa"

! 2x2 xp send, em_b_wave
      ims = -4
      ime = 27
      kms = 1
      kme = 65
      jms = 34
      jme = 85
      limit_4D_arrays(1:number_4D) = 2

      is = 18
      ie = 20
      ks = 1
      ke = 65
      js = 38
      je = 81
      param_first_scalar = 2

      call wrapper_timing_wrf( number_2D, number_3D, number_4D, &
        ims, ime, jms, jme, kms, kme, limit_4D_arrays, &
        is, ie, js, je, ks, ke, param_first_scalar, outer_loop, &
        inner_loop, filehandle_correctness, filehandle_debug, &
        testname, local_comm_pp )

! 4x4 xp send, em_b_wave
      ims = 4
      ime = 27
      kms = 1
      kme = 65
      jms = 14
      jme = 47
      limit_4D_arrays(1:number_4D) = 2

      is = 18
      ie = 20
      ks = 1
      ke = 65
      js = 18
      je = 43
      param_first_scalar = 2

      call wrapper_timing_wrf( number_2D, number_3D, number_4D, &
        ims, ime, jms, jme, kms, kme, limit_4D_arrays, &
        is, ie, js, je, ks, ke, param_first_scalar, outer_loop, &
        inner_loop, filehandle_correctness, filehandle_debug, &
        testname, local_comm_pp )

! 8x8 xp send, em_b_wave
      ims = 4
      ime = 22
      kms = 1
      kme = 65
      jms = 4
      jme = 27
      limit_4D_arrays(1:number_4D) = 2

      is = 13
      ie = 15
      ks = 1
      ke = 65
      js = 8
      je = 23
      param_first_scalar = 2

      call wrapper_timing_wrf( number_2D, number_3D, number_4D, &
        ims, ime, jms, jme, kms, kme, limit_4D_arrays, &
        is, ie, js, je, ks, ke, param_first_scalar, outer_loop, &
        inner_loop, filehandle_correctness, filehandle_debug, &
        testname, local_comm_pp )

      deallocate( limit_4D_arrays, stat=ier )

      write (testname(1),'(A)') ""
      write (testname(2),'(A)') ""

!> ====================================================================
!> ==================== MILC su3 zdown direction ======================
!> ====================================================================

      DIM1 = 16
      DIM2 = 16
      DIM3 = 16
      DIM4 = 16

      call wrapper_timing_milc_su3_zdown( DIM1, DIM2, DIM3, &
        DIM4, outer_loop, inner_loop, filehandle_correctness, &
        filehandle_debug, testname(1), local_comm_pp )
 
      DIM1 = 8
      DIM2 = 8
      DIM3 = 16
      DIM4 = 16

      call wrapper_timing_milc_su3_zdown( DIM1, DIM2, DIM3, &
        DIM4, outer_loop, inner_loop, filehandle_correctness, &
        filehandle_debug, testname(1), local_comm_pp )
  
      DIM1 = 8
      DIM2 = 8
      DIM3 = 8
      DIM4 = 16

      call wrapper_timing_milc_su3_zdown( DIM1, DIM2, DIM3, &
        DIM4, outer_loop, inner_loop, filehandle_correctness, &
        filehandle_debug, testname(1), local_comm_pp )

      DIM1 = 8
      DIM2 = 8
      DIM3 = 8
      DIM4 = 8

      call wrapper_timing_milc_su3_zdown( DIM1, DIM2, DIM3, &
        DIM4, outer_loop, inner_loop, filehandle_correctness, &
        filehandle_debug, testname(1), local_comm_pp )

!> ====================================================================
!> =============================== NAS_LU =============================
!> ====================================================================

      DIM2 = 12
      DIM3 = 12

      call wrapper_timing_nas_lu( DIM2, DIM3, outer_loop, &
        inner_loop, filehandle_correctness, filehandle_debug, &
        testname(1), local_comm_pp )

      DIM2 = 33
      DIM3 = 33

      call wrapper_timing_nas_lu( DIM2, DIM3, outer_loop, &
        inner_loop, filehandle_correctness, filehandle_debug, &
        testname(1), local_comm_pp )

      DIM2 = 64
      DIM3 = 64

      call wrapper_timing_nas_lu( DIM2, DIM3, outer_loop, &
        inner_loop, filehandle_correctness, filehandle_debug, &
        testname(1), local_comm_pp )

      DIM2 = 102
      DIM3 = 102

      call wrapper_timing_nas_lu( DIM2, DIM3, outer_loop, &
        inner_loop, filehandle_correctness, filehandle_debug, &
        testname(1), local_comm_pp )

      DIM2 = 162
      DIM3 = 162

      call wrapper_timing_nas_lu( DIM2, DIM3, outer_loop, &
        inner_loop, filehandle_correctness, filehandle_debug, &
        testname(1), local_comm_pp )

      DIM2 = 408
      DIM3 = 408

      call wrapper_timing_nas_lu( DIM2, DIM3, outer_loop, &
        inner_loop, filehandle_correctness, filehandle_debug, &
        testname(1), local_comm_pp )

      DIM2 = 1020
      DIM3 = 1020

      call wrapper_timing_nas_lu( DIM2, DIM3, outer_loop, &
        inner_loop, filehandle_correctness, filehandle_debug, &
        testname(1), local_comm_pp )

!> ====================================================================
!> =============================== NAS_MG =============================
!> ====================================================================

      DIM1 = 34
      DIM2 = 18
      DIM3 = 18

      call wrapper_timing_nas_mg( DIM1, DIM2, DIM3, &
        outer_loop, inner_loop, filehandle_correctness, &
        filehandle_debug, testname, local_comm_pp )

      DIM1 = 130
      DIM2 = 66
      DIM3 = 66

      call wrapper_timing_nas_mg( DIM1, DIM2, DIM3, &
        outer_loop, inner_loop, filehandle_correctness, &
        filehandle_debug, testname, local_comm_pp )

      DIM1 = 258
      DIM2 = 130
      DIM3 = 130

      call wrapper_timing_nas_mg( DIM1, DIM2, DIM3, &
        outer_loop, inner_loop, filehandle_correctness, &
        filehandle_debug, testname, local_comm_pp )

      DIM1 = 514
      DIM2 = 258
      DIM3 = 258

      call wrapper_timing_nas_mg( DIM1, DIM2, DIM3, &
        outer_loop, inner_loop, filehandle_correctness, &
        filehandle_debug, testname, local_comm_pp )

!> ====================================================================
!> ============================ LAMMPS_full ===========================
!> ====================================================================

!> peptide example with 2 process (maximum)
      icount1 = 3062
      DIM1 = 3534

      call wrapper_timing_lammps_full( DIM1, icount1, outer_loop, &
        inner_loop, filehandle_correctness, filehandle_debug, &
        testname(1), local_comm_pp )

!> peptide example with 4 process (maximum)
      icount1 = 2243
      DIM1 = 2597

      call wrapper_timing_lammps_full( DIM1, icount1, outer_loop, &
        inner_loop, filehandle_correctness, filehandle_debug, &
        testname(1), local_comm_pp )

!> peptide example with 8 process (maximum)
      icount1 = 1662
      DIM1 = 1907

      call wrapper_timing_lammps_full( DIM1, icount1, outer_loop, &
        inner_loop, filehandle_correctness, filehandle_debug, &
        testname(1), local_comm_pp )

!> ====================================================================
!> =========================== LAMMPS_atomic ==========================
!> ====================================================================

!> crack example with 2 process (maximum)
      icount1 = 243
      DIM1 = 4084

      call wrapper_timing_lammps_atomic( DIM1, icount1, outer_loop, &
        inner_loop, filehandle_correctness, filehandle_debug, &
        testname(1), local_comm_pp )

!> crack example with 4 process (maximum)
      icount1 = 145
      DIM1 = 2157

      call wrapper_timing_lammps_atomic( DIM1, icount1, outer_loop, &
        inner_loop, filehandle_correctness, filehandle_debug, &
        testname(1), local_comm_pp )

!> crack example with 8 process (maximum)
      icount1 = 114
      DIM1 = 1217

      call wrapper_timing_lammps_atomic( DIM1, icount1, outer_loop, &
        inner_loop, filehandle_correctness, filehandle_debug, &
        testname(1), local_comm_pp )

!> ====================================================================
!> ============================ SPECFEM3D_oc ==========================
!> ====================================================================

!> 10x10x6, c=4
      icount1 = 3225
      DIM1 = 88881

      call wrapper_timing_specfem3D_oc( DIM1, icount1, outer_loop, &
        inner_loop, filehandle_correctness, filehandle_debug, &
        testname(1), local_comm_pp )

!> 10x10x6, c=3
      icount1 = 1897
      DIM1 = 38585

      call wrapper_timing_specfem3D_oc( DIM1, icount1, outer_loop, &
        inner_loop, filehandle_correctness, filehandle_debug, &
        testname(1), local_comm_pp )

!> 10x10x6, c=2
      icount1 = 877
      DIM1 = 12857

      call wrapper_timing_specfem3D_oc( DIM1, icount1, outer_loop, &
        inner_loop, filehandle_correctness, filehandle_debug, &
        testname(1), local_comm_pp )

!> 10x10x6, c=1
      icount1 = 493
      DIM1 = 3697

      call wrapper_timing_specfem3D_oc( DIM1, icount1, outer_loop, &
        inner_loop, filehandle_correctness, filehandle_debug, &
        testname(1), local_comm_pp )

!> ====================================================================
!> ============================ SPECFEM3D_cm ==========================
!> ====================================================================

!> 10x10x6, c=4
      icount1 = 11797
      icount2 = 3009
      DIM1 = 834917
      DIM2 = 51153

      call wrapper_timing_specfem3D_cm( DIM1, DIM2, icount1, icount2, &
        outer_loop, inner_loop, filehandle_correctness, &
        filehandle_debug, testname(1), local_comm_pp )

!> 10x10x6, c=3
      icount1 = 7125
      icount2 = 1729
      DIM1 = 396849
      DIM2 = 22477

      call wrapper_timing_specfem3D_cm( DIM1, DIM2, icount1, icount2, &
        outer_loop, inner_loop, filehandle_correctness, &
        filehandle_debug, testname(1), local_comm_pp )

!> 10x10x6, c=2
      icount1 = 3801
      icount2 = 801
      DIM1 = 152001
      DIM2 = 7209

      call wrapper_timing_specfem3D_cm( DIM1, DIM2, icount1, icount2, &
        outer_loop, inner_loop, filehandle_correctness, &
        filehandle_debug, testname(1), local_comm_pp )

!> 10x10x6, c=1
      icount1 = 1957
      icount2 = 245
      DIM1 = 39929
      DIM2 = 1225

      call wrapper_timing_specfem3D_cm( DIM1, DIM2, icount1, icount2, &
        outer_loop, inner_loop, filehandle_correctness, &
        filehandle_debug, testname(1), local_comm_pp )

!> ====================================================================
!> ============================ SPECFEM3D_mt ==========================
!> ====================================================================

!> 10x10x6, c=1
      DIM1 = 3
      DIM2 = 2
      DIM3 = 7600

      call wrapper_timing_specfem3d_mt( DIM1, DIM2, DIM3, &
        outer_loop, inner_loop, filehandle_correctness, &
        filehandle_debug, testname(1), local_comm_pp )

!> 10x10x6, c=2
      DIM1 = 3
      DIM2 = 2
      DIM3 = 6400

      call wrapper_timing_specfem3d_mt( DIM1, DIM2, DIM3, &
        outer_loop, inner_loop, filehandle_correctness, &
        filehandle_debug, testname(1), local_comm_pp )

!> 10x10x6, c=3
      DIM1 = 3
      DIM2 = 2
      DIM3 = 5600

      call wrapper_timing_specfem3d_mt( DIM1, DIM2, DIM3, &
        outer_loop, inner_loop, filehandle_correctness, &
        filehandle_debug, testname(1), local_comm_pp )

!> 10x10x6, c=4
      DIM1 = 3
      DIM2 = 2
      DIM3 = 6200

      call wrapper_timing_specfem3d_mt( DIM1, DIM2, DIM3, &
        outer_loop, inner_loop, filehandle_correctness, &
        filehandle_debug, testname(1), local_comm_pp )


      endif !> of myrank < 2

!> ====================================================================
!> ================================ FFT ===============================
!> ====================================================================

      call MPI_Barrier( local_comm_all2all, ier )

      DIM1 = 256
      call wrapper_timing_fft( DIM1, outer_loop, inner_loop, &
        filehandle_correctness, filehandle_debug, testname(1), &
        local_comm_all2all )

      DIM1 = 512
      call wrapper_timing_fft( DIM1, outer_loop, inner_loop, &
        filehandle_correctness, filehandle_debug, testname(1), &
        local_comm_all2all )

      DIM1 = 768
      call wrapper_timing_fft( DIM1, outer_loop, inner_loop, &
        filehandle_correctness, filehandle_debug, testname(1), &
        local_comm_all2all )

      DIM1 = 1024
      call wrapper_timing_fft( DIM1, outer_loop, inner_loop, &
        filehandle_correctness, filehandle_debug, testname(1), &
        local_comm_all2all )

      DIM1 = 1536
      call wrapper_timing_fft( DIM1, outer_loop, inner_loop, &
        filehandle_correctness, filehandle_debug, testname(1), &
        local_comm_all2all )

!> ====================================================================
!> =============================== outro ==============================
!> ====================================================================

!> close the filehandle for the output of the timing
      if ( myrank .EQ. 0 ) then
        call timing_close_file()
#if TEST_TYPE > 1
        call cleanup_papi()
#endif
      endif

      call MPI_Comm_free( local_comm_all2all, ier )
      if ( myrank .LT. 2 ) then
        call MPI_Comm_free( local_comm_pp, ier )
      endif

      call MPI_Finalize( ier )

      end program ddtbench
