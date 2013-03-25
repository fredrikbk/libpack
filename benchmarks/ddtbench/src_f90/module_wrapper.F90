      module wrapper
!> contains a wrapper for each test, so that the main program don't need
!> to call every method by himself
!> also this is the part, where the parameter checking should be put
!> also the output for the correctness checking should be placed here

      private

      private :: wrapper_timing_wrf_vec
      private :: wrapper_timing_nas_lu_x
      private :: wrapper_timing_nas_lu_y
      private :: wrapper_timing_nas_mg_x
      private :: wrapper_timing_nas_mg_y
      private :: wrapper_timing_nas_mg_z
  
      public :: wrapper_timing_wrf     
      public :: wrapper_timing_milc_su3_zdown
      public :: wrapper_timing_nas_lu
      public :: wrapper_timing_nas_mg
      public :: wrapper_timing_fft
      public :: wrapper_timing_specfem3D_mt
      public :: wrapper_timing_lammps_full
      public :: wrapper_timing_lammps_atomic
      public :: wrapper_timing_specfem3D_oc
      public :: wrapper_timing_specfem3D_cm

      contains

      subroutine wrapper_timing_wrf_vec( number_2D, number_3D, &
        number_4D, ims, ime, jms, jme, kms, kme, limit_4D_arrays, &
        is, ie, js, je, ks, ke, param_first_scalar, outer_loop, &
        inner_loop, filehandle_correctness, filehandle_debug, &
        ptestname, local_communicator )

      use timing_basic, only: time_ping_pong_nelements
      use timing_wrf, only: timing_wrf_vec_ddt, timing_wrf_manual, &
        timing_wrf_vec_mpi_pack_ddt

      implicit none

      integer, intent(in) :: number_2D, number_3D, number_4D
      integer, intent(in) :: ims, ime, jms, jme, kms, kme
      integer, intent(in) :: is, ie, js, je, ks, ke

      integer, intent(in) :: param_first_scalar
      integer, dimension(number_4D), intent(in) :: limit_4D_arrays

      integer, intent(in) :: outer_loop, inner_loop

      integer, intent(in) :: filehandle_correctness, filehandle_debug

      character(50), intent(in) :: ptestname

      integer, intent(in) :: local_communicator

!> local variables
      logical :: correct_flag
      integer :: typesize
      character(50) :: testname
      integer :: nelements, loops
      integer :: m

!> just some statements to prevent compiler warnings of unused variables
!> those parameter are included for future features
      typesize = filehandle_correctness

      if ( len(trim(ptestname)) .EQ. 0 ) then
        write (testname,'(A)') "WRF_vec"
      else
        testname = ptestname
      endif

      call timing_wrf_vec_ddt ( number_2D, number_3D, number_4D, &
        ims, ime, jms, jme, kms, kme, limit_4D_arrays, &
        is, ie, js, je, ks, ke, param_first_scalar, &
        outer_loop, inner_loop, correct_flag, typesize, testname, &
        filehandle_debug, local_communicator )

      call timing_wrf_manual ( number_2D, number_3D, number_4D, &
        ims, ime, jms, jme, kms, kme, limit_4D_arrays, &
        is, ie, js, je, ks, ke, param_first_scalar, &
        outer_loop, inner_loop, correct_flag, typesize, testname, &
        filehandle_debug, local_communicator )

      call timing_wrf_vec_mpi_pack_ddt ( number_2D, number_3D, &
        number_4D, ims, ime, jms, jme, kms, kme, limit_4D_arrays, &
        is, ie, js, je, ks, ke, param_first_scalar, &
        outer_loop, inner_loop, correct_flag, typesize, testname, &
        filehandle_debug, local_communicator )

      nelements = number_2D * (ie-is+1) * (je-js+1) + &
        number_3D * (ie-is+1) * (je-js+1) * (ke-ks+1)
      do m=1,number_4D
        if (limit_4D_arrays(m) .GE. param_first_scalar) then
          nelements = nelements + (limit_4D_arrays(m) - &
            param_first_scalar+1) * (ie-is+1) * (je-js+1) * (ke-ks+1)
        endif
      enddo
      loops = outer_loop * inner_loop
      call time_ping_pong_nelements( nelements, loops, testname, &
        local_communicator )

      end subroutine wrapper_timing_wrf_vec

      subroutine wrapper_timing_wrf_sa( number_2D, number_3D, &
        number_4D, ims, ime, jms, jme, kms, kme, limit_4D_arrays, &
        is, ie, js, je, ks, ke, param_first_scalar, outer_loop, &
        inner_loop, filehandle_correctness, filehandle_debug, &
        ptestname, local_communicator )

      use timing_basic, only: time_ping_pong_nelements
      use timing_wrf, only: timing_wrf_sa_ddt, timing_wrf_manual, &
        timing_wrf_sa_mpi_pack_ddt

      implicit none

      integer, intent(in) :: number_2D, number_3D, number_4D
      integer, intent(in) :: ims, ime, jms, jme, kms, kme
      integer, intent(in) :: is, ie, js, je, ks, ke

      integer, intent(in) :: param_first_scalar
      integer, dimension(number_4D), intent(in) :: limit_4D_arrays

      integer, intent(in) :: outer_loop, inner_loop

      integer, intent(in) :: filehandle_correctness, filehandle_debug

      character(50), intent(in) :: ptestname

      integer, intent(in) :: local_communicator

!> local variables
      logical :: correct_flag
      integer :: typesize
      character(50) :: testname
      integer :: nelements, loops
      integer :: m

!> just some statements to prevent compiler warnings of unused variables
!> those parameter are included for future features
      typesize = filehandle_correctness

      if ( len(trim(ptestname)) .EQ. 0 ) then
        write (testname,'(A)') "WRF_sa"
      else
        testname = ptestname
      endif

      call timing_wrf_sa_ddt ( number_2D, number_3D, number_4D, &
        ims, ime, jms, jme, kms, kme, limit_4D_arrays, &
        is, ie, js, je, ks, ke, param_first_scalar, &
        outer_loop, inner_loop, correct_flag, typesize, testname, &
        filehandle_debug, local_communicator )

      call timing_wrf_manual ( number_2D, number_3D, number_4D, &
        ims, ime, jms, jme, kms, kme, limit_4D_arrays, &
        is, ie, js, je, ks, ke, param_first_scalar, &
        outer_loop, inner_loop, correct_flag, typesize, testname, &
        filehandle_debug, local_communicator )

      call timing_wrf_sa_mpi_pack_ddt ( number_2D, number_3D, &
        number_4D, ims, ime, jms, jme, kms, kme, limit_4D_arrays, &
        is, ie, js, je, ks, ke, param_first_scalar, &
        outer_loop, inner_loop, correct_flag, typesize, testname, &
        filehandle_debug, local_communicator )

      nelements = number_2D * (ie-is+1) * (je-js+1) + &
        number_3D * (ie-is+1) * (je-js+1) * (ke-ks+1)
      do m=1,number_4D
        if (limit_4D_arrays(m) .GE. param_first_scalar) then
          nelements = nelements + (limit_4D_arrays(m) - &
            param_first_scalar+1) * (ie-is+1) * (je-js+1) * (ke-ks+1)
        endif
      enddo
      loops = outer_loop * inner_loop
      call time_ping_pong_nelements( nelements, loops, testname, &
        local_communicator )

      end subroutine wrapper_timing_wrf_sa

      subroutine wrapper_timing_wrf( number_2D, number_3D, &
        number_4D, ims, ime, jms, jme, kms, kme, limit_4D_arrays, &
        is, ie, js, je, ks, ke, param_first_scalar, outer_loop, &
        inner_loop, filehandle_correctness, filehandle_debug, &
        ptestname, local_communicator )

      implicit none

      integer, intent(in) :: number_2D, number_3D, number_4D
      integer, intent(in) :: ims, ime, jms, jme, kms, kme
      integer, intent(in) :: is, ie, js, je, ks, ke

      integer, intent(in) :: param_first_scalar
      integer, dimension(number_4D), intent(in) :: limit_4D_arrays

      integer, intent(in) :: outer_loop, inner_loop

      integer, intent(in) :: filehandle_correctness, filehandle_debug

      character(50), intent(in) :: ptestname(2)

      integer, intent(in) :: local_communicator

      call wrapper_timing_wrf_vec( number_2D, number_3D, &
        number_4D, ims, ime, jms, jme, kms, kme, limit_4D_arrays, &
        is, ie, js, je, ks, ke, param_first_scalar, outer_loop, &
        inner_loop, filehandle_correctness, filehandle_debug, &
        ptestname(1), local_communicator )

      call wrapper_timing_wrf_sa( number_2D, number_3D, &
        number_4D, ims, ime, jms, jme, kms, kme, limit_4D_arrays, &
        is, ie, js, je, ks, ke, param_first_scalar, outer_loop, &
        inner_loop, filehandle_correctness, filehandle_debug, &
        ptestname(2), local_communicator )

      end subroutine wrapper_timing_wrf

      subroutine wrapper_timing_milc_su3_zdown( DIM2, DIM3, DIM4, &
        DIM5, outer_loop, inner_loop, filehandle_correctness, &
        filehandle_debug, ptestname, local_communicator )

      use timing_basic, only: time_ping_pong_nelements
      use timing_milc, only: timing_milc_su3_zdown_ddt, &
        timing_milc_su3_zdown_manual, timing_milc_su3_zdown_mpi_pack_ddt

      implicit none

      integer, intent(in) :: DIM2, DIM3, DIM4, DIM5
      integer, intent(in) :: outer_loop, inner_loop

      integer, intent(in) :: filehandle_correctness, filehandle_debug

      character(50), intent(in) :: ptestname

      integer, intent(in) :: local_communicator

!> local variables
      logical :: correct_flag
      integer :: typesize
      character(50) :: testname
      integer :: nelements, loops

!> just some statements to prevent compiler warnings of unused variables
!> those parameter are included for future features
      typesize = filehandle_correctness

      if ( len(trim(ptestname)) .EQ. 0 ) then
        write (testname,'(A)') "MILC_su3_zd"
      else
        testname = ptestname
      endif

      call timing_milc_su3_zdown_ddt( DIM2, DIM3, DIM4, DIM5, &
        outer_loop, inner_loop, correct_flag, typesize, testname, &
        filehandle_debug, local_communicator )

      call timing_milc_su3_zdown_manual( DIM2, DIM3, DIM4, DIM5, &
        outer_loop, inner_loop, correct_flag, typesize, testname, &
        filehandle_debug, local_communicator )

      call timing_milc_su3_zdown_mpi_pack_ddt( DIM2, DIM3, DIM4, DIM5, &
        outer_loop, inner_loop, correct_flag, typesize, testname, &
        filehandle_debug, local_communicator )

!> not necessarily correct, since it assumes that a complex uses twice
!> the bytes a real does
      nelements = DIM2*DIM3/2*DIM5*2*3*2
      loops = outer_loop * inner_loop
      call time_ping_pong_nelements( nelements, loops, testname, &
        local_communicator )

      end subroutine wrapper_timing_milc_su3_zdown

      subroutine wrapper_timing_nas_lu_x( DIM2, DIM3, outer_loop, &
        inner_loop, filehandle_correctness, filehandle_debug, &
        ptestname, local_communicator )

      use timing_basic, only: time_ping_pong_nelements
      use timing_nas, only: timing_nas_lu_x_ddt, timing_nas_lu_x_manual, &
        timing_nas_lu_x_mpi_pack_ddt

      implicit none

      integer, intent(in) :: DIM2, DIM3
      integer, intent(in) :: outer_loop, inner_loop

      integer, intent(in) :: filehandle_correctness, filehandle_debug

      character(50), intent(in) :: ptestname

      integer, intent(in) :: local_communicator
!> local variables
      logical :: correct_flag
      integer :: typesize
      character(50) :: testname
      integer :: nelements, loops

!> just some statements to prevent compiler warnings of unused variables
!> those parameter are included for future features
      typesize = filehandle_correctness

      if ( len(trim(ptestname)) .EQ. 0 ) then
        write (testname,'(A)') "NAS_LU_x"
      else
        testname = ptestname
      endif

      call timing_nas_lu_x_ddt( DIM2, DIM3, outer_loop, &
        inner_loop, correct_flag, typesize, testname, filehandle_debug, &
        local_communicator )

      call timing_nas_lu_x_manual( DIM2, DIM3, outer_loop, &
        inner_loop, correct_flag, typesize, testname, filehandle_debug, &
        local_communicator )

      call timing_nas_lu_x_mpi_pack_ddt( DIM2, DIM3, outer_loop, &
        inner_loop, correct_flag, typesize, testname, filehandle_debug, &
        local_communicator )
!> not necessarily correct, since it assumes that a double uses twice
!> the bytes a real does
      nelements = 5 * DIM2 * 2
      loops = outer_loop * inner_loop
      call time_ping_pong_nelements( nelements, loops, testname, &
        local_communicator )

      end subroutine wrapper_timing_nas_lu_x

      subroutine wrapper_timing_nas_lu_y( DIM2, DIM3, outer_loop, &
        inner_loop, filehandle_correctness, filehandle_debug, &
        ptestname, local_communicator )

      use timing_basic, only: time_ping_pong_nelements
      use timing_nas, only: timing_nas_lu_y_ddt, timing_nas_lu_y_manual, &
        timing_nas_lu_y_mpi_pack_ddt

      implicit none

      integer, intent(in) :: DIM2, DIM3
      integer, intent(in) :: outer_loop, inner_loop

      integer, intent(in) :: filehandle_correctness, filehandle_debug

      character(50), intent(in) :: ptestname

      integer, intent(in) :: local_communicator
!> local variables
      logical :: correct_flag
      integer :: typesize
      character(50) :: testname
      integer :: nelements, loops

!> just some statements to prevent compiler warnings of unused variables
!> those parameter are included for future features
      typesize = filehandle_correctness

      if ( len(trim(ptestname)) .EQ. 0 ) then
        write (testname,'(A)') "NAS_LU_y"
      else
        testname = ptestname
      endif

      call timing_nas_lu_y_ddt( DIM2, DIM3, outer_loop, &
        inner_loop, correct_flag, typesize, testname, filehandle_debug, &
        local_communicator )

      call timing_nas_lu_y_manual( DIM2, DIM3, outer_loop, &
        inner_loop, correct_flag, typesize, testname, filehandle_debug, &
        local_communicator )

      call timing_nas_lu_y_mpi_pack_ddt( DIM2, DIM3, outer_loop, &
        inner_loop, correct_flag, typesize, testname, filehandle_debug, &
        local_communicator )
!> not necessarily correct, since it assumes that a double uses twice
!> the bytes a real does
      nelements = 5 * DIM3 * 2
      loops = outer_loop * inner_loop
      call time_ping_pong_nelements( nelements, loops, testname, &
        local_communicator )

      end subroutine wrapper_timing_nas_lu_y

      subroutine wrapper_timing_nas_lu( DIM2, DIM3, outer_loop, &
        inner_loop, filehandle_correctness, filehandle_debug, &
        ptestname, local_communicator )

      implicit none

      integer, intent(in) :: DIM2, DIM3
      integer, intent(in) :: outer_loop, inner_loop

      integer, intent(in) :: filehandle_correctness, filehandle_debug

      character(50), intent(in) :: ptestname(2)

      integer, intent(in) :: local_communicator

      call wrapper_timing_nas_lu_x( DIM2, DIM3, outer_loop, &
        inner_loop, filehandle_correctness, filehandle_debug, &
        ptestname(1), local_communicator )

      call wrapper_timing_nas_lu_y( DIM2, DIM3, outer_loop, &
        inner_loop, filehandle_correctness, filehandle_debug, &
        ptestname(2), local_communicator )

      end subroutine wrapper_timing_nas_lu

      subroutine wrapper_timing_nas_mg_x( DIM1, DIM2, DIM3, outer_loop, &
        inner_loop, filehandle_correctness, filehandle_debug, &
        ptestname, local_communicator )

      use timing_basic, only: time_ping_pong_nelements
      use timing_nas, only: timing_nas_mg_x_ddt, timing_nas_mg_x_manual, &
        timing_nas_mg_x_mpi_pack_ddt

      implicit none

      integer, intent(in) :: DIM1, DIM2, DIM3
      integer, intent(in) :: outer_loop, inner_loop

      integer, intent(in) :: filehandle_correctness, filehandle_debug

      character(50), intent(in) :: ptestname

      integer, intent(in) :: local_communicator
!> local variables
      logical :: correct_flag
      integer :: typesize
      character(50) :: testname
      integer :: nelements, loops

!> just some statements to prevent compiler warnings of unused variables
!> those parameter are included for future features
      typesize = filehandle_correctness

      if ( len(trim(ptestname)) .EQ. 0 ) then
        write (testname,'(A)') "NAS_MG_x"
      else
        testname = ptestname
      endif

      call timing_nas_mg_x_ddt( DIM1, DIM2, DIM3, outer_loop, &
        inner_loop, correct_flag, typesize, testname, filehandle_debug, &
        local_communicator )

      call timing_nas_mg_x_manual( DIM1, DIM2, DIM3, outer_loop, &
        inner_loop, correct_flag, typesize, testname, filehandle_debug, &
        local_communicator )

      call timing_nas_mg_x_mpi_pack_ddt( DIM1, DIM2, DIM3, outer_loop, &
        inner_loop, correct_flag, typesize, testname, filehandle_debug, &
        local_communicator )
!> not necessarily correct, since it assumes that a double uses twice
!> the bytes a real does
      nelements = (DIM2-2)*(DIM3-2)*2
      loops = outer_loop * inner_loop
      call time_ping_pong_nelements( nelements, loops, testname, &
        local_communicator )

      end subroutine wrapper_timing_nas_mg_x

      subroutine wrapper_timing_nas_mg_y( DIM1, DIM2, DIM3, outer_loop, &
        inner_loop, filehandle_correctness, filehandle_debug, &
        ptestname, local_communicator )

      use timing_basic, only: time_ping_pong_nelements
      use timing_nas, only: timing_nas_mg_y_ddt, timing_nas_mg_y_manual, &
        timing_nas_mg_y_mpi_pack_ddt

      implicit none

      integer, intent(in) :: DIM1, DIM2, DIM3
      integer, intent(in) :: outer_loop, inner_loop

      integer, intent(in) :: filehandle_correctness, filehandle_debug

      character(50), intent(in) :: ptestname

      integer, intent(in) :: local_communicator
!> local variables
      logical :: correct_flag
      integer :: typesize
      character(50) :: testname
      integer :: nelements, loops

!> just some statements to prevent compiler warnings of unused variables
!> those parameter are included for future features
      typesize = filehandle_correctness

      if ( len(trim(ptestname)) .EQ. 0 ) then
        write (testname,'(A)') "NAS_MG_y"
      else
        testname = ptestname
      endif

      call timing_nas_mg_y_ddt( DIM1, DIM2, DIM3, outer_loop, &
        inner_loop, correct_flag, typesize, testname, filehandle_debug, &
        local_communicator )

      call timing_nas_mg_y_manual( DIM1, DIM2, DIM3, outer_loop, &
        inner_loop, correct_flag, typesize, testname, filehandle_debug, &
        local_communicator )

      call timing_nas_mg_y_mpi_pack_ddt( DIM1, DIM2, DIM3, outer_loop, &
        inner_loop, correct_flag, typesize, testname, filehandle_debug, &
        local_communicator )
!> not necessarily correct, since it assumes that a double uses twice
!> the bytes a real does
      nelements = (DIM1-2)*(DIM3-2)*2
      loops = outer_loop * inner_loop
      call time_ping_pong_nelements( nelements, loops, testname, &
        local_communicator )

      end subroutine wrapper_timing_nas_mg_y

      subroutine wrapper_timing_nas_mg_z( DIM1, DIM2, DIM3, outer_loop, &
        inner_loop, filehandle_correctness, filehandle_debug, &
        ptestname, local_communicator )

      use timing_basic, only: time_ping_pong_nelements
      use timing_nas, only: timing_nas_mg_z_ddt, timing_nas_mg_z_manual, &
        timing_nas_mg_z_mpi_pack_ddt

      implicit none

      integer, intent(in) :: DIM1, DIM2, DIM3
      integer, intent(in) :: outer_loop, inner_loop

      integer, intent(in) :: filehandle_correctness, filehandle_debug

      character(50), intent(in) :: ptestname

      integer, intent(in) :: local_communicator
!> local variables
      logical :: correct_flag
      integer :: typesize
      character(50) :: testname
      integer :: nelements, loops

!> just some statements to prevent compiler warnings of unused variables
!> those parameter are included for future features
      typesize = filehandle_correctness

      if ( len(trim(ptestname)) .EQ. 0 ) then
        write (testname,'(A)') "NAS_MG_z"
      else
        testname = ptestname
      endif

      call timing_nas_mg_z_ddt( DIM1, DIM2, DIM3, outer_loop, &
        inner_loop, correct_flag, typesize, testname, filehandle_debug, &
        local_communicator )

      call timing_nas_mg_z_manual( DIM1, DIM2, DIM3, outer_loop, &
        inner_loop, correct_flag, typesize, testname, filehandle_debug, &
        local_communicator )

      call timing_nas_mg_z_mpi_pack_ddt( DIM1, DIM2, DIM3, outer_loop, &
        inner_loop, correct_flag, typesize, testname, filehandle_debug, &
        local_communicator )
!> not necessarily correct, since it assumes that a double uses twice
!> the bytes a real does
      nelements = (DIM1-2)*(DIM2-2)*2
      loops = outer_loop * inner_loop
      call time_ping_pong_nelements( nelements, loops, testname, &
        local_communicator )

      end subroutine wrapper_timing_nas_mg_z

      subroutine wrapper_timing_nas_mg( DIM1, DIM2, DIM3, outer_loop, &
        inner_loop, filehandle_correctness, filehandle_debug, &
        ptestname, local_communicator )

      implicit none

      integer, intent(in) :: DIM1, DIM2, DIM3
      integer, intent(in) :: outer_loop, inner_loop

      integer, intent(in) :: filehandle_correctness, filehandle_debug

      character(50), intent(in) :: ptestname(3)

      integer, intent(in) :: local_communicator

      call wrapper_timing_nas_mg_x( DIM1, DIM2, DIM3, outer_loop, &
        inner_loop, filehandle_correctness, filehandle_debug, &
        ptestname(1), local_communicator )

      call wrapper_timing_nas_mg_y( DIM1, DIM2, DIM3, outer_loop, &
        inner_loop, filehandle_correctness, filehandle_debug, &
        ptestname(2), local_communicator )

      call wrapper_timing_nas_mg_z( DIM1, DIM2, DIM3, outer_loop, &
        inner_loop, filehandle_correctness, filehandle_debug, &
        ptestname(3), local_communicator )

      end subroutine wrapper_timing_nas_mg

      subroutine wrapper_timing_fft( DIM1, outer_loop, inner_loop, &
        filehandle_correctness, filehandle_debug, ptestname, &
        local_communicator )

      use timing_basic, only: time_alltoall_nelements
      use timing_fft2d, only: timing_fft2d_ddt, timing_fft2d_manual, &
        timing_fft2d_mpi_pack_ddt

      implicit none

      integer, intent(in) :: DIM1

      integer, intent(in) :: outer_loop, inner_loop

      integer, intent(in) :: filehandle_correctness, filehandle_debug

      character(50), intent(in) :: ptestname

      integer, intent(in) :: local_communicator

      !> local variables
      logical :: correct_flag
      integer :: typesize
      character(50) :: testname
      integer :: nelements, loops
      integer :: myrank, procs, ier, new_dim

!> just some statements to prevent compiler warnings of unused variables
!> those parameter are included for future features
      typesize = filehandle_correctness

      call MPI_Comm_rank( local_communicator, myrank, ier )
      call MPI_Comm_size( local_communicator, procs, ier )

      if ( len(trim(ptestname)) .EQ. 0 ) then
        write (testname,'(A,I0)') "FFT", procs
      else
        testname = ptestname
      endif

      new_dim = DIM1 + mod(DIM1,procs)

      call timing_fft2d_ddt( new_dim, procs, outer_loop, inner_loop, &
        correct_flag, typesize, testname, filehandle_debug, &
        local_communicator )

      call timing_fft2d_manual( new_dim, procs, outer_loop, inner_loop, &
        correct_flag, typesize, testname, filehandle_debug, &
        local_communicator )

      call timing_fft2d_mpi_pack_ddt( new_dim, procs, outer_loop, &
        inner_loop, correct_flag, typesize, testname, &
        filehandle_debug, local_communicator )

!> not necessarily correct, since it assumes that a double complex uses forth times
!> the bytes a real does
      nelements = new_dim/procs*new_dim/procs*4
      loops = outer_loop * inner_loop
      call time_alltoall_nelements( nelements, procs, loops, testname, &
        local_communicator )

      end subroutine wrapper_timing_fft
      
      subroutine wrapper_timing_specfem3d_mt( DIM1, DIM2, DIM3, &
        outer_loop, inner_loop, filehandle_correctness, &
        filehandle_debug, ptestname, local_communicator )

      use timing_basic, only: time_ping_pong_nelements
      use timing_specfem3D, only: timing_specfem3d_mt_ddt, &
        timing_specfem3d_mt_manual, timing_specfem3d_mt_mpi_pack_ddt

      implicit none

      integer, intent(in) :: DIM1, DIM2, DIM3

      integer, intent(in) :: outer_loop, inner_loop

      integer, intent(in) :: filehandle_correctness, filehandle_debug

      character(50), intent(in) :: ptestname

      integer, intent(in) :: local_communicator

!> local variables
      logical :: correct_flag
      integer :: typesize
      character(50) :: testname
      integer :: nelements, loops

!> just some statements to prevent compiler warnings of unused variables
!> those parameter are included for future features
      typesize = filehandle_correctness

      if ( len(trim(ptestname)) .EQ. 0 ) then
        write (testname,'(A)') "SPECFEM3D_mt"
      else
        testname = ptestname
      endif

      call timing_specfem3D_mt_ddt( DIM1, DIM2, DIM3, outer_loop, &
        inner_loop, correct_flag, typesize, testname, &
        filehandle_debug, local_communicator )

      call timing_specfem3D_mt_manual( DIM1, DIM2, DIM3, outer_loop, &
        inner_loop, correct_flag, typesize, testname, &
        filehandle_debug, local_communicator )

      call timing_specfem3D_mt_mpi_pack_ddt( DIM1, DIM2, DIM3, &
        outer_loop, inner_loop, correct_flag, typesize, testname, &
        filehandle_debug, local_communicator )

      nelements = DIM1 * DIM3 
      loops = outer_loop * inner_loop
      call time_ping_pong_nelements( nelements, loops, testname, &
        local_communicator )

      end subroutine wrapper_timing_specfem3d_mt

      subroutine wrapper_timing_lammps_full( DIM1, icount, outer_loop,&
        inner_loop, filehandle_correctness, filehandle_debug, &
        ptestname, local_communicator )

      use timing_basic, only: time_ping_pong_nelements
      use timing_lammps, only: timing_lammps_full_ddt, &
        timing_lammps_full_manual, timing_lammps_full_mpi_pack_ddt
      use utilities, only: random_array_shuffle

      implicit none

      integer, intent(in) :: DIM1
      integer, intent(in) :: icount

      integer, intent(in) :: outer_loop, inner_loop

      integer, intent(in) :: filehandle_correctness, filehandle_debug

      character(50), intent(in) :: ptestname

      integer, intent(in) :: local_communicator

!> local variables
      logical :: correct_flag
      integer :: typesize
      character(50) :: testname
      integer :: nelements, loops, i
      integer, dimension(icount,outer_loop) :: list

!> just some statements to prevent compiler warnings of unused variables
!> those parameter are included for future features
      typesize = filehandle_correctness

      if ( len(trim(ptestname)) .EQ. 0 ) then
        write (testname,'(A)') "LAMMPS_full"
      else
        testname = ptestname
      endif

      do i=1,outer_loop
        call random_array_shuffle( list(1,i), icount, DIM1 )
      enddo

      call timing_lammps_full_ddt( DIM1, icount, list, &
        outer_loop, inner_loop, correct_flag, typesize, testname, &
        filehandle_debug, local_communicator )

      call timing_lammps_full_manual( DIM1, icount, list, &
        outer_loop, inner_loop, correct_flag, typesize, testname, &
        filehandle_debug, local_communicator )

      call timing_lammps_full_mpi_pack_ddt( DIM1, icount, list, &
        outer_loop, inner_loop, correct_flag, typesize, testname, &
        filehandle_debug, local_communicator )

!> not necessarily correct, since it assumes that a double uses twice
!> the bytes a real does
      nelements = icount * 8 * 2
      loops = outer_loop * inner_loop
      call time_ping_pong_nelements( nelements, loops, testname, &
        local_communicator )

      end subroutine wrapper_timing_lammps_full

      subroutine wrapper_timing_lammps_atomic( DIM1, icount, outer_loop,&
        inner_loop, filehandle_correctness, filehandle_debug, &
        ptestname, local_communicator )

      use timing_basic, only: time_ping_pong_nelements
      use timing_lammps, only: timing_lammps_atomic_ddt, &
        timing_lammps_atomic_manual, timing_lammps_atomic_mpi_pack_ddt
      use utilities, only: random_array_shuffle

      implicit none

      integer, intent(in) :: DIM1
      integer, intent(in) :: icount

      integer, intent(in) :: outer_loop, inner_loop

      integer, intent(in) :: filehandle_correctness, filehandle_debug

      character(50), intent(in) :: ptestname

      integer, intent(in) :: local_communicator

!> local variables
      logical :: correct_flag
      integer :: typesize
      character(50) :: testname
      integer :: nelements, loops, i
      integer, dimension(icount,outer_loop) :: list

!> just some statements to prevent compiler warnings of unused variables
!> those parameter are included for future features
      typesize = filehandle_correctness

      if ( len(trim(ptestname)) .EQ. 0 ) then
        write (testname,'(A)') "LAMMPS_atomic"
      else
        testname = ptestname
      endif

      do i=1,outer_loop
        call random_array_shuffle( list(1,i), icount, DIM1 )
      enddo

      call timing_lammps_atomic_ddt( DIM1, icount, list, &
        outer_loop, inner_loop, correct_flag, typesize, testname, &
        filehandle_debug, local_communicator )

      call timing_lammps_atomic_manual( DIM1, icount, list, &
        outer_loop, inner_loop, correct_flag, typesize, testname, &
        filehandle_debug, local_communicator )

      call timing_lammps_atomic_mpi_pack_ddt( DIM1, icount, list, &
        outer_loop, inner_loop, correct_flag, typesize, testname, &
        filehandle_debug, local_communicator )

!> not necessarily correct, since it assumes that a double uses twice
!> the bytes a real does
      nelements = icount * 6 * 2
      loops = outer_loop * inner_loop
      call time_ping_pong_nelements( nelements, loops, testname, &
        local_communicator )

      end subroutine wrapper_timing_lammps_atomic

      subroutine wrapper_timing_specfem3D_oc( DIM1, icount, outer_loop,&
        inner_loop, filehandle_correctness, filehandle_debug, &
        ptestname, local_communicator )

      use timing_basic, only: time_ping_pong_nelements
      use timing_specfem3D, only: timing_specfem3D_oc_ddt, &
        timing_specfem3D_oc_manual, timing_specfem3D_oc_mpi_pack_ddt
      use utilities, only: random_array_shuffle

      implicit none

      integer, intent(in) :: DIM1
      integer, intent(in) :: icount

      integer, intent(in) :: outer_loop, inner_loop

      integer, intent(in) :: filehandle_correctness, filehandle_debug

      character(50), intent(in) :: ptestname

      integer, intent(in) :: local_communicator

!> local variables
      logical :: correct_flag
      integer :: typesize
      character(50) :: testname
      integer :: nelements, loops, i
      integer, dimension(icount,outer_loop) :: list

!> just some statements to prevent compiler warnings of unused variables
!> those parameter are included for future features
      typesize = filehandle_correctness

      if ( len(trim(ptestname)) .EQ. 0 ) then
        write (testname,'(A)') "SPECFEM3D_oc"
      else
        testname = ptestname
      endif

      do i=1,outer_loop
        call random_array_shuffle( list(1,i), icount, DIM1 )
      enddo

      call timing_specfem3D_oc_ddt( DIM1, icount, list, &
        outer_loop, inner_loop, correct_flag, typesize, testname, &
        filehandle_debug, local_communicator )

      call timing_specfem3D_oc_manual( DIM1, icount, list, &
        outer_loop, inner_loop, correct_flag, typesize, testname, &
        filehandle_debug, local_communicator )

      call timing_specfem3D_oc_mpi_pack_ddt( DIM1, icount, list, &
        outer_loop, inner_loop, correct_flag, typesize, testname, &
        filehandle_debug, local_communicator )

      nelements = icount 
      loops = outer_loop * inner_loop
      call time_ping_pong_nelements( nelements, loops, testname, &
        local_communicator )

      end subroutine wrapper_timing_specfem3D_oc

      subroutine wrapper_timing_specfem3D_cm( DIM2_cm, DIM2_ic, &
        icount_cm, icount_ic, outer_loop, inner_loop, &
        filehandle_correctness, filehandle_debug, ptestname, &
        local_communicator )

      use timing_basic, only: time_ping_pong_nelements
      use timing_specfem3D, only: timing_specfem3D_cm_ddt, &
        timing_specfem3D_cm_manual, timing_specfem3D_cm_mpi_pack_ddt
      use utilities, only: random_array_shuffle

      implicit none

      integer, intent(in) :: DIM2_cm, DIM2_ic
      integer, intent(in) :: icount_cm, icount_ic

      integer, intent(in) :: outer_loop, inner_loop

      integer, intent(in) :: filehandle_correctness, filehandle_debug

      character(50), intent(in) :: ptestname

      integer, intent(in) :: local_communicator

!> local variables
      logical :: correct_flag
      integer :: typesize
      character(50) :: testname
      integer :: nelements, loops, i
      integer, dimension(icount_cm,outer_loop) :: list_cm
      integer, dimension(icount_ic,outer_loop) :: list_ic

!> just some statements to prevent compiler warnings of unused variables
!> those parameter are included for future features
      typesize = filehandle_correctness

      if ( len(trim(ptestname)) .EQ. 0 ) then
        write (testname,'(A)') "SPECFEM3D_cm"
      else
        testname = ptestname
      endif

      do i=1,outer_loop
        call random_array_shuffle( list_cm(1,i), icount_cm, DIM2_cm )
        call random_array_shuffle( list_ic(1,i), icount_ic, DIM2_ic )
      enddo

      call timing_specfem3D_cm_ddt( DIM2_cm, DIM2_ic, icount_cm, &
        icount_ic, list_cm, list_ic, outer_loop, inner_loop, &
        correct_flag, typesize, testname, filehandle_debug, &
        local_communicator )

      call timing_specfem3D_cm_manual( DIM2_cm, DIM2_ic, icount_cm, &
        icount_ic, list_cm, list_ic, outer_loop, inner_loop, &
        correct_flag, typesize, testname, filehandle_debug, &
        local_communicator )

      call timing_specfem3D_cm_mpi_pack_ddt( DIM2_cm, DIM2_ic, &
        icount_cm, icount_ic, list_cm, list_ic, outer_loop, &
        inner_loop, correct_flag, typesize, testname, &
        filehandle_debug, local_communicator )

      nelements = (icount_cm+icount_ic) * 3 
      loops = outer_loop * inner_loop
      call time_ping_pong_nelements( nelements, loops, testname, &
        local_communicator )

      end subroutine wrapper_timing_specfem3D_cm

      end module wrapper
