      module timing_nas
!> module contains the timing benchmarks for NAS tests (MG/LU)
      
      private

      public :: timing_nas_lu_x_ddt
      public :: timing_nas_lu_x_manual
      public :: timing_nas_lu_x_mpi_pack_ddt
      public :: timing_nas_lu_y_ddt
      public :: timing_nas_lu_y_manual
      public :: timing_nas_lu_y_mpi_pack_ddt
      public :: timing_nas_mg_x_ddt
      public :: timing_nas_mg_x_manual
      public :: timing_nas_mg_x_mpi_pack_ddt
      public :: timing_nas_mg_y_ddt
      public :: timing_nas_mg_y_manual
      public :: timing_nas_mg_y_mpi_pack_ddt
      public :: timing_nas_mg_z_ddt
      public :: timing_nas_mg_z_manual
      public :: timing_nas_mg_z_mpi_pack_ddt

      contains

      subroutine timing_nas_lu_y_ddt( DIM2, DIM3, outer_loop, &
        inner_loop, correct_flag, ptypesize, testname, &
        filehandle_debug, local_communicator )

!      use timing, only: timing_init, timing_record, timing_print
      use utilities, only: fill_unique_array_3D_double

      implicit none

      include 'mpif.h'

      integer, intent(in) :: DIM2, DIM3
      integer, intent(in) :: inner_loop, outer_loop

      logical, intent(out) :: correct_flag
      integer, intent(out) :: ptypesize

      character(50), intent(in) :: testname

      integer, intent(in) :: filehandle_debug

      integer, intent(in) :: local_communicator

      double precision, dimension(5,DIM2+2,DIM3+2) :: array

      integer :: myrank, ier
      integer :: i, j, base, bytes, typesize
      integer, parameter :: itag = 0

      integer :: dtype_y_t, dtype_temp_t

      character(50) :: method

!> just some statements to prevent compiler warnings of unused variables      
!> those parameter are included for future features
      correct_flag = .false.
      ptypesize = 0
      typesize = filehandle_debug
      
      call MPI_Comm_rank( local_communicator, myrank, ier )

      base = myrank * 5 * (DIM2+2) * (DIM3+2) + 1
      call fill_unique_array_3D_double( array, 5, (DIM2+2), (DIM3+2), &
        base )

      if ( myrank .EQ. 0 ) then
        write (method,'(A)') "mpi_ddt"

        call MPI_Type_size( MPI_DOUBLE_PRECISION, typesize, ier )
        bytes = 5 * DIM3 * typesize

        call timing_init( testname, method, bytes )
      endif

      do i=1,outer_loop

        call MPI_Type_contiguous(5, MPI_DOUBLE_PRECISION, &
          dtype_temp_t, ier )

        call MPI_Type_vector( DIM3, 1, DIM2+2, dtype_temp_t, &
          dtype_y_t, ier )
        call MPI_Type_commit( dtype_y_t, ier )

        call MPI_Type_free( dtype_temp_t, ier )

        if ( myrank .EQ. 0 ) then
          call timing_record(1)
        endif

        do j=1,inner_loop

          if ( myrank .EQ. 0 ) then
            call MPI_Send( array(1,DIM2+1,2), 1, dtype_y_t, 1, itag, &
              local_communicator, ier )
            call MPI_Recv( array(1,1,2),1,dtype_y_t, 1, itag, &
              local_communicator, MPI_STATUS_IGNORE, ier )
            call timing_record(3)
          else
            call MPI_Recv( array(1,1,2),1,dtype_y_t, 0, itag, &
              local_communicator, MPI_STATUS_IGNORE, ier )
            call MPI_Send( array(1,DIM2+1,2), 1, dtype_y_t, 0, itag, &
              local_communicator, ier )
          endif

        enddo

        call MPI_Type_free( dtype_y_t, ier )

        if ( myrank .EQ. 0 ) then
          call timing_record(5)
        endif

      enddo 

      if ( myrank .EQ. 0 ) then
        call timing_print( .true. )
      endif

      end subroutine timing_nas_lu_y_ddt

      subroutine timing_nas_lu_y_manual( DIM2, DIM3, outer_loop, &
        inner_loop, correct_flag, ptypesize, testname, &
        filehandle_debug, local_communicator )

!      use timing, only: timing_init, timing_record, timing_print
      use utilities, only: fill_unique_array_3D_double

      implicit none

      include 'mpif.h'

      integer, intent(in) :: DIM2, DIM3
      integer, intent(in) :: inner_loop, outer_loop

      logical, intent(out) :: correct_flag
      integer, intent(out) :: ptypesize

      character(50), intent(in) :: testname

      integer, intent(in) :: filehandle_debug

      integer, intent(in) :: local_communicator

      double precision, dimension(5,DIM2+2,DIM3+2) :: array
      double precision, dimension(:), allocatable :: buffer

      integer :: myrank, ier
      integer :: i, j, k, l,  base, bytes, typesize
      integer, parameter :: itag = 0

      character(50) :: method

!> just some statements to prevent compiler warnings of unused variables      
!> those parameter are included for future features
      correct_flag = .false.
      ptypesize = 0
      typesize = filehandle_debug

      call MPI_Comm_rank( local_communicator, myrank, ier )

      base = myrank * 5 * (DIM2+2) * (DIM3+2) + 1
      call fill_unique_array_3D_double( array, 5, (DIM2+2), (DIM3+2), &
        base )
      
      if ( myrank .EQ. 0 ) then
        write (method,'(A)') "manual"

        call MPI_Type_size( MPI_DOUBLE_PRECISION, typesize, ier )
        bytes = 5 * DIM3 * typesize

        call timing_init( testname, method, bytes )
      endif

      do i=1,outer_loop

        allocate( buffer(5*DIM3), stat=ier )
        
        if ( myrank .EQ. 0 ) then
          call timing_record(1)
        endif

        do j=1,inner_loop
          if ( myrank .EQ. 0 ) then
!> pack the data
            base = 1
            do k=2,DIM3+1
              do l=1,5
                buffer(base) = array(l,DIM2+1,k)
                base = base + 1
              enddo
            enddo
            call timing_record(2)
            call MPI_Send( buffer, 5*DIM3, MPI_DOUBLE_PRECISION, 1, itag, &
              local_communicator, ier )
            call MPI_Recv( buffer, 5*DIM3, MPI_DOUBLE_PRECISION, 1, itag, &
              local_communicator, MPI_STATUS_IGNORE, ier )
            call timing_record(3)
!> unpack the data
            base = 1
            do k=2,DIM3+1
              do l=1,5
                array(l,1,k) = buffer(base)
                base = base + 1
              enddo
            enddo
            call timing_record(4)
          else
            call MPI_Recv( buffer, 5*DIM3, MPI_DOUBLE_PRECISION, 0, itag, &
              local_communicator, MPI_STATUS_IGNORE, ier )
!> unpack the data
            base = 1
            do k=2,DIM3+1
              do l=1,5
                array(l,1,k) = buffer(base)
                base = base + 1
              enddo
            enddo
!> pack the data
            base = 1
            do k=2,DIM3+1
              do l=1,5
                buffer(base) = array(l,DIM2+1,k)
                base = base + 1
              enddo
            enddo
            call MPI_Send( buffer, 5*DIM3, MPI_DOUBLE_PRECISION, 0, itag, &
              local_communicator, ier )
          endif

        enddo

        deallocate( buffer, stat=ier )

        if ( myrank .EQ. 0 ) then
          call timing_record(5)
        endif

      enddo 

      if ( myrank .EQ. 0 ) then
        call timing_print( .true. )
      endif

      end subroutine timing_nas_lu_y_manual

      subroutine timing_nas_lu_y_mpi_pack_ddt( DIM2, DIM3, outer_loop, &
        inner_loop, correct_flag, ptypesize, testname, &
        filehandle_debug, local_communicator )

!      use timing, only: timing_init, timing_record, timing_print
      use utilities, only: fill_unique_array_3D_double

      implicit none

      include 'mpif.h'

      integer, intent(in) :: DIM2, DIM3
      integer, intent(in) :: inner_loop, outer_loop

      logical, intent(out) :: correct_flag
      integer, intent(out) :: ptypesize

      character(50), intent(in) :: testname

      integer, intent(in) :: filehandle_debug

      integer, intent(in) :: local_communicator

      double precision, dimension(5,DIM2+2,DIM3+2) :: array
      double precision, dimension(:), allocatable :: buffer

      integer :: myrank, ier
      integer :: i, j, base, bytes, typesize, pos
      integer, parameter :: itag = 0

      integer :: dtype_y_t, dtype_temp_t

      character(50) :: method

!> just some statements to prevent compiler warnings of unused variables      
!> those parameter are included for future features
      correct_flag = .false.
      ptypesize = 0
      typesize = filehandle_debug
     
      call MPI_Comm_rank( local_communicator, myrank, ier )

      base = myrank * 5 * (DIM2+2) * (DIM3+2) + 1
      call fill_unique_array_3D_double( array, 5, (DIM2+2), (DIM3+2), base)

      if ( myrank .EQ. 0 ) then
        write (method,'(A)') "mpi_pack_ddt"

        call MPI_Type_size( MPI_DOUBLE_PRECISION, typesize, ier )
        bytes = 5 * DIM3 * typesize
         
        call timing_init( testname, method, bytes )
      endif

      do i=1,outer_loop

        allocate( buffer(5*DIM3), stat=ier )
        call MPI_Type_size( MPI_DOUBLE_PRECISION, typesize, ier )
        bytes = 5 * DIM3 * typesize

        call MPI_Type_contiguous(5, MPI_DOUBLE_PRECISION, &
          dtype_temp_t, ier )

        call MPI_Type_vector( DIM3, 1, DIM2+2, dtype_temp_t, &
          dtype_y_t, ier )
        call MPI_Type_commit( dtype_y_t, ier )

        call MPI_Type_free( dtype_temp_t, ier )
       
        if ( myrank .EQ. 0 ) then
          call timing_record(1)
        endif

        do j=1,inner_loop
          if ( myrank .EQ. 0 ) then
!> pack the data
            pos = 0
            call MPI_Pack( array(1,DIM2+1,2), 1, dtype_y_t, &
              buffer(1), bytes, pos, local_communicator, ier )
            call timing_record(2)
            call MPI_Send( buffer, pos, MPI_PACKED, 1, itag, &
              local_communicator, ier )
            call MPI_Recv( buffer, bytes, MPI_PACKED, 1, itag, &
              local_communicator, MPI_STATUS_IGNORE, ier )
            call timing_record(3)
!> unpack the data
            pos = 0
            call MPI_Unpack( buffer, bytes, pos, array(1,1,2), 1, &
              dtype_y_t, local_communicator, ier )
            call timing_record(4)
          else
            call MPI_Recv( buffer, bytes, MPI_PACKED, 0, itag, &
              local_communicator, MPI_STATUS_IGNORE, ier )
!> unpack the data
            pos = 0
            call MPI_Unpack( buffer, bytes, pos, array(1,1,2), 1, &
              dtype_y_t, local_communicator, ier )
!> pack the data
            pos = 0
            call MPI_Pack( array(1,DIM2+1,2), 1, dtype_y_t, &
              buffer(1), bytes, pos, local_communicator, ier )
            call MPI_Send( buffer, pos, MPI_PACKED, 0, itag, &
              local_communicator, ier )
          endif

        enddo

        deallocate( buffer, stat=ier )

        call MPI_Type_free( dtype_y_t, ier )

        if ( myrank .EQ. 0 ) then
          call timing_record(5)
        endif

      enddo 

      if ( myrank .EQ. 0 ) then
        call timing_print( .true. )
      endif

      end subroutine timing_nas_lu_y_mpi_pack_ddt

      subroutine timing_nas_lu_x_ddt( DIM2, DIM3, outer_loop, &
        inner_loop, correct_flag, ptypesize, testname, &
        filehandle_debug, local_communicator )

!      use timing, only: timing_init, timing_record, timing_print
      use utilities, only: fill_unique_array_3D_double

      implicit none

      include 'mpif.h'

      integer, intent(in) :: DIM2, DIM3
      integer, intent(in) :: inner_loop, outer_loop

      logical, intent(out) :: correct_flag
      integer, intent(out) :: ptypesize

      character(50), intent(in) :: testname

      integer, intent(in) :: filehandle_debug

      integer, intent(in) :: local_communicator

      double precision, dimension(5,DIM2+2,DIM3+2) :: array

      integer :: myrank, ier
      integer :: i, j, base, bytes, typesize
      integer, parameter :: itag = 0

      integer :: dtype_x_t

      character(50) :: method

!> just some statements to prevent compiler warnings of unused variables      
!> those parameter are included for future features
      correct_flag = .false.
      ptypesize = 0
      typesize = filehandle_debug

      call MPI_Comm_rank( local_communicator, myrank, ier )

      base = myrank * 5 * (DIM2+2) * (DIM3+2) + 1
      call fill_unique_array_3D_double( array, 5, (DIM2+2), (DIM3+2), base )

      if ( myrank .EQ. 0 ) then
        write (method,'(A)') "mpi_ddt"

        call MPI_Type_size( MPI_DOUBLE_PRECISION, typesize, ier )
        bytes = 5 * DIM2 * typesize

        call timing_init( testname, method, bytes )
      endif

      do i=1,outer_loop

        call MPI_Type_contiguous(5*DIM2, MPI_DOUBLE_PRECISION, &
          dtype_x_t, ier )

        call MPI_Type_commit( dtype_x_t, ier )

        if ( myrank .EQ. 0 ) then
          call timing_record(1)
        endif

        do j=1,inner_loop

          if ( myrank .EQ. 0 ) then
            call MPI_Send( array(1,2,DIM3+1), 1, dtype_x_t, 1, itag, &
              local_communicator, ier )
            call MPI_Recv( array(1,2,1), 1, dtype_x_t, 1, itag, &
              local_communicator, MPI_STATUS_IGNORE, ier )
            call timing_record(3)
          else
            call MPI_Recv( array(1,2,1), 1, dtype_x_t, 0, itag, &
              local_communicator, MPI_STATUS_IGNORE, ier )
            call MPI_Send( array(1,2,DIM3+1), 1, dtype_x_t, 0, itag, &
              local_communicator, ier )
          endif

        enddo

        call MPI_Type_free( dtype_x_t, ier )

        if ( myrank .EQ. 0 ) then
          call timing_record(5)
        endif

      enddo 

      if ( myrank .EQ. 0 ) then
        call timing_print( .true. )
      endif

      end subroutine timing_nas_lu_x_ddt

      subroutine timing_nas_lu_x_manual( DIM2, DIM3, outer_loop, &
        inner_loop, correct_flag, ptypesize, testname, &
        filehandle_debug, local_communicator )

!      use timing, only: timing_init, timing_record, timing_print
      use utilities, only: fill_unique_array_3D_double

      implicit none

      include 'mpif.h'

      integer, intent(in) :: DIM2, DIM3
      integer, intent(in) :: inner_loop, outer_loop

      logical, intent(out) :: correct_flag
      integer, intent(out) :: ptypesize

      character(50), intent(in) :: testname

      integer, intent(in) :: filehandle_debug

      integer, intent(in) :: local_communicator

      double precision, dimension(5,DIM2+2,DIM3+2) :: array
      double precision, dimension(:), allocatable :: buffer

      integer :: myrank, ier
      integer :: i, j, k, l,  base, bytes, typesize
      integer, parameter :: itag = 0

      character(50) :: method

!> just some statements to prevent compiler warnings of unused variables      
!> those parameter are included for future features
      correct_flag = .false.
      ptypesize = 0
      typesize = filehandle_debug

      call MPI_Comm_rank( local_communicator, myrank, ier )

      base = myrank * 5 * (DIM2+2) * (DIM3+2) + 1
      call fill_unique_array_3D_double( array, 5, (DIM2+2), (DIM3+2), base )

      if ( myrank .EQ. 0 ) then
        write (method,'(A)') "manual"          

        call MPI_Type_size( MPI_DOUBLE_PRECISION, typesize, ier )
        bytes = 5 * DIM2 * typesize

        call timing_init( testname, method, bytes )
      endif

      do i=1,outer_loop

        allocate( buffer(5*DIM2), stat=ier )
        
        if ( myrank .EQ. 0 ) then
          call timing_record(1)
        endif

        do j=1,inner_loop
          if ( myrank .EQ. 0 ) then
!> pack the data
            base = 1
            do k=2,DIM2+1
              do l=1,5
                buffer(base) = array(l,k,DIM3+1)
                base = base + 1
              enddo
            enddo
            call timing_record(2)
            call MPI_Send( buffer, 5*DIM2, MPI_DOUBLE_PRECISION, 1, itag, &
              local_communicator, ier )
            call MPI_Recv( buffer, 5*DIM2, MPI_DOUBLE_PRECISION, 1, itag, &
              local_communicator, MPI_STATUS_IGNORE, ier )
            call timing_record(3)
!> unpack the data
            base = 1
            do k=2,DIM2+1
              do l=1,5
                array(l,k,1) = buffer(base)
                base = base + 1
              enddo
            enddo
            call timing_record(4)
          else
            call MPI_Recv( buffer, 5*DIM2, MPI_DOUBLE_PRECISION, 0, itag, &
              local_communicator, MPI_STATUS_IGNORE, ier )
!> unpack the data
            base = 1
            do k=2,DIM2+1
              do l=1,5
                array(l,k,1) = buffer(base)
                base = base + 1
              enddo
            enddo
!> pack the data
            base = 1
            do k=2,DIM2+1
              do l=1,5
                buffer(base) = array(l,k,DIM3+1)
                base = base + 1
              enddo
            enddo
            call MPI_Send( buffer, 5*DIM2, MPI_DOUBLE_PRECISION, 0, itag, &
              local_communicator, ier )
          endif

        enddo

        deallocate( buffer, stat=ier )

        if ( myrank .EQ. 0 ) then
          call timing_record(5)
        endif

      enddo 

      if ( myrank .EQ. 0 ) then
        call timing_print( .true. )
      endif

      end subroutine timing_nas_lu_x_manual

      subroutine timing_nas_lu_x_mpi_pack_ddt( DIM2, DIM3, outer_loop, &
        inner_loop, correct_flag, ptypesize, testname, &
        filehandle_debug, local_communicator )

!      use timing, only: timing_init, timing_record, timing_print
      use utilities, only: fill_unique_array_3D_double

      implicit none

      include 'mpif.h'

      integer, intent(in) :: DIM2, DIM3
      integer, intent(in) :: inner_loop, outer_loop

      logical, intent(out) :: correct_flag
      integer, intent(out) :: ptypesize

      character(50), intent(in) :: testname

      integer, intent(in) :: filehandle_debug

      integer, intent(in) :: local_communicator

      double precision, dimension(5,DIM2+2,DIM3+2) :: array
      double precision, dimension(:), allocatable :: buffer

      integer :: myrank, ier
      integer :: i, j, base, bytes, typesize, pos
      integer, parameter :: itag = 0

      integer :: dtype_x_t

      character(50) :: method

!> just some statements to prevent compiler warnings of unused variables      
!> those parameter are included for future features
      correct_flag = .false.
      ptypesize = 0
      typesize = filehandle_debug
       
      call MPI_Comm_rank( local_communicator, myrank, ier )

      base = myrank * 5 * (DIM2+2) * (DIM3+2) + 1
      call fill_unique_array_3D_double( array, 5, (DIM2+2), (DIM3+2), base )

      if ( myrank .EQ. 0 ) then
        write (method,'(A)') "mpi_pack_ddt"
        
        call MPI_Type_size( MPI_DOUBLE_PRECISION, typesize, ier )
        bytes = 5 * DIM2 * typesize

        call timing_init( testname, method, bytes )
      endif

      do i=1,outer_loop

        allocate( buffer(5*DIM2), stat=ier )
        call MPI_Type_size( MPI_DOUBLE_PRECISION, typesize, ier )
        bytes = 5 * DIM2 * typesize

        call MPI_Type_contiguous(5*DIM2, MPI_DOUBLE_PRECISION, &
          dtype_x_t, ier )
        call MPI_Type_commit( dtype_x_t, ier )

        if ( myrank .EQ. 0 ) then
          call timing_record(1)
        endif

        do j=1,inner_loop
          if ( myrank .EQ. 0 ) then
!> pack the data
            pos = 0
            call MPI_Pack( array(1,2,DIM3+1), 1, dtype_x_t, &
              buffer(1), bytes, pos, local_communicator, ier )
            call timing_record(2)
            call MPI_Send( buffer, pos, MPI_PACKED, 1, itag, &
              local_communicator, ier )
            call MPI_Recv( buffer, bytes, MPI_PACKED, 1, itag, &
              local_communicator, MPI_STATUS_IGNORE, ier )
            call timing_record(3)
!> unpack the data
            pos = 0
            call MPI_Unpack( buffer, bytes, pos, array(1,2,1), 1, &
              dtype_x_t, local_communicator, ier )
            call timing_record(4)
          else
            call MPI_Recv( buffer, bytes, MPI_PACKED, 0, itag, &
              local_communicator, MPI_STATUS_IGNORE, ier )
!> unpack the data
            pos = 0
            call MPI_Unpack( buffer, bytes, pos, array(1,2,1), 1, &
              dtype_x_t, local_communicator, ier )
!> pack the data
            pos = 0
            call MPI_Pack( array(1,2,DIM3+1), 1, dtype_x_t, &
              buffer(1), bytes, pos, local_communicator, ier )
            call MPI_Send( buffer, pos, MPI_PACKED, 0, itag, &
              local_communicator, ier )
          endif

        enddo

        deallocate( buffer, stat=ier )

        call MPI_Type_free( dtype_x_t, ier )

        if ( myrank .EQ. 0 ) then
          call timing_record(5)
        endif

      enddo 

      if ( myrank .EQ. 0 ) then
        call timing_print( .true. )
      endif

      end subroutine timing_nas_lu_x_mpi_pack_ddt

      subroutine timing_nas_mg_x_ddt( DIM1, DIM2, DIM3, outer_loop, &
        inner_loop, correct_flag, ptypesize, testname, filehandle_debug, &
        local_communicator )

!      use timing, only: timing_init, timing_record, timing_print
      use utilities, only: fill_unique_array_3D_double

      implicit none

      include 'mpif.h'

      integer, intent(in) :: DIM1, DIM2, DIM3
      integer, intent(in) :: outer_loop, inner_loop

      logical, intent(out) :: correct_flag
      integer, intent(out) :: ptypesize

      character(50), intent(in) :: testname

      integer, intent(in) :: filehandle_debug

      integer, intent(in) :: local_communicator

!> local variables
      double precision, dimension(DIM1,DIM2,DIM3) :: array

      integer :: myrank, ier
      integer :: base, i, j, typesize, bytes

      integer(kind=MPI_ADDRESS_KIND) :: stride
      integer :: dtype_face_x_t, dtype_temp_t

      integer, parameter :: itag = 0

      character(50) :: method

!> just some statements to prevent compiler warnings of unused variables      
!> those parameter are included for future features
      correct_flag = .false.
      ptypesize = 0
      typesize = filehandle_debug

      call MPI_Comm_rank( local_communicator, myrank, ier )

      base = myrank * DIM1 * DIM2 * DIM3 + 1
      call fill_unique_array_3D_double( array, DIM1, DIM2, DIM3, base)

      if ( myrank .EQ. 0 ) then
         write (method,'(A)') "mpi_ddt"
  
         call MPI_Type_size( MPI_DOUBLE_PRECISION, typesize, ier )
         bytes = (DIM2-2)*(DIM3-2) * typesize

         call timing_init( testname, method, bytes )
      endif

      do i=1,outer_loop
  
        call MPI_Type_vector( DIM2-2, 1, DIM1, MPI_DOUBLE_PRECISION, &
          dtype_temp_t, ier )

        call MPI_Type_size( MPI_DOUBLE_PRECISION, typesize, ier )
        stride = DIM1*DIM2*typesize

        call MPI_Type_create_hvector( DIM3-2, 1, stride, &
          dtype_temp_t, dtype_face_x_t, ier )
        call MPI_Type_commit( dtype_face_x_t, ier )

        call MPI_Type_free( dtype_temp_t, ier )

        if ( myrank .EQ. 0 ) then
          call timing_record(1)
        endif

        do j=1,inner_loop

          if ( myrank .EQ. 0 ) then
            call MPI_Send( array(DIM1-1,2,2), 1, dtype_face_x_t, &
              1, itag, local_communicator, ier )
            call MPI_Recv( array(DIM1,2,2), 1, dtype_face_x_t, &
              1, itag, local_communicator, MPI_STATUS_IGNORE, ier )
            call timing_record(3)
          else
            call MPI_Recv( array(DIM1,2,2), 1, dtype_face_x_t, &
              0, itag, local_communicator, MPI_STATUS_IGNORE, ier )
            call MPI_Send( array(DIM1-1,2,2), 1, dtype_face_x_t, &
              0, itag, local_communicator, ier )
          endif

        enddo

        call MPI_Type_free( dtype_face_x_t, ier )

        if ( myrank .EQ. 0 ) then
          call timing_record(5)
        endif

      enddo

      if ( myrank .EQ. 0 ) then
        call timing_print( .true. )
      endif

      end subroutine timing_nas_mg_x_ddt

      subroutine timing_nas_mg_x_manual( DIM1, DIM2, DIM3, outer_loop, &
        inner_loop, correct_flag, ptypesize, testname, filehandle_debug, &
        local_communicator )

!      use timing, only: timing_init, timing_record, timing_print
      use utilities, only: fill_unique_array_3D_double

      implicit none

      include 'mpif.h'

      integer, intent(in) :: DIM1, DIM2, DIM3
      integer, intent(in) :: outer_loop, inner_loop

      logical, intent(out) :: correct_flag
      integer, intent(out) :: ptypesize

      character(50), intent(in) :: testname

      integer, intent(in) :: filehandle_debug

      integer, intent(in) :: local_communicator

!> local variables
      double precision, dimension(DIM1,DIM2,DIM3) :: array
      double precision, dimension(:), allocatable :: buffer

      integer :: myrank, ier
      integer :: base, i, j, k, l, typesize, bytes, psize

      integer, parameter :: itag = 0

      character(50) :: method

!> just some statements to prevent compiler warnings of unused variables      
!> those parameter are included for future features
      correct_flag = .false.
      ptypesize = 0
      typesize = filehandle_debug

      call MPI_Comm_rank( local_communicator, myrank, ier )

      base = myrank * DIM1 * DIM2 * DIM3 + 1
      call fill_unique_array_3D_double( array, DIM1, DIM2, DIM3, base)

      if ( myrank .EQ. 0 ) then
         write (method,'(A)') "manual"

         call MPI_Type_size( MPI_DOUBLE_PRECISION, typesize, ier )
         bytes = (DIM2-2)*(DIM3-2) * typesize

         call timing_init( testname, method, bytes )
      endif

      do i=1,outer_loop

        psize = (DIM2-2)*(DIM3-2)
        allocate( buffer(psize), stat=ier )

        if ( myrank .EQ. 0 ) then
          call timing_record(1)
        endif

        do j=1,inner_loop

          if ( myrank .EQ. 0 ) then
            base =1
            do k=2,DIM3-1
              do l=2,DIM2-1
                buffer(base) = array(DIM1-1,l,k)
                base = base + 1
              enddo
            enddo
            call timing_record(2)
            call MPI_Send( buffer, psize, MPI_DOUBLE_PRECISION, &
              1, itag, local_communicator, ier )
            call MPI_Recv( buffer, psize, MPI_DOUBLE_PRECISION, &
              1, itag, local_communicator, MPI_STATUS_IGNORE, ier )
            call timing_record(3)
            base = 1
            do k=2,DIM3-1
              do l=2,DIM2-1
                array(DIM1,l,k) = buffer(base)
                base = base + 1
              enddo
            enddo
            call timing_record(4)
          else
            call MPI_Recv( buffer, psize, MPI_DOUBLE_PRECISION, &
              0, itag, local_communicator, MPI_STATUS_IGNORE, ier )
            base = 1
            do k=2,DIM3-1
              do l=2,DIM2-1
                array(DIM1,l,k) = buffer(base)
                base = base + 1
              enddo
            enddo
            base =1
            do k=2,DIM3-1
              do l=2,DIM2-1
                buffer(base) = array(DIM1-1,l,k)
                base = base + 1
              enddo
            enddo
            call MPI_Send( buffer, psize, MPI_DOUBLE_PRECISION, &
              0, itag, local_communicator, ier )
          endif

        enddo

        deallocate( buffer, stat=ier )

        if ( myrank .EQ. 0 ) then
          call timing_record(5)
        endif

      enddo

      if ( myrank .EQ. 0 ) then
        call timing_print( .true. )
      endif

      end subroutine timing_nas_mg_x_manual
      
      subroutine timing_nas_mg_x_mpi_pack_ddt( DIM1, DIM2, DIM3, &
        outer_loop, inner_loop, correct_flag, ptypesize, testname, &
        filehandle_debug, local_communicator )

!      use timing, only: timing_init, timing_record, timing_print
      use utilities, only: fill_unique_array_3D_double

      implicit none

      include 'mpif.h'

      integer, intent(in) :: DIM1, DIM2, DIM3
      integer, intent(in) :: outer_loop, inner_loop

      logical, intent(out) :: correct_flag
      integer, intent(out) :: ptypesize

      character(50), intent(in) :: testname

      integer, intent(in) :: filehandle_debug

      integer, intent(in) :: local_communicator

!> local variables
      double precision, dimension(DIM1,DIM2,DIM3) :: array
      double precision, dimension(:), allocatable :: buffer

      integer :: myrank, ier
      integer :: base, i, j, typesize, bytes, pos

      integer(kind=MPI_ADDRESS_KIND) :: stride
      integer :: dtype_face_x_t, dtype_temp_t

      integer, parameter :: itag = 0

      character(50) :: method

!> just some statements to prevent compiler warnings of unused variables      
!> those parameter are included for future features
      correct_flag = .false.
      ptypesize = 0
      typesize = filehandle_debug

      call MPI_Comm_rank( local_communicator, myrank, ier )

      base = myrank * DIM1 * DIM2 * DIM3 + 1
      call fill_unique_array_3D_double( array, DIM1, DIM2, DIM3, base)

      if ( myrank .EQ. 0 ) then
         write (method,'(A)') "mpi_pack_ddt"

         call MPI_Type_size( MPI_DOUBLE_PRECISION, typesize, ier )
         bytes = (DIM2-2)*(DIM3-2) * typesize

         call timing_init( testname, method, bytes )
      endif

      do i=1,outer_loop

        allocate( buffer((DIM2-2)*(DIM3-2)), stat=ier )
        call MPI_Type_size( MPI_DOUBLE_PRECISION, typesize, ier )
        bytes = (DIM2-2)*(DIM3-2) * typesize
  
        call MPI_Type_vector( DIM2-2, 1, DIM1, MPI_DOUBLE_PRECISION, &
          dtype_temp_t, ier )

        call MPI_Type_size( MPI_DOUBLE_PRECISION, typesize, ier )
        stride = DIM1*DIM2*typesize

        call MPI_Type_create_hvector( DIM3-2, 1, stride, &
          dtype_temp_t, dtype_face_x_t, ier )
        call MPI_Type_commit( dtype_face_x_t, ier )

        call MPI_Type_free( dtype_temp_t, ier )

        if ( myrank .EQ. 0 ) then
          call timing_record(1)
        endif

        do j=1,inner_loop

          if ( myrank .EQ. 0 ) then
            pos = 0
            call MPI_Pack( array(DIM1-1,2,2), 1, dtype_face_x_t, &
              buffer, bytes, pos, local_communicator, ier )
            call timing_record(2)
            call MPI_Send( buffer, pos, MPI_PACKED, &
              1, itag, local_communicator, ier )
            call MPI_Recv( buffer, bytes, MPI_PACKED, &
              1, itag, local_communicator, MPI_STATUS_IGNORE, ier )
            call timing_record(3)
            pos = 0
            call MPI_Unpack( buffer, bytes, pos, array(DIM1,2,2), 1, &
              dtype_face_x_t, local_communicator, ier )
            call timing_record(4)
          else
            call MPI_Recv( buffer, bytes, MPI_PACKED, &
              0, itag, local_communicator, MPI_STATUS_IGNORE, ier )
            pos = 0
            call MPI_Unpack( buffer, bytes, pos, array(DIM1,2,2), 1, &
              dtype_face_x_t, local_communicator, ier )
            pos = 0
            call MPI_Pack( array(DIM1-1,2,2), 1, dtype_face_x_t, &
              buffer, bytes, pos, local_communicator, ier )
            call MPI_Send( buffer, pos, MPI_PACKED, &
              0, itag, local_communicator, ier )
          endif

        enddo

        deallocate( buffer, stat=ier )

        call MPI_Type_free( dtype_face_x_t, ier )

        if ( myrank .EQ. 0 ) then
          call timing_record(5)
        endif

      enddo

      if ( myrank .EQ. 0 ) then
        call timing_print( .true. )
      endif

      end subroutine timing_nas_mg_x_mpi_pack_ddt

      subroutine timing_nas_mg_y_ddt( DIM1, DIM2, DIM3, outer_loop, &
        inner_loop, correct_flag, ptypesize, testname, &
        filehandle_debug, local_communicator )

!      use timing, only: timing_init, timing_record, timing_print
      use utilities, only: fill_unique_array_3D_double

      implicit none

      include 'mpif.h'

      integer, intent(in) :: DIM1, DIM2, DIM3
      integer, intent(in) :: outer_loop, inner_loop

      logical, intent(out) :: correct_flag
      integer, intent(out) :: ptypesize

      character(50), intent(in) :: testname

      integer, intent(in) :: filehandle_debug

      integer, intent(in) :: local_communicator

!> local variables
      double precision, dimension(DIM1,DIM2,DIM3) :: array

      integer :: myrank, ier
      integer :: base, i, j, typesize, bytes

      integer :: dtype_face_y_t

      integer, parameter :: itag = 0

      character(50) :: method

!> just some statements to prevent compiler warnings of unused variables      
!> those parameter are included for future features
      correct_flag = .false.
      ptypesize = 0
      typesize = filehandle_debug

      call MPI_Comm_rank( local_communicator, myrank, ier )

      base = myrank * DIM1 * DIM2 * DIM3 + 1
      call fill_unique_array_3D_double( array, DIM1, DIM2, DIM3, base)

      if ( myrank .EQ. 0 ) then
        write (method,'(A)') "mpi_ddt"

        call MPI_Type_size( MPI_DOUBLE_PRECISION, typesize, ier )
        bytes = (DIM1-2)*(DIM3-2) * typesize

        call timing_init( testname, method, bytes )
      endif

      do i=1,outer_loop
  
        call MPI_Type_vector( DIM3-2, DIM1-2, DIM1*DIM2, &
          MPI_DOUBLE_PRECISION, dtype_face_y_t, ier )

        call MPI_Type_commit( dtype_face_y_t, ier )

        if ( myrank .EQ. 0 ) then
          call timing_record(1)
        endif

        do j=1,inner_loop

          if ( myrank .EQ. 0 ) then
            call MPI_Send( array(2,DIM2-1,2), 1, dtype_face_y_t, &
              1, itag, local_communicator, ier )
            call MPI_Recv( array(2,DIM2,2), 1, dtype_face_y_t, &
              1, itag, local_communicator, MPI_STATUS_IGNORE, ier )
            call timing_record(3)
          else
            call MPI_Recv( array(2,DIM2,2), 1, dtype_face_y_t, &
              0, itag, local_communicator, MPI_STATUS_IGNORE, ier )
            call MPI_Send( array(2,DIM2-1,2), 1, dtype_face_y_t, &
              0, itag, local_communicator, ier )
          endif

        enddo

        call MPI_Type_free( dtype_face_y_t, ier )

        if ( myrank .EQ. 0 ) then
          call timing_record(5)
        endif

      enddo

      if ( myrank .EQ. 0 ) then
        call timing_print( .true. )
      endif

      end subroutine timing_nas_mg_y_ddt

      subroutine timing_nas_mg_y_manual( DIM1, DIM2, DIM3, &
        outer_loop, inner_loop, correct_flag, ptypesize, testname, &
        filehandle_debug, local_communicator )

!      use timing, only: timing_init, timing_record, timing_print
      use utilities, only: fill_unique_array_3D_double

      implicit none

      include 'mpif.h'

      integer, intent(in) :: DIM1, DIM2, DIM3
      integer, intent(in) :: outer_loop, inner_loop

      logical, intent(out) :: correct_flag
      integer, intent(out) :: ptypesize

      character(50), intent(in) :: testname

      integer, intent(in) :: filehandle_debug

      integer, intent(in) :: local_communicator

!> local variables
      double precision, dimension(DIM1,DIM2,DIM3) :: array
      double precision, dimension(:), allocatable :: buffer

      integer :: myrank, ier
      integer :: base, i, j, k, l, typesize, bytes, psize

      integer, parameter :: itag = 0

      character(50) :: method

!> just some statements to prevent compiler warnings of unused variables      
!> those parameter are included for future features
      correct_flag = .false.
      ptypesize = 0
      typesize = filehandle_debug

      call MPI_Comm_rank( local_communicator, myrank, ier )

      base = myrank * DIM1 * DIM2 * DIM3 + 1
      call fill_unique_array_3D_double( array, DIM1, DIM2, DIM3, base)

      if ( myrank .EQ. 0 ) then
        write (method,'(A)') "manual"

        call MPI_Type_size( MPI_DOUBLE_PRECISION, typesize, ier )
        bytes = (DIM1-2)*(DIM3-2) * typesize

        call timing_init( testname, method, bytes )
      endif

      do i=1,outer_loop

        psize = (DIM1-2)*(DIM3-2)
        allocate( buffer(psize), stat=ier )

        if ( myrank .EQ. 0 ) then
          call timing_record(1)
        endif

        do j=1,inner_loop

          if ( myrank .EQ. 0 ) then
            base =1
            do k=2,DIM3-1
              do l=2,DIM1-1
                buffer(base) = array(l,DIM2-1,k)
                base = base + 1
              enddo
            enddo
            call timing_record(2)
            call MPI_Send( buffer, psize, MPI_DOUBLE_PRECISION, &
              1, itag, local_communicator, ier )
            call MPI_Recv( buffer, psize, MPI_DOUBLE_PRECISION, &
              1, itag, local_communicator, MPI_STATUS_IGNORE, ier )
            call timing_record(3)
            base = 1
            do k=2,DIM3-1
              do l=2,DIM1-1
                array(l,DIM2,k) = buffer(base)
                base = base + 1
              enddo
            enddo
            call timing_record(4)
          else
            call MPI_Recv( buffer, psize, MPI_DOUBLE_PRECISION, &
              0, itag, local_communicator, MPI_STATUS_IGNORE, ier )
            base = 1
            do k=2,DIM3-1
              do l=2,DIM1-1
                array(l,DIM2,k) = buffer(base)
                base = base + 1
              enddo
            enddo
            base =1
            do k=2,DIM3-1
              do l=2,DIM1-1
                buffer(base) = array(l,DIM2-1,k)
                base = base + 1
              enddo
            enddo
            call MPI_Send( buffer, psize, MPI_DOUBLE_PRECISION, &
              0, itag, local_communicator, ier )
          endif

        enddo

        deallocate( buffer, stat=ier )

        if ( myrank .EQ. 0 ) then
          call timing_record(5)
        endif

      enddo

      if ( myrank .EQ. 0 ) then
        call timing_print( .true. )
      endif

      end subroutine timing_nas_mg_y_manual
      
      subroutine timing_nas_mg_y_mpi_pack_ddt( DIM1, DIM2, DIM3, &
        outer_loop, inner_loop, correct_flag, ptypesize, testname, &
        filehandle_debug, local_communicator )

!      use timing, only: timing_init, timing_record, timing_print
      use utilities, only: fill_unique_array_3D_double

      implicit none

      include 'mpif.h'

      integer, intent(in) :: DIM1, DIM2, DIM3
      integer, intent(in) :: outer_loop, inner_loop

      logical, intent(out) :: correct_flag
      integer, intent(out) :: ptypesize

      character(50), intent(in) :: testname

      integer, intent(in) :: filehandle_debug

      integer, intent(in) :: local_communicator

!> local variables
      double precision, dimension(DIM1,DIM2,DIM3) :: array
      double precision, dimension(:), allocatable :: buffer

      integer :: myrank, ier
      integer :: base, i, j, typesize, bytes, pos

      integer :: dtype_face_y_t

      integer, parameter :: itag = 0

      character(50) :: method

!> just some statements to prevent compiler warnings of unused variables      
!> those parameter are included for future features
      correct_flag = .false.
      ptypesize = 0
      typesize = filehandle_debug

      call MPI_Comm_rank( local_communicator, myrank, ier )

      base = myrank * DIM1 * DIM2 * DIM3 + 1
      call fill_unique_array_3D_double( array, DIM1, DIM2, DIM3, base)

      if ( myrank .EQ. 0 ) then
        write (method,'(A)') "mpi_pack_ddt"

        call MPI_Type_size( MPI_DOUBLE_PRECISION, typesize, ier )
        bytes = (DIM1-2)*(DIM3-2) * typesize

        call timing_init( testname, method, bytes )
      endif

      do i=1,outer_loop

        allocate( buffer((DIM1-2)*(DIM3-2)), stat=ier )
        call MPI_Type_size( MPI_DOUBLE_PRECISION, typesize, ier )
        bytes = (DIM1-2)*(DIM3-2) * typesize
  
        call MPI_Type_vector( DIM3-2, DIM1-2, DIM1*DIM2, &
          MPI_DOUBLE_PRECISION, dtype_face_y_t, ier )
        call MPI_Type_commit( dtype_face_y_t, ier )

        if ( myrank .EQ. 0 ) then
          call timing_record(1)
        endif

        do j=1,inner_loop

          if ( myrank .EQ. 0 ) then
            pos = 0
            call MPI_Pack( array(2,DIM2-1,2), 1, dtype_face_y_t, &
              buffer, bytes, pos, local_communicator, ier )
            call timing_record(2)
            call MPI_Send( buffer, pos, MPI_PACKED, &
              1, itag, local_communicator, ier )
            call MPI_Recv( buffer, bytes, MPI_PACKED, &
              1, itag, local_communicator, MPI_STATUS_IGNORE, ier )
            call timing_record(3)
            pos = 0
            call MPI_Unpack( buffer, bytes, pos, array(2,DIM2,2), 1, &
              dtype_face_y_t, local_communicator, ier )
            call timing_record(4)
          else
            call MPI_Recv( buffer, bytes, MPI_PACKED, &
              0, itag, local_communicator, MPI_STATUS_IGNORE, ier )
            pos = 0
            call MPI_Unpack( buffer, bytes, pos, array(2,DIM2,2), 1, &
              dtype_face_y_t, local_communicator, ier )
            pos = 0
            call MPI_Pack( array(2,DIM2-1,2), 1, dtype_face_y_t, &
              buffer, bytes, pos, local_communicator, ier )
            call MPI_Send( buffer, pos, MPI_PACKED, &
              0, itag, local_communicator, ier )
          endif

        enddo

        deallocate( buffer, stat=ier )

        call MPI_Type_free( dtype_face_y_t, ier )

        if ( myrank .EQ. 0 ) then
          call timing_record(5)
        endif

      enddo

      if ( myrank .EQ. 0 ) then
        call timing_print( .true. )
      endif

      end subroutine timing_nas_mg_y_mpi_pack_ddt

      subroutine timing_nas_mg_z_ddt( DIM1, DIM2, DIM3, outer_loop, &
        inner_loop, correct_flag, ptypesize, testname, filehandle_debug, &
        local_communicator )

!      use timing, only: timing_init, timing_record, timing_print
      use utilities, only: fill_unique_array_3D_double

      implicit none

      include 'mpif.h'

      integer, intent(in) :: DIM1, DIM2, DIM3
      integer, intent(in) :: outer_loop, inner_loop

      logical, intent(out) :: correct_flag
      integer, intent(out) :: ptypesize

      character(50), intent(in) :: testname

      integer, intent(in) :: filehandle_debug

      integer, intent(in) :: local_communicator

!> local variables
      double precision, dimension(DIM1,DIM2,DIM3) :: array

      integer :: myrank, ier
      integer :: base, i, j, typesize, bytes

      integer :: dtype_face_z_t

      integer, parameter :: itag = 0

      character(50) :: method

!> just some statements to prevent compiler warnings of unused variables      
!> those parameter are included for future features
      correct_flag = .false.
      ptypesize = 0
      typesize = filehandle_debug

      call MPI_Comm_rank( local_communicator, myrank, ier )

      base = myrank * DIM1 * DIM2 * DIM3 + 1
      call fill_unique_array_3D_double( array, DIM1, DIM2, DIM3, base)

      if ( myrank .EQ. 0 ) then
        write (method,'(A)') "mpi_ddt"

        call MPI_Type_size( MPI_DOUBLE_PRECISION, typesize, ier )
        bytes = (DIM1-2)*(DIM2-2) * typesize

        call timing_init( testname, method, bytes )
      endif

      do i=1,outer_loop
  
        call MPI_Type_vector( DIM2-2, DIM1-2, DIM1, &
          MPI_DOUBLE_PRECISION, dtype_face_z_t, ier )
        call MPI_Type_commit( dtype_face_z_t, ier )

        if ( myrank .EQ. 0 ) then
          call timing_record(1)
        endif

        do j=1,inner_loop

          if ( myrank .EQ. 0 ) then
            call MPI_Send( array(2,2,2), 1, dtype_face_z_t, &
              1, itag, local_communicator, ier )
            call MPI_Recv( array(2,2,1), 1, dtype_face_z_t, &
              1, itag, local_communicator, MPI_STATUS_IGNORE, ier )
            call timing_record(3)
          else
            call MPI_Recv( array(2,2,1), 1, dtype_face_z_t, &
              0, itag, local_communicator, MPI_STATUS_IGNORE, ier )
            call MPI_Send( array(2,2,2), 1, dtype_face_z_t, &
              0, itag, local_communicator, ier )
          endif

        enddo

        call MPI_Type_free( dtype_face_z_t, ier )

        if ( myrank .EQ. 0 ) then
          call timing_record(5)
        endif

      enddo

      if ( myrank .EQ. 0 ) then
        call timing_print( .true. )
      endif

      end subroutine timing_nas_mg_z_ddt

      subroutine timing_nas_mg_z_manual( DIM1, DIM2, DIM3, &
        outer_loop, inner_loop, correct_flag, ptypesize, testname, &
        filehandle_debug, local_communicator )

!      use timing, only: timing_init, timing_record, timing_print
      use utilities, only: fill_unique_array_3D_double

      implicit none

      include 'mpif.h'

      integer, intent(in) :: DIM1, DIM2, DIM3
      integer, intent(in) :: outer_loop, inner_loop

      logical, intent(out) :: correct_flag
      integer, intent(out) :: ptypesize

      character(50), intent(in) :: testname

      integer, intent(in) :: filehandle_debug

      integer, intent(in) :: local_communicator

!> local variables
      double precision, dimension(DIM1,DIM2,DIM3) :: array
      double precision, dimension(:), allocatable :: buffer

      integer :: myrank, ier
      integer :: base, i, j, k, l, typesize, bytes, psize

      integer, parameter :: itag = 0

      character(50) :: method

!> just some statements to prevent compiler warnings of unused variables      
!> those parameter are included for future features
      correct_flag = .false.
      ptypesize = 0
      typesize = filehandle_debug

      call MPI_Comm_rank( local_communicator, myrank, ier )

      base = myrank * DIM1 * DIM2 * DIM3 + 1
      call fill_unique_array_3D_double( array, DIM1, DIM2, DIM3, base )

      if ( myrank .EQ. 0 ) then
        write (method,'(A)') "manual"

        call MPI_Type_size( MPI_DOUBLE_PRECISION, typesize, ier )
        bytes = (DIM1-2)*(DIM2-2) * typesize
         
        call timing_init( testname, method, bytes )
      endif

      do i=1,outer_loop

        psize = (DIM1-2)*(DIM2-2)
        allocate( buffer(psize), stat=ier )

        if ( myrank .EQ. 0 ) then
          call timing_record(1)
        endif

        do j=1,inner_loop

          if ( myrank .EQ. 0 ) then
            base =1
            do k=2,DIM2-1
              do l=2,DIM1-1
                buffer(base) = array(l,k,2)
                base = base + 1
              enddo
            enddo
            call timing_record(2)
            call MPI_Send( buffer, psize, MPI_DOUBLE_PRECISION, &
              1, itag, local_communicator, ier )
            call MPI_Recv( buffer, psize, MPI_DOUBLE_PRECISION, &
              1, itag, local_communicator, MPI_STATUS_IGNORE, ier )
            call timing_record(3)
            base = 1
            do k=2,DIM2-1
              do l=2,DIM1-1
                array(l,k,1) = buffer(base)
                base = base + 1
              enddo
            enddo
            call timing_record(4)
          else
            call MPI_Recv( buffer, psize, MPI_DOUBLE_PRECISION, &
              0, itag, local_communicator, MPI_STATUS_IGNORE, ier )
            base = 1
            do k=2,DIM2-1
              do l=2,DIM1-1
                array(l,k,1) = buffer(base)
                base = base + 1
              enddo
            enddo
            base =1
            do k=2,DIM2-1
              do l=2,DIM1-1
                buffer(base) = array(l,k,2)
                base = base + 1
              enddo
            enddo
            call MPI_Send( buffer, psize, MPI_DOUBLE_PRECISION, &
              0, itag, local_communicator, ier )
          endif

        enddo

        deallocate( buffer, stat=ier )

        if ( myrank .EQ. 0 ) then
          call timing_record(5)
        endif

      enddo

      if ( myrank .EQ. 0 ) then
        call timing_print( .true. )
      endif

      end subroutine timing_nas_mg_z_manual
      
      subroutine timing_nas_mg_z_mpi_pack_ddt( DIM1, DIM2, DIM3, &
        outer_loop, inner_loop, correct_flag, ptypesize, testname, &
        filehandle_debug, local_communicator )

!      use timing, only: timing_init, timing_record, timing_print
      use utilities, only: fill_unique_array_3D_double

      implicit none

      include 'mpif.h'

      integer, intent(in) :: DIM1, DIM2, DIM3
      integer, intent(in) :: outer_loop, inner_loop

      logical, intent(out) :: correct_flag
      integer, intent(out) :: ptypesize

      character(50), intent(in) :: testname

      integer, intent(in) :: filehandle_debug

      integer, intent(in) :: local_communicator

!> local variables
      double precision, dimension(DIM1,DIM2,DIM3) :: array
      double precision, dimension(:), allocatable :: buffer

      integer :: myrank, ier
      integer :: base, i, j, typesize, bytes, pos

      integer :: dtype_face_z_t

      integer, parameter :: itag = 0

      character(50) :: method

!> just some statements to prevent compiler warnings of unused variables      
!> those parameter are included for future features
      correct_flag = .false.
      ptypesize = 0
      typesize = filehandle_debug

      call MPI_Comm_rank( local_communicator, myrank, ier )

      base = myrank * DIM1 * DIM2 * DIM3 + 1
      call fill_unique_array_3D_double( array, DIM1, DIM2, DIM3, base )

      if ( myrank .EQ. 0 ) then
        write (method,'(A)') "mpi_pack_ddt"

        call MPI_Type_size( MPI_DOUBLE_PRECISION, typesize, ier )
        bytes = (DIM1-2)*(DIM2-2) * typesize

        call timing_init( testname, method, bytes )
      endif

      do i=1,outer_loop

        allocate( buffer((DIM1-2)*(DIM2-2)), stat=ier )
        call MPI_Type_size( MPI_DOUBLE_PRECISION, typesize, ier )
        bytes = (DIM1-2)*(DIM2-2) * typesize
  
        call MPI_Type_vector( DIM2-2, DIM1-2, DIM1, &
          MPI_DOUBLE_PRECISION, dtype_face_z_t, ier )
        call MPI_Type_commit( dtype_face_z_t, ier )

        if ( myrank .EQ. 0 ) then
          call timing_record(1)
        endif

        do j=1,inner_loop

          if ( myrank .EQ. 0 ) then
            pos = 0
            call MPI_Pack( array(2,2,2), 1, dtype_face_z_t, &
              buffer, bytes, pos, local_communicator, ier )
            call timing_record(2)
            call MPI_Send( buffer, pos, MPI_PACKED, &
              1, itag, local_communicator, ier )
            call MPI_Recv( buffer, bytes, MPI_PACKED, &
              1, itag, local_communicator, MPI_STATUS_IGNORE, ier )
            call timing_record(3)
            pos = 0
            call MPI_Unpack( buffer, bytes, pos, array(2,2,1), 1, &
              dtype_face_z_t, local_communicator, ier )
            call timing_record(4)
          else
            call MPI_Recv( buffer, bytes, MPI_PACKED, &
              0, itag, local_communicator, MPI_STATUS_IGNORE, ier )
            pos = 0
            call MPI_Unpack( buffer, bytes, pos, array(2,2,1), 1, &
              dtype_face_z_t, local_communicator, ier )
            pos = 0
            call MPI_Pack( array(2,2,2), 1, dtype_face_z_t, &
              buffer, bytes, pos, local_communicator, ier )
            call MPI_Send( buffer, pos, MPI_PACKED, &
              0, itag, local_communicator, ier )
          endif

        enddo

        deallocate( buffer, stat=ier )

        call MPI_Type_free( dtype_face_z_t, ier )

        if ( myrank .EQ. 0 ) then
          call timing_record(5)
        endif

      enddo

      if ( myrank .EQ. 0 ) then
        call timing_print( .true. )
      endif

      end subroutine timing_nas_mg_z_mpi_pack_ddt

      end module timing_nas
