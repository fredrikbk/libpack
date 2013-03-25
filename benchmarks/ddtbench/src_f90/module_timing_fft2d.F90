      module timing_fft2d

      private

      public :: timing_fft2d_ddt
      public :: timing_fft2d_manual
      public :: timing_fft2d_mpi_pack_ddt

      contains

      subroutine timing_fft2d_ddt( DIM1, procs, outer_loop, &
        inner_loop, correct_flag, ptypesize, testname, &
        filehandle_debug, local_communicator )

!      use timing, only: timing_init, timing_record, timing_print
      use utilities, only: fill_unique_array_2D_double_complex

      implicit none

      include 'mpif.h'

      integer, intent(in) :: DIM1
      integer, intent(in) :: procs

      integer, intent(in) :: outer_loop, inner_loop

      logical, intent(out) :: correct_flag
      integer, intent(out) :: ptypesize

      character(50), intent(in) :: testname

      integer, intent(in) :: filehandle_debug

      integer, intent(in) :: local_communicator

!> local variables
      double complex, dimension(DIM1, DIM1/procs) :: matrix, recv_array

      integer :: myrank, ier
      integer :: i, j , base, typesize, bytes

!> variables for the datatype construction
      integer :: dtype_vector_t, dtype_resize_t
      integer :: dtype_scatter_t, dtype_gather_t
      integer :: t_complex_size
      integer(kind=MPI_ADDRESS_KIND) :: lb, extent

      character(50) :: method

!> just some statements to prevent compiler warnings of unused variables
!> those parameter are included for future features
      correct_flag = .false.
      ptypesize = 0
      typesize = filehandle_debug

      call MPI_Comm_rank( local_communicator, myrank, ier )

!> ================= initialize the arrays =================

      base = myrank * DIM1 * DIM1/procs * 2 + 1
      call fill_unique_array_2D_double_complex( matrix, DIM1, &
        DIM1/procs, base )
 
      if ( myrank .EQ. 0 ) then
        write (method,'(A)') "mpi_ddt"

        call MPI_Type_size( MPI_DOUBLE_COMPLEX, typesize, ier )
        bytes = DIM1/procs * DIM1 * typesize
 
        call timing_init( testname, method, bytes )
      endif
     
      do i=1,outer_loop

        call MPI_Type_vector( DIM1/procs, 1, DIM1, MPI_DOUBLE_COMPLEX, &
          dtype_vector_t, ier )
        call MPI_Type_size( MPI_DOUBLE_COMPLEX, t_complex_size, ier )
        lb = 0
        extent = t_complex_size
        call MPI_Type_create_resized( dtype_vector_t, lb, extent, &
          dtype_resize_t, ier )
        call MPI_Type_contiguous( DIM1/procs, dtype_resize_t, &
          dtype_scatter_t, ier )
        call MPI_Type_commit( dtype_scatter_t, ier )

        call MPI_Type_free( dtype_vector_t, ier )
        call MPI_Type_free( dtype_resize_t, ier )

        call MPI_Type_vector( DIM1/procs, DIM1/procs, DIM1, &
          MPI_DOUBLE_COMPLEX, dtype_vector_t, ier )
        lb = 0
        extent = DIM1/procs * t_complex_size
        call MPI_Type_create_resized( dtype_vector_t, lb, &
          extent, dtype_gather_t, ier )
        call MPI_Type_commit( dtype_gather_t, ier )

        call MPI_Type_free( dtype_vector_t, ier )

        if ( myrank .EQ. 0 ) then
          call timing_record(1)
        endif

        do j=1,inner_loop

          call MPI_Alltoall( matrix(1,1), 1, dtype_gather_t, &
            recv_array(1,1), 1, dtype_scatter_t, local_communicator, &
            ier )
          call MPI_Alltoall( recv_array(1,1), 1, dtype_gather_t, &
            matrix(1,1), 1, dtype_scatter_t, local_communicator, &
            ier )
          if ( myrank .EQ. 0 ) then
            call timing_record(3)
          endif

        enddo

        call MPI_Type_free( dtype_gather_t, ier )
        call MPI_Type_free( dtype_scatter_t, ier )

        if ( myrank .EQ. 0 ) then
          call timing_record(5)
        endif

      enddo !> outer loop

      if ( myrank .EQ. 0 ) then
        call timing_print( .true. )
      endif

      end subroutine timing_fft2d_ddt

      subroutine timing_fft2d_manual( DIM1, procs, outer_loop, &
        inner_loop, correct_flag, ptypesize, testname, &
        filehandle_debug, local_communicator )

!      use timing, only: timing_init, timing_record, timing_print
      use utilities, only: fill_unique_array_2D_double_complex

      implicit none

      include 'mpif.h'

      integer, intent(in) :: DIM1
      integer, intent(in) :: procs
      
      integer, intent(in) :: outer_loop, inner_loop

      logical, intent(out) :: correct_flag
      integer, intent(out) :: ptypesize

      character(50), intent(in) :: testname

      integer, intent(in) :: filehandle_debug

      integer, intent(in) :: local_communicator

!> local variables
      double complex, dimension(DIM1,DIM1/procs) :: matrix
      double complex, dimension(DIM1*DIM1/procs) :: recv_buffer
      double complex, dimension(:), allocatable :: buffer

      integer :: myrank, ier
      integer :: i, j, k, l, base, base2, typesize, bytes

      character(50) :: method

!> just some statements to prevent compiler warnings of unused variables
!> those parameter are included for future features
      correct_flag = .false.
      ptypesize = 0
      typesize = filehandle_debug

      call MPI_Comm_rank( local_communicator, myrank, ier )

!> ================= initialize the arrays =================

      base = myrank * DIM1 * DIM1/procs * 2 + 1
      call fill_unique_array_2D_double_complex( matrix, DIM1, &
        DIM1/procs, base )

      if ( myrank .EQ. 0 ) then
        write (method,'(A)') "manual"

        call MPI_Type_size( MPI_DOUBLE_COMPLEX, typesize, ier )
        bytes = DIM1/procs * DIM1 * typesize

        call timing_init( testname, method, bytes )
      endif
     
      do i=1,outer_loop

        allocate( buffer(DIM1*DIM1/procs), stat=ier )

        if ( myrank .EQ. 0 ) then
          call timing_record(1)
        endif

        do j=1,inner_loop
!> pack the data
          do k=0,procs-1
            do l=1,DIM1/procs
              base = k * DIM1/procs * DIM1/procs + (l-1) * DIM1/procs + 1
              base2 = k * DIM1/procs + 1
              buffer(base:base+DIM1/procs-1) = &
                matrix(base2:base2+DIM1/procs-1,l)
            enddo
          enddo

          if ( myrank .EQ. 0 ) then
            call timing_record(2)            
          endif

          call MPI_Alltoall( buffer(1), DIM1/procs*DIM1/procs, &
            MPI_DOUBLE_COMPLEX, recv_buffer(1), DIM1/procs*DIM1/procs, &
            MPI_DOUBLE_COMPLEX, local_communicator, ier )

!> unpack the data
          do k=1,DIM1/procs
            do l=1,DIM1
              base = (k-1) + (l-1) * DIM1/procs + 1
              matrix(l,k) = recv_buffer(base)
            enddo
          enddo

!> pack the data            
          do k=0,procs-1
            do l=1,DIM1/procs
              base = k * DIM1/procs * DIM1/procs + (l-1) * DIM1/procs + 1
              base2 = k * DIM1/procs + 1
              buffer(base:base+DIM1/procs-1) = &
                matrix(base2:base2+DIM1/procs-1,l)
            enddo
          enddo

          call MPI_Alltoall( buffer(1), DIM1/procs*DIM1/procs, &
            MPI_DOUBLE_COMPLEX, recv_buffer(1), DIM1/procs*DIM1/procs, &
            MPI_DOUBLE_COMPLEX, local_communicator, ier )

          if ( myrank .EQ. 0 ) then
            call timing_record(3)
          endif

!> unpack the data
          do k=1,DIM1/procs
            do l=1,DIM1
              base = (k-1) + (l-1) * DIM1/procs + 1
              matrix(l,k) = recv_buffer(base)
            enddo
          enddo          

          if ( myrank .EQ. 0 ) then
            call timing_record(4)
          endif

        enddo !> inner loop

        deallocate( buffer, stat=ier )

        if ( myrank .EQ. 0 ) then
          call timing_record(5)
        endif

      enddo !> outer loop

      if ( myrank .EQ. 0 ) then
        call timing_print( .true. )
      endif

      end subroutine timing_fft2d_manual

      subroutine timing_fft2d_mpi_pack_ddt( DIM1, procs, outer_loop, &
        inner_loop, correct_flag, ptypesize, testname, &
        filehandle_debug, local_communicator )

!      use timing, only: timing_init, timing_record, timing_print
      use utilities, only: fill_unique_array_2D_double_complex

      implicit none

      include 'mpif.h'

      integer, intent(in) :: DIM1
      integer, intent(in) :: procs

      integer, intent(in) :: outer_loop, inner_loop

      logical, intent(out) :: correct_flag
      integer, intent(out) :: ptypesize

      character(50), intent(in) :: testname

      integer, intent(in) :: filehandle_debug

      integer, intent(in) :: local_communicator

!> local variables
      double complex, dimension(DIM1,DIM1/procs) :: matrix
      double complex, dimension(DIM1*DIM1/procs) :: recv_buffer
      double complex, dimension(:), allocatable :: buffer

      integer :: myrank, ier
      integer :: i, j, base, typesize, bytes
      integer :: pos

!> variables for the datatype construction
      integer :: dtype_vector_t, dtype_resize_t
      integer :: dtype_scatter_t, dtype_gather_t
      integer(kind=MPI_ADDRESS_KIND) :: lb, extent

      character(50) ::  method

!> just some statements to prevent compiler warnings of unused variables
!> those parameter are included for future features
      correct_flag = .false.
      ptypesize = 0
      typesize = filehandle_debug

      call MPI_Comm_rank( local_communicator, myrank, ier )

!> ================= initialize the arrays =================

      base = myrank * DIM1 * DIM1/procs * 2 + 1
      call fill_unique_array_2D_double_complex( matrix, DIM1, &
        DIM1/procs, base )
 
      if ( myrank .EQ. 0 ) then
        write (method,'(A)') "mpi_pack_ddt"

        call MPI_Type_size( MPI_DOUBLE_COMPLEX, typesize, ier )
        bytes = DIM1/procs * DIM1 * typesize

        call timing_init( testname, method, bytes )
      endif
     
      do i=1,outer_loop

        allocate( buffer(DIM1*DIM1/procs), stat=ier )
        call MPI_Type_size( MPI_DOUBLE_COMPLEX, typesize, ier )
        bytes = DIM1/procs*DIM1/procs * typesize

        call MPI_Type_vector( DIM1/procs, 1, DIM1, MPI_DOUBLE_COMPLEX, &
          dtype_vector_t, ier )
        lb = 0
        extent = typesize
        call MPI_Type_create_resized( dtype_vector_t, lb, extent, &
          dtype_resize_t, ier )
        call MPI_Type_contiguous( DIM1/procs, dtype_resize_t, &
          dtype_scatter_t, ier )
        call MPI_Type_commit( dtype_scatter_t, ier )

        call MPI_Type_free( dtype_vector_t, ier )
        call MPI_Type_free( dtype_resize_t, ier )

        call MPI_Type_vector( DIM1/procs, DIM1/procs, DIM1, &
          MPI_DOUBLE_COMPLEX, dtype_vector_t, ier )
        lb = 0
        extent = DIM1/procs * typesize
        call MPI_Type_create_resized( dtype_vector_t, lb, &
          extent, dtype_gather_t, ier )
        call MPI_Type_commit( dtype_gather_t, ier )

        call MPI_Type_free( dtype_vector_t, ier )

        if ( myrank .EQ. 0 ) then
          call timing_record(1)
        endif

        do j=1,inner_loop
!> pack the data
          pos = 0
          call MPI_Pack( matrix(1,1), 1, dtype_gather_t, &
            buffer(1), bytes, pos, local_communicator, ier )
          if ( myrank .EQ. 0 ) then
            call timing_record(2)            
          endif

          call MPI_AllToAll( buffer(1), bytes, MPI_PACKED, &
            recv_buffer(1), bytes, MPI_PACKED, &
            local_communicator, ier )

!> unpack the data
          pos = 0
          call MPI_Unpack( recv_buffer, bytes, pos, matrix(1,1), 1, &
            dtype_scatter_t, local_communicator, ier )

!> pack the data
          pos = 0
          call MPI_Pack( matrix(1,1), 1, dtype_gather_t, &
            buffer(1), bytes, pos, local_communicator, ier )

          call MPI_AllToAll( buffer(1), bytes, MPI_PACKED, &
            recv_buffer(1), bytes, MPI_PACKED, &
            local_communicator, ier )

          if ( myrank .EQ. 0 ) then
            call timing_record(3)
          endif
!> unpack the data
          pos = 0
          call MPI_Unpack( recv_buffer, bytes, pos, matrix(1,1), 1, &
            dtype_scatter_t, local_communicator, ier )
          if ( myrank .EQ. 0 ) then
            call timing_record(4)
          endif
        enddo !> inner loop

        call MPI_Type_free( dtype_gather_t, ier )
        call MPI_Type_free( dtype_scatter_t, ier )

        deallocate( buffer, stat=ier )

        if ( myrank .EQ. 0 ) then
          call timing_record(5)
        endif

      enddo !> outer loop

      if ( myrank .EQ. 0 ) then
        call timing_print( .true. )
      endif

      end subroutine timing_fft2d_mpi_pack_ddt

      end module timing_fft2d
