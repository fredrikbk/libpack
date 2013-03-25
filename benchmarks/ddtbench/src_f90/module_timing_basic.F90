      module timing_basic

      private 

      public :: time_ping_pong_nelements
      public :: time_alltoall_nelements

      contains

      subroutine time_ping_pong_nelements(DIM1, loop, testname, &
        local_communicator)

!      use timing, only: timing_init, timing_record, timing_print
      use utilities, only: fill_unique_array_1D_real

      implicit none

      include 'mpif.h'

      integer, intent(in) :: DIM1
      integer, intent(in) :: loop

      character(50), intent(in) :: testname

      integer, intent(in) :: local_communicator

! local variables
      real, dimension(DIM1) :: array
      integer :: myrank, ier
      integer :: base, typesize, bytes, i
      integer, parameter :: itag = 0
      character(50) :: method

      call MPI_Comm_rank( local_communicator, myrank, ier )

      base = myrank * DIM1 + 1
      call fill_unique_array_1D_real( array, DIM1, base )

      if ( myrank .EQ. 0 ) then
        write (method,'(A)') "reference"

        call MPI_Type_size(MPI_REAL, typesize, ier )
        bytes = typesize * DIM1

        call timing_init( testname, method, bytes )
      endif

      do i=1,loop
        if ( myrank .EQ. 0 ) then
          call MPI_Send( array(1), DIM1, MPI_REAL, 1, itag, &
            local_communicator, ier )
          call MPI_Recv( array(1), DIM1, MPI_REAL, 1, itag, &
            local_communicator, MPI_STATUS_IGNORE, ier )
          call timing_record(3)
        else
          call MPI_Recv( array(1), DIM1, MPI_REAL, 0, itag, &
            local_communicator, MPI_STATUS_IGNORE, ier )
          call MPI_Send( array(1), DIM1, MPI_REAL, 0, itag, &
            local_communicator, ier )
        endif
      enddo

      if ( myrank .EQ. 0 ) then
        call timing_print(.true.)
      endif
     
      end subroutine time_ping_pong_nelements

      subroutine time_alltoall_nelements(DIM1, procs, loop, testname, &
        local_communicator)

!      use timing, only: timing_init, timing_record, timing_print
      use utilities, only: fill_unique_array_1D_real

      implicit none

      include 'mpif.h'

      integer, intent(in) :: DIM1
      integer, intent(in) :: procs
      integer, intent(in) :: loop

      character(50), intent(in) :: testname

      integer, intent(in) :: local_communicator

! local variables
      real, dimension(DIM1*procs) :: send_array, recv_array
      integer :: myrank, ier
      integer :: base, typesize, bytes, i
      integer, parameter :: itag = 0
      character(50) :: method

      call MPI_Comm_rank( local_communicator, myrank, ier )

      base = myrank * DIM1 + 1
      call fill_unique_array_1D_real( send_array, DIM1, base )

      if ( myrank .EQ. 0 ) then
        write (method,'(A)') "reference"
        
        call MPI_Type_size(MPI_REAL, typesize, ier )
        bytes = typesize * DIM1 * procs

        call timing_init( testname, method, bytes )
      endif

      do i=1,loop
        call MPI_Alltoall(send_array(1), DIM1, MPI_REAL, &
          recv_array(1), DIM1, MPI_REAL, local_communicator, ier )

        call MPI_Alltoall(recv_array(1), DIM1, MPI_REAL, &
          send_array(1), DIM1, MPI_REAL, local_communicator, ier )

        if ( myrank .EQ. 0 ) then
          call timing_record(3)
        endif
      enddo

      if ( myrank .EQ. 0 ) then
        call timing_print(.true.)
      endif
     
      end subroutine time_alltoall_nelements

      end module timing_basic
