!> handles all the time measurement
      module timing

      private
      double precision, dimension(0:1024), save :: timings
      integer, dimension(0:1024), save :: ids
      logical, dimension(5) :: epoch_used
      integer, save :: counter

      character(50), save :: testname, method
      character(19), dimension(5), save :: epoch_name
      integer, save :: bytes
      integer, save :: filehandle_values

      real, save :: max_tests
      real, save :: current_tests

      public :: timing_init
      public :: timing_record
      public :: timing_print
      public :: timing_open_file
      public :: timing_close_file
      public :: timing_set_max_tests

      contains

!> sets some init parameter like the testname, method and number of
!> bytes
!> also sets static the epoch names, the epoch_used array and the 
!> initial timing with a special id (zero)
      subroutine timing_init( ptestname, pmethod, pbytes )

      implicit none

      include 'mpif.h'

      character(50), intent(in) :: ptestname, pmethod
      integer, intent(in) :: pbytes

      counter = 1

      testname = ptestname
      method = pmethod
      bytes = pbytes

      epoch_name(1) = "ddt_create_overhead"
      epoch_name(2) = "pack"
      epoch_name(3) = "communication"
      epoch_name(4) = "unpack"
      epoch_name(5) = "ddt_free_overhead"

      epoch_used = .false.

      ids(0) = 0
      timings(0) = MPI_Wtime()

      end subroutine timing_init

!> records the timing for the given epoch id
!> if number of recorded timings exceed 1024, then the function empties
!> the buffer and they will be printed
      subroutine timing_record( id )

      implicit none

      include 'mpif.h'

      integer, intent(in) :: id

      timings(counter) = MPI_Wtime()
      ids(counter) = id
      counter = counter + 1
      if ( counter .GT. 1024 ) then
        call timing_print(.false.)
        ids(0) = 0
        timings(0) = MPI_Wtime()
        counter = 1
      endif

      end subroutine timing_record

!> writes the timing values to the output file
!> if the last parameter is set, then all the unused epochs are also
!> printed
      subroutine timing_print(last)

      use utilities, only: write_mpi

      implicit none

      logical, intent(in) :: last

!> local variables
      character(256) :: line
      character :: newline
      integer :: i

      newline = char(10)

      do i=1,counter-1
!> if a epoch id is not in 1..5 is found, only the id is printed (instead of
!> the the epoch name)
        if ( (ids(i) .GT. 5) .OR. (ids(i) .LT. 0) ) then
          write (line,'(A30,A30,I15,I20,F20.4,A)') &
            trim(testname), trim(method), bytes, ids(i), &
            (timings(i)-timings(i-1))*10**6, newline
        else
          write (line,'(A30,A30,I15,A20,F20.4,A)') &
            trim(testname), trim(method), bytes, &
            trim(epoch_name(ids(i))), &
            (timings(i)-timings(i-1))*10**6, newline
!> keep tracking if a epoch is used or not
          if (epoch_used(ids(i)) .EQV. .false.) then
            epoch_used(ids(i)) = .true.
          endif
        endif
        call write_mpi( line, filehandle_values)
      enddo

      if (last) then
!> each epoch must appear for all tests to make the analysis afterwards
!> easier
!> at the end of the test, we print for each unused epoch a line with a 
!> timing of zero micro seconds
        do i=1,5
          if (epoch_used(i) .EQV. .false.) then
            write (line,'(A30,A30,I15,A20,F20.4,A)') &
              trim(testname), trim(method), bytes, trim(epoch_name(i)), &
              real(0), newline
            call write_mpi( line, filehandle_values)
            current_tests = current_tests + 1
          endif
        enddo
      endif

!> prints some progression information to stdout
      current_tests = current_tests + counter - 1
      write (*,'(A,F7.3,A)') " Finished ", current_tests/max_tests*100, &
        "% of all tests"

      end subroutine timing_print

!> opens a file handle, to where the timing values are written
      subroutine timing_open_file(filename)

      use utilities, only: write_mpi

      implicit none

      include 'mpif.h'

      character(50), intent(in) :: filename

      integer :: ier
      character :: newline
      character(256) :: line

      newline = char(10)

      call MPI_File_open( MPI_COMM_SELF, filename, &
        MPI_MODE_EXCL + MPI_MODE_CREATE + MPI_MODE_WRONLY, &
        MPI_INFO_NULL, filehandle_values, ier )

      if ( ier .NE. MPI_SUCCESS ) then
        write (*,'(A,A,A,A,A,A)') "Error at open file ", trim(filename), &
          " for writing the timing values. The file probably already", &
          " exists.", newline, "Will now abort.", newline
        call MPI_Abort( MPI_COMM_WORLD, 1, ier )
      endif

      write (line,'(A30,A30,A15,A20,A20,A)') "testname", &
        "method", "bytes", "id", "time", newline
      call write_mpi( line, filehandle_values)

      end subroutine timing_open_file

!> closes the above mentioned filehandle
      subroutine timing_close_file

      implicit none
      
      include 'mpif.h'

      integer :: ier

      call MPI_File_close( filehandle_values, ier )

      end subroutine timing_close_file

!> sets the initial parameter for the progression bar
      subroutine timing_set_max_tests( value )

      implicit none

      integer, intent(in) :: value

      max_tests = value
      current_tests = 0      

      end subroutine

      end module timing
