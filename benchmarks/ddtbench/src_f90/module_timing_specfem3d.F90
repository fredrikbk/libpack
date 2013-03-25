      module timing_specfem3D

      private

      public timing_specfem3D_oc_ddt
      public timing_specfem3D_oc_manual
      public timing_specfem3D_oc_mpi_pack_ddt
      public timing_specfem3D_cm_ddt
      public timing_specfem3D_cm_manual
      public timing_specfem3D_cm_mpi_pack_ddt
      public timing_specfem3D_mt_ddt
      public timing_specfem3D_mt_manual
      public timing_specfem3D_mt_mpi_pack_ddt

      contains

!> case: SPECFEM3D GLOBE - scalar with MPI derived datatypes

      subroutine timing_specfem3D_oc_ddt( DIM1, icount, list, &
        outer_loop, inner_loop, correct_flag, ptypesize, testname, &
        filehandle_debug, local_communicator )
      
!      use timing, only: timing_init, timing_record, timing_print
      use utilities, only: fill_unique_array_1D_real

      implicit none

      include 'mpif.h'

!> parameters for the dimensions of the arrays
      integer, intent(in) :: DIM1

!> parameters for the number of elements, that are communicated
      integer, intent(in) :: icount

      integer, intent(in) :: outer_loop, inner_loop

      integer, dimension(icount,outer_loop), intent(in) :: list

      logical, intent(out) :: correct_flag
      integer, intent(out) :: ptypesize

      character(50), intent(in) :: testname

      integer, intent(in) :: filehandle_debug
      integer, intent(in) :: local_communicator

      real, dimension(DIM1) :: array

      integer :: ier, i, j, bytes, base, typesize
      integer :: myrank

!> variables for the MPI derived datatypes
      integer, dimension(:), allocatable :: displacement
      integer :: dtype_indexed_t

      integer :: itag
      parameter (itag = 0)

      character(50) :: method

!> just some statements to prevent compiler warnings of unused variables
!> those parameter are included for future features
      correct_flag = .false.
      ptypesize = 0
      typesize = filehandle_debug

      call MPI_Comm_rank( local_communicator, myrank, ier )

!> ================= initialize the array ==================

      base = myrank * DIM1 + 1
      call fill_unique_array_1D_real( array, DIM1, base)

      if ( myrank .EQ. 0 ) then
        write (method,'(A)') "mpi_ddt"

        call MPI_Type_size( MPI_REAL, typesize, ier )
        bytes = icount * typesize

        call timing_init( testname, method, bytes )
      endif

      do i=1,outer_loop

!> ========== building the MPI derived datatype ============
        allocate( displacement(icount), stat=ier )
        displacement(1:icount) = list(1:icount,i)-1

        call MPI_Type_create_indexed_block(icount, 1, displacement, &
          MPI_REAL, dtype_indexed_t, ier )
        call MPI_Type_commit( dtype_indexed_t, ier )

        deallocate( displacement, stat=ier )

        if ( myrank .EQ. 0 ) then
          call timing_record(1)
        endif

        do j=1,inner_loop
!> =============== ping pong communication =================

          if ( myrank .EQ. 0 ) then
!> send the data from rank 0 to rank 1          
            call MPI_Send( array(1), 1, dtype_indexed_t, &
              1, itag, local_communicator, ier )
!> receive the data from rank 1 back
            call MPI_Recv( array(1), 1, dtype_indexed_t, &
              1, itag, local_communicator, mpi_status_ignore, ier )
            call timing_record(3)
!> now for rank 1
          else
            call MPI_Recv( array(1), 1, dtype_indexed_t, &
              0, itag, local_communicator, MPI_STATUS_IGNORE, ier )
!> send the received data back
            call MPI_Send( array(1), 1, dtype_indexed_t, &
              0, itag, local_communicator, ier )
          endif

        enddo !> inner loop

!> ======================= clean up ========================

        call MPI_Type_free( dtype_indexed_t, ier )

        if ( myrank .EQ. 0 ) then
          call timing_record(5)
        endif

      enddo !> outer loop

      if ( myrank .EQ. 0 ) then
        call timing_print( .true. )
      endif

      end subroutine timing_specfem3D_oc_ddt

      subroutine timing_specfem3D_oc_manual( DIM1, icount, list, &
        outer_loop, inner_loop, correct_flag, ptypesize, testname, &
        filehandle_debug, local_communicator )

!      use timing, only: timing_init, timing_record, timing_print
      use utilities, only: fill_unique_array_1D_real

      implicit none

      include 'mpif.h'

!> parameters for the dimensions of the arrays
      integer, intent(in) :: DIM1

!> parameters for the number of elements, that are communicated
      integer, intent(in) :: icount

      integer, intent(in) :: outer_loop, inner_loop

      integer, dimension(icount,outer_loop), intent(in) :: list

      logical, intent(out) :: correct_flag
      integer, intent(out) :: ptypesize

      character(50), intent(in) :: testname

      integer, intent(in) :: filehandle_debug
      integer, intent(in) :: local_communicator

      real, dimension(DIM1) :: array
      real, dimension(:), allocatable :: buffer

      integer :: ier, i, j, k, bytes, base, typesize
      integer :: myrank

      integer :: itag
      parameter (itag = 0)

      character(50) :: method

!> just some statements to prevent compiler warnings of unused variables
!> those parameter are included for future features
      correct_flag = .false.
      ptypesize = 0
      typesize = filehandle_debug

      call MPI_Comm_rank( local_communicator, myrank, ier )

!> ================= initialize the array ==================

      base = myrank * DIM1 + 1
      call fill_unique_array_1D_real( array, DIM1, base)

      if ( myrank .EQ. 0 ) then
        write (method,'(A)') "manual"

        call MPI_Type_size( MPI_REAL, typesize, ier )
        bytes = icount * typesize

        call timing_init( testname, method, bytes )
      endif

      do i=1,outer_loop

        allocate( buffer(icount), stat=ier )

        if ( myrank .EQ. 0 ) then
          call timing_record(1)
        endif

        do j=1,inner_loop
!> =============== ping pong communication =================

          if ( myrank .EQ. 0 ) then
!> pack of the data
            do k=1,icount
              buffer(k) = array(list(k,i))
            enddo
            call timing_record(2)
!> send the data from rank 0 to rank 1          
            call MPI_Send( buffer, icount, MPI_REAL, &
              1, itag, local_communicator, ier )
!> receive the data from rank 1 back
            call MPI_Recv( buffer, icount, MPI_REAL, &
              1, itag, local_communicator, mpi_status_ignore, ier )
            call timing_record(3)
!> unpack of the data
            do k=1,icount
              array(list(k,i)) = buffer(k)
            enddo
            call timing_record(4)
!> now for rank 1
          else
            call MPI_Recv( buffer, icount, MPI_REAL, &
              0, itag, local_communicator, MPI_STATUS_IGNORE, ier )
            do k=1,icount
              array(list(k,i)) = buffer(k)
            enddo
            do k=1,icount
              buffer(k) = array(list(k,i))
            enddo
!> send the received data back
            call MPI_Send( buffer, icount, MPI_REAL, &
              0, itag, local_communicator, ier )
          endif

        enddo !> inner loop

!> ======================= clean up ========================

        deallocate( buffer, stat=ier )

        if ( myrank .EQ. 0 ) then
          call timing_record(5)
        endif

      enddo !> outer loop

      if ( myrank .EQ. 0 ) then
        call timing_print( .true. )
      endif

      end subroutine timing_specfem3D_oc_manual

      subroutine timing_specfem3D_oc_mpi_pack_ddt( DIM1, icount, list, &
        outer_loop, inner_loop, correct_flag, ptypesize, testname, &
        filehandle_debug, local_communicator )

!      use timing, only: timing_init, timing_record, timing_print
      use utilities, only: fill_unique_array_1D_real

      implicit none

      include 'mpif.h'

!> parameters for the dimensions of the arrays
      integer, intent(in) :: DIM1

!> parameters for the number of elements, that are communicated
      integer, intent(in) :: icount

      integer, intent(in) :: outer_loop, inner_loop

      integer, dimension(icount,outer_loop), intent(in) :: list

      logical, intent(out) :: correct_flag
      integer, intent(out) :: ptypesize

      character(50), intent(in) :: testname

      integer, intent(in) :: filehandle_debug
      integer, intent(in) :: local_communicator

      real, dimension(DIM1) :: array
      real, dimension(:), allocatable :: buffer

      integer :: ier, i, j, bytes, base, typesize, pos
      integer :: myrank

!> variables for the MPI derived datatypes
      integer, dimension(:), allocatable :: displacement
      integer :: dtype_indexed_t

      integer :: itag
      parameter (itag = 0)
      
      character(50) :: method

!> just some statements to prevent compiler warnings of unused variables
!> those parameter are included for future features
      correct_flag = .false.
      ptypesize = 0
      typesize = filehandle_debug

      call MPI_Comm_rank( local_communicator, myrank, ier )

!> ================= initialize the array ==================

      base = myrank * DIM1 + 1
      call fill_unique_array_1D_real( array, DIM1, base)

      if ( myrank .EQ. 0 ) then
        write (method,'(A)') "mpi_pack_ddt"

        call MPI_Type_size( MPI_REAL, typesize, ier )
        bytes = icount * typesize

        call timing_init( testname, method, bytes )
      endif

      do i=1,outer_loop

        allocate( buffer(icount), stat=ier )

        call MPI_Type_size( MPI_REAL, typesize, ier )
        bytes = icount * typesize

!> ========== building the MPI derived datatype ============
        allocate( displacement(icount), stat=ier )
        displacement(1:icount) = list(1:icount,i)-1

        call MPI_Type_create_indexed_block(icount, 1, displacement, &
          MPI_REAL, dtype_indexed_t, ier )
        call MPI_Type_commit( dtype_indexed_t, ier )

        deallocate( displacement, stat=ier )

        if ( myrank .EQ. 0 ) then
          call timing_record(1)
        endif

        do j=1,inner_loop
!> =============== ping pong communication =================

          if ( myrank .EQ. 0 ) then
!> pack of the data
            pos = 0
            call MPI_Pack( array(1), 1, dtype_indexed_t, &
              buffer, bytes, pos, local_communicator, ier )
            call timing_record(2)
!> send the data from rank 0 to rank 1          
            call MPI_Send( buffer, pos, MPI_PACKED, &
              1, itag, local_communicator, ier )
!> receive the data from rank 1 back
            call MPI_Recv( buffer, bytes, MPI_PACKED, &
              1, itag, local_communicator, mpi_status_ignore, ier )
            call timing_record(3)
!> unpack of the data
            pos = 0
            call MPI_Unpack( buffer, bytes, pos, array(1), 1, &
              dtype_indexed_t, local_communicator, ier )
            call timing_record(4)
!> now for rank 1
          else
            call MPI_Recv( buffer, bytes, MPI_PACKED, &
              0, itag, local_communicator, MPI_STATUS_IGNORE, ier )
            pos = 0
            call MPI_Unpack( buffer, bytes, pos, array(1), 1, &
              dtype_indexed_t, local_communicator, ier )
            pos = 0
            call MPI_Pack( array(1), 1, dtype_indexed_t, &
              buffer, bytes, pos, local_communicator, ier )
!> send the received data back
            call MPI_Send( buffer, pos, MPI_PACKED, &
              0, itag, local_communicator, ier )
          endif

        enddo !> inner loop

!> ======================= clean up ========================

        deallocate( buffer, stat=ier )

        call MPI_Type_free( dtype_indexed_t, ier )

        if ( myrank .EQ. 0 ) then
          call timing_record(5)
        endif

      enddo !> outer loop

      if ( myrank .EQ. 0 ) then
        call timing_print( .true. )
      endif

      end subroutine timing_specfem3D_oc_mpi_pack_ddt

      subroutine timing_specfem3D_cm_ddt( DIM2_cm, DIM2_ic, &
        icount_cm, icount_ic, list_cm, list_ic, outer_loop, &
        inner_loop, correct_flag, ptypesize, testname, &
        filehandle_debug, local_communicator )

!      use timing, only: timing_init, timing_record, timing_print
      use utilities, only: fill_unique_array_2D_real

      implicit none

      include 'mpif.h'

! parameters for the dimensions of the arrays
      integer, intent(in) :: DIM2_cm, DIM2_ic

! parameters for the number of columns, that are communicated
      integer, intent(in) :: icount_cm, icount_ic

      integer, intent(in) :: outer_loop, inner_loop

      integer, dimension(icount_cm,outer_loop), intent(in) :: list_cm
      integer, dimension(icount_ic,outer_loop), intent(in) :: list_ic

      logical, intent(out) :: correct_flag
      integer, intent(out) :: ptypesize

      character(50), intent(in) :: testname

      integer, intent(in) :: filehandle_debug
      integer, intent(in) :: local_communicator

      real, dimension(3, DIM2_cm) :: array_cm
      real, dimension(3, DIM2_ic) :: array_ic

      integer :: ier, i, j, base, bytes, typesize
      integer :: myrank, maximum

      integer, parameter :: itag = 0

      character(50) :: method

! variables for the MPI derived datatypes
      integer(KIND=MPI_ADDRESS_KIND), dimension(2) :: struct_displacement
      integer, dimension(2) :: dtype_temp_t, blocklength
      integer, dimension(:), allocatable :: displacement
      integer :: dtype_indexed_t

!> just some statements to prevent compiler warnings of unused variables
!> those parameter are included for future features
      correct_flag = .false.
      ptypesize = 0
      typesize = filehandle_debug

      call MPI_Comm_rank( local_communicator, myrank, ier )

! ================= initialize the arrays =================
      base = myrank * (DIM2_cm+DIM2_ic) * 3 + 1
      call fill_unique_array_2D_real( array_cm, 3, DIM2_cm, base)
      base = base + DIM2_cm * 3
      call fill_unique_array_2D_real( array_ic, 3, DIM2_ic, base)

      if ( myrank .EQ. 0 ) then
        write (method,'(A)') "mpi_ddt"

        call MPI_Type_size( MPI_REAL, typesize, ier )
        bytes = (icount_cm + icount_ic) * 3 * typesize
 
        call timing_init( testname, method, bytes )
      endif

      do i=1,outer_loop
! ========== building the MPI derived datatype ============
        if (icount_cm .LT. icount_ic) then
          maximum = icount_ic
        else
          maximum = icount_cm
        endif

        allocate( displacement(maximum), stat=ier )

        call MPI_Get_address( array_cm(1,1), struct_displacement(1), &
          ier )
        call MPI_Get_address( array_ic(1,1), struct_displacement(2), &
          ier )

        blocklength(1:2) = 1

        displacement(1:icount_cm) = (list_cm(1:icount_cm,i)-1) * 3
        call MPI_Type_create_indexed_block( icount_cm, 3, displacement, &
          MPI_REAL, dtype_temp_t(1), ier )

        displacement(1:icount_ic) = (list_ic(1:icount_ic,i)-1) * 3
        call MPI_Type_create_indexed_block( icount_ic, 3, displacement, &
          MPI_REAL, dtype_temp_t(2), ier )

        call MPI_Type_create_struct( 2, blocklength, &
          struct_displacement, dtype_temp_t, dtype_indexed_t, ier )
        call MPI_Type_commit( dtype_indexed_t, ier )
        call MPI_Type_free( dtype_temp_t(1), ier )
        call MPI_Type_free( dtype_temp_t(2), ier )

        deallocate( displacement, stat=ier )

        if ( myrank .EQ. 0 ) then
          call timing_record(1)
        endif

! =============== ping pong communication =================

        do j=1,inner_loop
          if ( myrank .EQ. 0 ) then
! send the data from rank 0 to rank 1          
            call MPI_Send( MPI_BOTTOM, 1, dtype_indexed_t, &
              1, itag, local_communicator, ier )
! receive the data from rank 1 back
            call MPI_Recv( MPI_BOTTOM, 1, dtype_indexed_t, &
              1, itag, local_communicator, MPI_STATUS_IGNORE, ier )
            call timing_record(3)
! now for rank 1
          else
! receive the data from rank 0
            call MPI_Recv( MPI_BOTTOM, 1, dtype_indexed_t, &
              0, itag, local_communicator, MPI_STATUS_IGNORE, ier )
! send the data back to rank 0
            call MPI_Send( MPI_BOTTOM, 1, dtype_indexed_t, &
              0, itag, local_communicator, ier )
          endif
       enddo ! inner loop

! ======================= clean up ========================

        call MPI_Type_free( dtype_indexed_t, ier )

        if ( myrank .EQ. 0 ) then
          call timing_record(5)
        endif

      enddo ! outer loop

      if ( myrank .EQ. 0 ) then
        call timing_print( .true. )
      endif

      end subroutine timing_specfem3D_cm_ddt

      subroutine timing_specfem3D_cm_manual( DIM2_cm, DIM2_ic, &
        icount_cm, icount_ic, list_cm, list_ic, outer_loop, &
        inner_loop, correct_flag, ptypesize, testname, &
        filehandle_debug, local_communicator )

!      use timing, only: timing_init, timing_record, timing_print
      use utilities, only: fill_unique_array_2D_real

      implicit none

      include 'mpif.h'

! parameters for the dimensions of the arrays
      integer, intent(in) :: DIM2_cm, DIM2_ic

! parameters for the number of columns, that are communicated
      integer, intent(in) :: icount_cm, icount_ic

      integer, intent(in) :: outer_loop, inner_loop

      integer, dimension(icount_cm,outer_loop), intent(in) :: list_cm
      integer, dimension(icount_ic,outer_loop), intent(in) :: list_ic

      logical, intent(out) :: correct_flag
      integer, intent(out) :: ptypesize

      character(50), intent(in) :: testname

      integer, intent(in) :: filehandle_debug
      integer, intent(in) :: local_communicator

      real, dimension(3, DIM2_cm) :: array_cm
      real, dimension(3, DIM2_ic) :: array_ic

      real, dimension(:), allocatable :: buffer

      integer :: ier, i, j, base, bytes, typesize, counter, isize, k
      integer :: myrank 

      integer, parameter :: itag = 0

      character(50) :: method

!> just some statements to prevent compiler warnings of unused variables
!> those parameter are included for future features
      correct_flag = .false.
      ptypesize = 0
      typesize = filehandle_debug

      call MPI_Comm_rank( local_communicator, myrank, ier )

! ================= initialize the arrays =================
      base = myrank * (DIM2_cm+DIM2_ic) * 3 + 1
      call fill_unique_array_2D_real( array_cm, 3, DIM2_cm, base )
      base = base + DIM2_cm * 3
      call fill_unique_array_2D_real( array_ic, 3, DIM2_ic, base )

      if ( myrank .EQ. 0 ) then
        write (method,'(A)') "manual"

        call MPI_Type_size( MPI_REAL, typesize, ier )
        bytes = (icount_cm + icount_ic) * 3 * typesize

        call timing_init( testname, method, bytes )
      endif

      do i=1,outer_loop
        isize = (icount_cm+icount_ic) * 3
        allocate( buffer(isize), stat=ier )        
        
        if ( myrank .EQ. 0 ) then
          call timing_record(1)
        endif

! =============== ping pong communication =================

        do j=1,inner_loop
          if ( myrank .EQ. 0 ) then
! pack the buffer
            counter = 1
            do k=1,icount_cm
              buffer(counter  ) = array_cm(1,list_cm(k,i))
              buffer(counter+1) = array_cm(2,list_cm(k,i))
              buffer(counter+2) = array_cm(3,list_cm(k,i))
              counter = counter + 3
            enddo
            do k=1,icount_ic
              buffer(counter  ) = array_ic(1,list_ic(k,i))
              buffer(counter+1) = array_ic(2,list_ic(k,i))
              buffer(counter+2) = array_ic(3,list_ic(k,i))
              counter = counter + 3
            enddo
            call timing_record(2)
! send the data from rank 0 to rank 1          
            call MPI_Send( buffer, isize, MPI_REAL, &
              1, itag, local_communicator, ier )
! receive the data from rank 1 back
            call MPI_Recv( buffer, isize, MPI_REAL, &
              1, itag, local_communicator, MPI_STATUS_IGNORE, ier )
            call timing_record(3)
! unpack the data
            counter = 1
            do k=1,icount_cm
              array_cm(1,list_cm(k,i)) = buffer(counter)
              array_cm(2,list_cm(k,i)) = buffer(counter+1)
              array_cm(3,list_cm(k,i)) = buffer(counter+2)
              counter = counter + 3
            enddo
            do k=1,icount_ic
              array_ic(1,list_ic(k,i)) = buffer(counter)
              array_ic(2,list_ic(k,i)) = buffer(counter+1)
              array_ic(3,list_ic(k,i)) = buffer(counter+2)
              counter = counter + 3
            enddo
            call timing_record(4)
! now for rank 1
          else
! receive the data from rank 0
            call MPI_Recv( buffer, isize, MPI_REAL, &
              0, itag, local_communicator, MPI_STATUS_IGNORE, ier )
! unpack the data
            counter = 1
            do k=1,icount_cm
              array_cm(1,list_cm(k,i)) = buffer(counter)
              array_cm(2,list_cm(k,i)) = buffer(counter+1)
              array_cm(3,list_cm(k,i)) = buffer(counter+2)
              counter = counter + 3
            enddo
            do k=1,icount_ic
              array_ic(1,list_ic(k,i)) = buffer(counter)
              array_ic(2,list_ic(k,i)) = buffer(counter+1)
              array_ic(3,list_ic(k,i)) = buffer(counter+2)
              counter = counter + 3
            enddo
! pack the data
            counter = 1
            do k=1,icount_cm
              buffer(counter  ) = array_cm(1,list_cm(k,i))
              buffer(counter+1) = array_cm(2,list_cm(k,i))
              buffer(counter+2) = array_cm(3,list_cm(k,i))
              counter = counter + 3
            enddo
            do k=1,icount_ic
              buffer(counter  ) = array_ic(1,list_ic(k,i))
              buffer(counter+1) = array_ic(2,list_ic(k,i))
              buffer(counter+2) = array_ic(3,list_ic(k,i))
              counter = counter + 3
            enddo
! send the data back to rank 0
            call MPI_Send( buffer, isize, MPI_REAL, &
              0, itag, local_communicator, ier )
          endif
       enddo ! inner loop

! ======================= clean up ========================
        deallocate( buffer, stat=ier )

        if ( myrank .EQ. 0 ) then
          call timing_record(5)
        endif

      enddo ! outer loop

      if ( myrank .EQ. 0 ) then
        call timing_print( .true. )
      endif

      end subroutine timing_specfem3D_cm_manual

      subroutine timing_specfem3D_cm_mpi_pack_ddt( DIM2_cm, DIM2_ic, &
        icount_cm, icount_ic, list_cm, list_ic, outer_loop, &
        inner_loop, correct_flag, ptypesize, testname, &
        filehandle_debug, local_communicator )

!      use timing, only: timing_init, timing_record, timing_print
      use utilities, only: fill_unique_array_2D_real

      implicit none

      include 'mpif.h'

! parameters for the dimensions of the arrays
      integer, intent(in) :: DIM2_cm, DIM2_ic

! parameters for the number of columns, that are communicated
      integer, intent(in) :: icount_cm, icount_ic

      integer, intent(in) :: outer_loop, inner_loop

      integer, dimension(icount_cm,outer_loop), intent(in) :: list_cm
      integer, dimension(icount_ic,outer_loop), intent(in) :: list_ic

      logical, intent(out) :: correct_flag
      integer, intent(out) :: ptypesize

      character(50), intent(in) :: testname

      integer, intent(in) :: filehandle_debug
      integer, intent(in) :: local_communicator

      real, dimension(3, DIM2_cm) :: array_cm
      real, dimension(3, DIM2_ic) :: array_ic

      real, dimension(:), allocatable :: buffer

      integer :: ier, i, j, base, bytes, typesize, isize
      integer :: myrank, maximum, pos

      integer, parameter :: itag = 0

      character(50) :: method

! variables for the MPI derived datatypes
      integer(KIND=MPI_ADDRESS_KIND), dimension(2) :: &
        struct_displacement
      integer, dimension(2) :: dtype_temp_t, blocklength
      integer, dimension(:), allocatable :: displacement
      integer :: dtype_indexed_t

!> just some statements to prevent compiler warnings of unused variables
!> those parameter are included for future features
      correct_flag = .false.
      ptypesize = 0
      typesize = filehandle_debug

      call MPI_Comm_rank( local_communicator, myrank, ier )

! ================= initialize the arrays =================
      base = myrank * (DIM2_cm+DIM2_ic) * 3 + 1
      call fill_unique_array_2D_real( array_cm, 3, DIM2_cm, base )
      base = base + DIM2_cm * 3
      call fill_unique_array_2D_real( array_ic, 3, DIM2_ic, base )

      if ( myrank .EQ. 0 ) then
        write (method,'(A)') "mpi_pack_ddt"

        call MPI_Type_size( MPI_REAL, typesize, ier )
        bytes = (icount_cm + icount_ic) * 3 * typesize

        call timing_init( testname, method, bytes )
      endif

      do i=1,outer_loop
        isize = (DIM2_cm+DIM2_ic) * 3
        call MPI_Type_size( MPI_REAL, typesize, ier )
        bytes = isize * typesize

        allocate( buffer(isize), stat=ier )        

! ========== building the MPI derived datatype ============
        if (icount_cm .LT. icount_ic) then
          maximum = icount_ic
        else
          maximum = icount_cm
        endif

        allocate( displacement(maximum), stat=ier )

        call MPI_Get_address( array_cm(1,1), struct_displacement(1), &
          ier )
        call MPI_Get_address( array_ic(1,1), struct_displacement(2), &
          ier )

        blocklength(1:2) = 1

        displacement(1:icount_cm) = (list_cm(1:icount_cm,i)-1) * 3
        call MPI_Type_create_indexed_block( icount_cm, 3, &
          displacement, MPI_REAL, dtype_temp_t(1), ier )

        displacement(1:icount_ic) = (list_ic(1:icount_ic,i)-1) * 3
        call MPI_Type_create_indexed_block(icount_ic, 3, &
          displacement, MPI_REAL, dtype_temp_t(2), ier )

        call MPI_Type_create_struct( 2, blocklength, &
          struct_displacement, dtype_temp_t, dtype_indexed_t, ier )
        call MPI_Type_commit( dtype_indexed_t, ier )
        call MPI_Type_free( dtype_temp_t(1), ier )
        call MPI_Type_free( dtype_temp_t(2), ier )

        deallocate( displacement, stat=ier )
 
        if ( myrank .EQ. 0 ) then
          call timing_record(1)
        endif

! =============== ping pong communication =================

        do j=1,inner_loop
          if ( myrank .EQ. 0 ) then
! pack the buffer
            pos = 0
            call MPI_Pack( MPI_BOTTOM, 1, dtype_indexed_t, &
              buffer, bytes, pos, local_communicator, ier )
            call timing_record(2)
! send the data from rank 0 to rank 1          
            call MPI_Send( buffer, pos, MPI_PACKED, &
              1, itag, local_communicator, ier )
! receive the data from rank 1 back
            call MPI_Recv( buffer, bytes, MPI_PACKED, &
              1, itag, local_communicator, MPI_STATUS_IGNORE, ier )
            call timing_record(3)
! unpack the data
            pos = 0
            call MPI_Unpack( buffer, bytes, pos, MPI_BOTTOM, 1, &
              dtype_indexed_t, local_communicator, ier )
            call timing_record(4)
! now for rank 1
          else
! receive the data from rank 0
            call MPI_Recv( buffer, bytes, MPI_PACKED, &
              0, itag, local_communicator, MPI_STATUS_IGNORE, ier )
! unpack the data
            pos = 0
            call MPI_Unpack( buffer, bytes, pos, MPI_BOTTOM, 1, &
              dtype_indexed_t, local_communicator, ier )
! pack the data
            pos = 0
            call MPI_Pack( MPI_BOTTOM, 1, dtype_indexed_t, &
              buffer, bytes, pos, local_communicator, ier )
! send the data back to rank 0
            call MPI_Send( buffer, pos, MPI_PACKED, &
              0, itag, local_communicator, ier )
          endif
       enddo ! inner loop

! ======================= clean up ========================
        call MPI_Type_free( dtype_indexed_t, ier )

        deallocate( buffer, stat=ier )

        if ( myrank .EQ. 0 ) then
          call timing_record(5)
        endif

      enddo ! outer loop

      if ( myrank .EQ. 0 ) then
        call timing_print( .true. )
      endif

      end subroutine timing_specfem3D_cm_mpi_pack_ddt

      subroutine timing_specfem3d_mt_ddt( DIM1, DIM2, DIM3, &
        outer_loop, inner_loop, correct_flag, ptypesize, testname, &
        filehandle_debug, local_communicator )

!      use timing, only: timing_init, timing_record, timing_print
      use utilities, only: fill_unique_array_3D_real

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
      real, dimension(DIM1,DIM2,DIM3) :: send_array
      real, dimension(DIM1,DIM3) :: recv_array

      integer :: myrank, ier
      integer :: i, j
      integer, parameter :: itag = 0
      integer :: base, bytes, typesize

      integer :: dtype_temp_t, dtype_send_t, dtype_recv_t

      character(50) :: method

!> just some statements to prevent compiler warnings of unused variables
!> those parameter are included for future features
      correct_flag = .false.
      ptypesize = 0
      typesize = filehandle_debug

      call MPI_Comm_rank( local_communicator, myrank, ier )

      base = myrank * DIM1 * DIM2 * DIM3 + 1
      call fill_unique_array_3D_real( send_array, DIM1, DIM2, DIM3, &
        base )

      if ( myrank .EQ. 0 ) then
        write (method,'(A)') "mpi_ddt"

        call MPI_Type_size( MPI_REAL, typesize, ier )
        bytes = DIM1 * DIM3 * typesize

        call timing_init( testname, method, bytes )
      endif

      do i=1,outer_loop
        
        call MPI_Type_contiguous( DIM1, MPI_REAL, dtype_temp_t, ier )
        call MPI_Type_vector( DIM3, 1, DIM2, dtype_temp_t, &
        dtype_send_t, ier )
        call MPI_Type_commit( dtype_send_t, ier )
        call MPI_Type_free( dtype_temp_t, ier )

        call MPI_Type_contiguous( DIM3*DIM1, MPI_REAL, dtype_recv_t, &
          ier )
        call MPI_Type_commit( dtype_recv_t, ier )

        if ( myrank .EQ. 0 ) then
          call timing_record(1)
        endif

        do j=1,inner_loop
          
          if ( myrank .EQ. 0 ) then
            call MPI_Send( send_array(1,1,1), 1, dtype_send_t, 1, &
              itag, local_communicator, ier )
            call MPI_Recv( recv_array(1,1), 1, dtype_recv_t, 1, itag, &
              local_communicator, MPI_STATUS_IGNORE, ier )
            call timing_record(3)
          else
            call MPI_Recv( recv_array(1,1), 1, dtype_recv_t, 0, itag, &
              local_communicator, MPI_STATUS_IGNORE, ier )
            call MPI_Send( send_array(1,1,1), 1, dtype_send_t, 0, &
              itag, local_communicator, ier )
          endif
    
        enddo !> inner loop

        call MPI_Type_free( dtype_send_t, ier )
        call MPI_Type_free( dtype_recv_t, ier )

        if ( myrank .EQ. 0 ) then
          call timing_record(5)
        endif

      enddo !> outer loop

      if ( myrank .EQ. 0 ) then
        call timing_print( .true. )
      endif

      end subroutine timing_specfem3d_mt_ddt

      subroutine timing_specfem3d_mt_manual( DIM1, DIM2, DIM3, &
        outer_loop, inner_loop, correct_flag, ptypesize, testname, &
        filehandle_debug, local_communicator )

!      use timing, only: timing_init, timing_record, timing_print
      use utilities, only: fill_unique_array_3D_real

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
      real, dimension(DIM1,DIM2,DIM3) :: send_array
      real, dimension(DIM1,DIM3) :: recv_array

      real, dimension(:,:), allocatable :: buffer

      integer :: myrank, ier
      integer :: i, j
      integer, parameter :: itag = 0
      integer :: base, bytes, typesize

      character(50) :: method

!> just some statements to prevent compiler warnings of unused variables
!> those parameter are included for future features
      correct_flag = .false.
      ptypesize = 0
      typesize = filehandle_debug

      call MPI_Comm_rank( local_communicator, myrank, ier )

      base = myrank * DIM1 * DIM2 * DIM3 + 1
      call fill_unique_array_3D_real( send_array, DIM1, DIM2, DIM3, &
        base )

      if ( myrank .EQ. 0 ) then
        write (method,'(A)') "manual"

        call MPI_Type_size( MPI_REAL, typesize, ier )
        bytes = DIM1 * DIM3 * typesize

        call timing_init( testname, method, bytes )
      endif

      do i=1,outer_loop

        allocate( buffer(DIM1,DIM3))

        call MPI_Type_size( MPI_REAL, typesize, ier )
        bytes = DIM1 * DIM3 * typesize
        
        if ( myrank .EQ. 0 ) then
          call timing_record(1)
        endif

        do j=1,inner_loop
          
          if ( myrank .EQ. 0 ) then
            buffer(1:DIM1,1:DIM3) = send_array(1:DIM1,1,1:DIM3)
            call timing_record(2)
            call MPI_Send( buffer, DIM1*DIM3, MPI_REAL, 1, &
              itag, local_communicator, ier )
            call MPI_Recv( recv_array, DIM1*DIM3, MPI_REAL, 1, itag, &
              local_communicator, MPI_STATUS_IGNORE, ier )
            call timing_record(3)
          else
            call MPI_Recv( recv_array, DIM1*DIM3, MPI_REAL, 0, itag, &
              local_communicator, MPI_STATUS_IGNORE, ier )
            buffer(1:DIM1,1:DIM3) = send_array(1:DIM1,1,1:DIM3)
            call MPI_Send( buffer, DIM1*DIM3, MPI_REAL, 0, &
              itag, local_communicator, ier )
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

      end subroutine timing_specfem3d_mt_manual

      subroutine timing_specfem3d_mt_mpi_pack_ddt( DIM1, DIM2, DIM3, &
        outer_loop, inner_loop, correct_flag, ptypesize, testname, &
        filehandle_debug, local_communicator )

!      use timing, only: timing_init, timing_record, timing_print
      use utilities, only: fill_unique_array_3D_real

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
      real, dimension(DIM1,DIM2,DIM3) :: send_array
      real, dimension(DIM1,DIM3) :: recv_array

      real, dimension(:), allocatable :: buffer

      integer :: myrank, ier
      integer :: i, j
      integer, parameter :: itag = 0
      integer :: base, bytes, typesize, pos

      integer :: dtype_temp_t, dtype_send_t, dtype_recv_t

      character(50) :: method

!> just some statements to prevent compiler warnings of unused variables
!> those parameter are included for future features
      correct_flag = .false.
      ptypesize = 0
      typesize = filehandle_debug

      call MPI_Comm_rank( local_communicator, myrank, ier )

      base = myrank * DIM1 * DIM2 * DIM3 + 1
      call fill_unique_array_3D_real( send_array, DIM1, DIM2, DIM3, &
        base )

      if ( myrank .EQ. 0 ) then
        write (method,'(A)') "mpi_pack_ddt"

        call MPI_Type_size( MPI_REAL, typesize, ier )
        bytes = DIM1 * DIM3 * typesize

        call timing_init( testname, method, bytes )
      endif

      do i=1,outer_loop

        allocate( buffer(DIM1*DIM3), stat=ier )
        call MPI_Type_size( MPI_REAL, typesize, ier )
        bytes = DIM1 * DIM3 * typesize
        
        call MPI_Type_contiguous( DIM1, MPI_REAL, dtype_temp_t, ier )
        call MPI_Type_vector( DIM3, 1, DIM2, dtype_temp_t, &
        dtype_send_t, ier )
        call MPI_Type_commit( dtype_send_t, ier )
        call MPI_Type_free( dtype_temp_t, ier )

        call MPI_Type_contiguous( DIM3*DIM1, MPI_REAL, dtype_recv_t, &
          ier )
        call MPI_Type_commit( dtype_recv_t, ier )

        if ( myrank .EQ. 0 ) then
          call timing_record(1)
        endif

        do j=1,inner_loop
          
          if ( myrank .EQ. 0 ) then
            pos = 0
            call MPI_Pack( send_array(1,1,1), 1, dtype_send_t, buffer, &
              bytes, pos, local_communicator, ier )
            call timing_record(2)
            call MPI_Send( buffer, pos, MPI_PACKED, 1, &
              itag, local_communicator, ier )
            call MPI_Recv( buffer, bytes, MPI_PACKED, 1, itag, &
              local_communicator, MPI_STATUS_IGNORE, ier )
            call timing_record(3)
            pos = 0
            call MPI_Unpack( buffer, bytes, pos, recv_array(1,1), 1, &
              dtype_recv_t, local_communicator, ier )
            call timing_record(4)
          else
            call MPI_Recv( buffer, bytes, MPI_PACKED, 0, itag, &
              local_communicator, MPI_STATUS_IGNORE, ier )
            pos = 0
            call MPI_Unpack( buffer, bytes, pos, recv_array(1,1), 1, &
              dtype_recv_t, local_communicator, ier )
            pos = 0
            call MPI_Pack( send_array(1,1,1), 1, dtype_send_t, buffer, &
              bytes, pos, local_communicator, ier )
            call MPI_Send( buffer, pos, MPI_PACKED, 0, itag, &
              local_communicator, ier )
          endif
    
        enddo !> inner loop

        deallocate( buffer, stat=ier )

        call MPI_Type_free( dtype_send_t, ier )
        call MPI_Type_free( dtype_recv_t, ier )

        if ( myrank .EQ. 0 ) then
          call timing_record(5)
        endif

      enddo !> outer loop

      if ( myrank .EQ. 0 ) then
        call timing_print( .true. )
      endif

      end subroutine timing_specfem3d_mt_mpi_pack_ddt

      end module timing_specfem3D
