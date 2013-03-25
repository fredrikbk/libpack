      module timing_lammps

        private 

        public :: timing_lammps_full_ddt
        public :: timing_lammps_full_manual
        public :: timing_lammps_full_mpi_pack_ddt
        public :: timing_lammps_atomic_ddt
        public :: timing_lammps_atomic_manual
        public :: timing_lammps_atomic_mpi_pack_ddt

      contains

      subroutine timing_lammps_full_ddt( DIM1, icount, list, &
        outer_loop, inner_loop, correct_flag, ptypesize, testname, &
        filehandle_debug, local_communicator)

!      use timing, only: timing_init, timing_record, timing_print
      use utilities, only: fill_unique_array_1D_double, &
        fill_unique_array_2D_double

      implicit none

      include 'mpif.h'

      integer, intent(in) :: DIM1
      integer, intent(in) :: icount

      integer, intent(in) :: outer_loop, inner_loop

      integer, intent(in), dimension(icount,outer_loop) :: list

      logical, intent(out) :: correct_flag
      integer, intent(out) :: ptypesize

      character(50), intent(in) :: testname

      integer, intent(in) :: filehandle_debug

      integer, intent(in) :: local_communicator

      double precision, dimension(DIM1+icount) :: atag, atype, amask, &
        amolecule, aq
      double precision, dimension(3,DIM1+icount) :: ax

      integer :: myrank, ier
      integer, parameter :: itag = 0
      integer :: i, j
      integer :: typesize, bytes, base

      integer :: dtype_indexed1_t, dtype_indexed3_t, &
        dtype_send_t, dtype_cont1_t, dtype_cont3_t, dtype_recv_t
      integer(kind=MPI_ADDRESS_KIND), dimension(6) :: &
        address_displacement
      integer, dimension(6) :: blocklength, oldtype
      integer, dimension(:), allocatable :: index_displacement

      character(50) :: method

!> just some statements to prevent compiler warnings of unused variables
!> those parameter are included for future features
      correct_flag = .false.
      ptypesize = 0
      typesize = filehandle_debug

      call MPI_Comm_rank( local_communicator, myrank, ier )

      base = myrank * (8*(DIM1+icount)) + 1
      call fill_unique_array_2D_double( ax, 3, DIM1+icount, base )
      base = base + 3*(DIM1+icount)
      call fill_unique_array_1D_double( atag, DIM1+icount, base )
      base = base + DIM1 + icount
      call fill_unique_array_1D_double( atype, DIM1+icount, base )
      base = base + DIM1 + icount
      call fill_unique_array_1D_double( amask, DIM1+icount, base )
      base = base + DIM1 + icount
      call fill_unique_array_1D_double( aq, DIM1+icount, base )
      base = base + DIM1 + icount
      call fill_unique_array_1D_double( amolecule, DIM1+icount, base )

      if ( myrank .EQ. 0 ) then
        write (method,'(A)') "mpi_ddt"

        call MPI_Type_size( MPI_DOUBLE_PRECISION, typesize, ier )
        bytes = icount * 8 * typesize
 
        call timing_init( testname, method, bytes )
      endif

      do i=1,outer_loop

        allocate(index_displacement(icount), stat=ier)

        index_displacement(1:icount) = list(1:icount,i) - 1
        call MPI_Type_create_indexed_block( icount, 1, index_displacement, &
          MPI_DOUBLE_PRECISION, dtype_indexed1_t, ier )
      
        index_displacement(1:icount) = 3 * index_displacement(1:icount)
        call MPI_Type_create_indexed_block( icount, 3, index_displacement, &
          MPI_DOUBLE_PRECISION, dtype_indexed3_t, ier )

        call MPI_Get_address( ax(1,1), address_displacement(1), ier )
        call MPI_Get_address( atag(1), address_displacement(2), ier )
        call MPI_Get_address( atype(1), address_displacement(3), ier )
        call MPI_Get_address( amask(1), address_displacement(4), ier )
        call MPI_Get_address( aq(1), address_displacement(5), ier )
        call MPI_Get_address( amolecule(1), address_displacement(6), ier )

        oldtype(1) = dtype_indexed3_t
        oldtype(2:6) = dtype_indexed1_t
 
        blocklength = 1
 
        call MPI_Type_create_struct( 6, blocklength, address_displacement, &
          oldtype, dtype_send_t, ier )
        call MPI_Type_commit( dtype_send_t, ier )

        call MPI_Type_contiguous( icount, MPI_DOUBLE_PRECISION, &
          dtype_cont1_t, ier );
        call MPI_Type_contiguous( 3*icount, MPI_DOUBLE_PRECISION, &
          dtype_cont3_t, ier );

        call MPI_Get_address( ax(1,DIM1+1), address_displacement(1), ier )
        call MPI_Get_address( atag(DIM1+1), address_displacement(2), ier )
        call MPI_Get_address( atype(DIM1+1), address_displacement(3), ier )
        call MPI_Get_address( amask(DIM1+1), address_displacement(4), ier )
        call MPI_Get_address( aq(DIM1+1), address_displacement(5), ier )
        call MPI_Get_address( amolecule(DIM1+1), address_displacement(6), ier )

        oldtype(1) = dtype_cont3_t
        oldtype(2:6) = dtype_cont1_t

        call MPI_Type_create_struct( 6, blocklength, address_displacement, &
          oldtype, dtype_recv_t, ier )
        call MPI_Type_commit( dtype_recv_t, ier )

        call MPI_Type_free( dtype_indexed1_t, ier )
        call MPI_Type_free( dtype_indexed3_t, ier )
        call MPI_Type_free( dtype_cont1_t, ier )
        call MPI_Type_free( dtype_cont3_t, ier )

        deallocate(index_displacement, stat=ier )

        if ( myrank .EQ. 0 ) then
          call timing_record(1)
        endif

        do j=1,inner_loop

          if ( myrank .EQ. 0 ) then
            call MPI_Send( MPI_BOTTOM, 1, dtype_send_t, 1, itag, &
              local_communicator, ier )
            call MPI_Recv( MPI_BOTTOM, 1, dtype_recv_t, 1, itag, &
              local_communicator, MPI_STATUS_IGNORE, ier )
            call timing_record(3)
          else
            call MPI_Recv( MPI_BOTTOM, 1, dtype_send_t, 0, itag, &
              local_communicator, MPI_STATUS_IGNORE, ier )
            call MPI_Send( MPI_BOTTOM, 1, dtype_recv_t, 0, itag, &
              local_communicator, ier )
          endif

        enddo ! inner loop

        call MPI_Type_free( dtype_send_t, ier )
        call MPI_Type_free( dtype_recv_t, ier )

        if ( myrank .EQ. 0 ) then
          call timing_record(5)
        endif

      enddo ! outer loop

      if ( myrank .EQ. 0 ) then
        call timing_print( .true. )
      endif

      end subroutine timing_lammps_full_ddt

      subroutine timing_lammps_full_manual( DIM1, icount, list, &
        outer_loop, inner_loop, correct_flag, ptypesize, testname, &
        filehandle_debug, local_communicator)

!      use timing, only: timing_init, timing_record, timing_print
      use utilities, only: fill_unique_array_1D_double, &
        fill_unique_array_2D_double

      implicit none

      include 'mpif.h'

      integer, intent(in) :: DIM1
      integer, intent(in) :: icount

      integer, intent(in) :: outer_loop, inner_loop

      integer, dimension(icount,outer_loop) :: list

      logical, intent(out) :: correct_flag
      integer, intent(out) :: ptypesize

      character(50), intent(in) :: testname

      integer, intent(in) :: filehandle_debug

      integer, intent(in) :: local_communicator

      double precision, dimension(DIM1+icount) :: atag, atype, amask, &
        amolecule, aq
      double precision, dimension(3,DIM1+icount) :: ax
      double precision, dimension(:), allocatable :: buffer

      integer :: myrank, ier
      integer, parameter :: itag = 0
      integer :: i, j, k, l
      integer :: typesize, bytes, base, pos, isize

      character(50) :: method

!> just some statements to prevent compiler warnings of unused variables
!> those parameter are included for future features
      correct_flag = .false.
      ptypesize = 0
      typesize = filehandle_debug

      call MPI_Comm_rank( local_communicator, myrank, ier )

      base = myrank * (8*(DIM1+icount)) + 1
      call fill_unique_array_2D_double( ax, 3, DIM1+icount, base )
      base = base + 3*(DIM1+icount)
      call fill_unique_array_1D_double( atag, DIM1+icount, base )
      base = base + DIM1 + icount
      call fill_unique_array_1D_double( atype, DIM1+icount, base )
      base = base + DIM1 + icount
      call fill_unique_array_1D_double( amask, DIM1+icount, base )
      base = base + DIM1 + icount
      call fill_unique_array_1D_double( aq, DIM1+icount, base )
      base = base + DIM1 + icount
      call fill_unique_array_1D_double( amolecule, DIM1+icount, base )
     
      if ( myrank .EQ. 0 ) then
        write (method,'(A)') "manual"

        call MPI_Type_size( MPI_DOUBLE_PRECISION, typesize, ier )
        bytes = icount * 8 * typesize
 
        call timing_init( testname, method, bytes )
      endif

      do k=1,outer_loop

        isize = 8*icount
        call MPI_Type_size( MPI_DOUBLE_PRECISION, typesize, ier )
        bytes = isize * typesize
        allocate( buffer(isize), stat=ier ) 

        if ( myrank .EQ. 0 ) then
          call timing_record(1)
        endif

        do l=1,inner_loop

          if ( myrank .EQ. 0 ) then
            pos = 1
            do i=1,icount
              j=list(i,k)
              buffer(pos) = ax(1,j)
              buffer(pos+1) = ax(2,j)
              buffer(pos+2) = ax(3,j)
              buffer(pos+3) = atag(j)
              buffer(pos+4) = atype(j)
              buffer(pos+5) = amask(j)
              buffer(pos+6) = amolecule(j)
              buffer(pos+7) = aq(j)
              pos = pos + 8
            enddo
            call timing_record(2)
            call MPI_Send( buffer, isize, MPI_DOUBLE_PRECISION, 1, itag, &
              local_communicator, ier )
            call MPI_Recv( buffer, isize, MPI_DOUBLE_PRECISION, 1, itag, &
              local_communicator, MPI_STATUS_IGNORE, ier )
            call timing_record(3)
            pos = 1
            do i=1,icount
              j=DIM1+i
              ax(1,j) = buffer(pos)
              ax(2,j) = buffer(pos+1)
              ax(3,j) = buffer(pos+2)
              atag(j) = buffer(pos+3)
              atype(j) = buffer(pos+4)
              amask(j) = buffer(pos+5)
              amolecule(j) = buffer(pos+6)
              aq(j) = buffer(pos+7)
              pos = pos + 8
            enddo
            call timing_record(4)
          else
            call MPI_Recv( buffer, isize, MPI_DOUBLE_PRECISION, 0, itag, &
              local_communicator, MPI_STATUS_IGNORE, ier )
            pos = 1
            do i=1,icount
              j=DIM1+i
              ax(1,j) = buffer(pos)
              ax(2,j) = buffer(pos+1)
              ax(3,j) = buffer(pos+2)
              atag(j) = buffer(pos+3)
              atype(j) = buffer(pos+4)
              amask(j) = buffer(pos+5)
              amolecule(j) = buffer(pos+6)
              aq(j) = buffer(pos+7)
              pos = pos + 8
            enddo
            pos = 1
            do i=1,icount
              j=list(i,k)
              buffer(pos) = ax(1,j)
              buffer(pos+1) = ax(2,j)
              buffer(pos+2) = ax(3,j)
              buffer(pos+3) = atag(j)
              buffer(pos+4) = atype(j)
              buffer(pos+5) = amask(j)
              buffer(pos+6) = amolecule(j)
              buffer(pos+7) = aq(j)
              pos = pos + 8
            enddo
            call MPI_Send( buffer, isize, MPI_DOUBLE_PRECISION, 0, itag, &
              local_communicator, ier )
          endif

        enddo ! inner loop

        deallocate( buffer, stat=ier )

        if ( myrank .EQ. 0 ) then
          call timing_record(5)
        endif

      enddo ! outer loop

      if ( myrank .EQ. 0 ) then
        call timing_print( .true. )
      endif

      end subroutine timing_lammps_full_manual

      subroutine timing_lammps_full_mpi_pack_ddt( DIM1, icount, list, &
        outer_loop, inner_loop, correct_flag, ptypesize, testname, &
        filehandle_debug, local_communicator)

!      use timing, only: timing_init, timing_record, timing_print
      use utilities, only: fill_unique_array_1D_double, &
        fill_unique_array_2D_double

      implicit none

      include 'mpif.h'

      integer, intent(in) :: DIM1
      integer, intent(in) :: icount

      integer, intent(in) :: outer_loop, inner_loop

      integer, dimension(icount,outer_loop) :: list

      logical, intent(out) :: correct_flag
      integer, intent(out) :: ptypesize

      character(50), intent(in) :: testname

      integer, intent(in) :: filehandle_debug

      integer, intent(in) :: local_communicator

      double precision, dimension(DIM1+icount) :: atag, atype, amask, &
        amolecule, aq
      double precision, dimension(3,DIM1+icount) :: ax
      double precision, dimension(:), allocatable :: buffer

      integer :: myrank, ier
      integer, parameter :: itag = 0
      integer :: i, j
      integer :: typesize, bytes, base, pos

      integer :: dtype_indexed1_t, dtype_indexed3_t, &
        dtype_send_t, dtype_cont1_t, dtype_cont3_t, dtype_recv_t
      integer(kind=MPI_ADDRESS_KIND), dimension(6) :: address_displacement
      integer, dimension(6) :: blocklength, oldtype
      integer, dimension(:), allocatable :: index_displacement

      character(50) :: method

!> just some statements to prevent compiler warnings of unused variables
!> those parameter are included for future features
      correct_flag = .false.
      ptypesize = 0
      typesize = filehandle_debug

      call MPI_Comm_rank( local_communicator, myrank, ier )

      base = myrank * (8*(DIM1+icount)) + 1
      call fill_unique_array_2D_double( ax, 3, DIM1+icount, base )
      base = base + 3*(DIM1+icount)
      call fill_unique_array_1D_double( atag, DIM1+icount, base )
      base = base + DIM1 + icount
      call fill_unique_array_1D_double( atype, DIM1+icount, base )
      base = base + DIM1 + icount
      call fill_unique_array_1D_double( amask, DIM1+icount, base )
      base = base + DIM1 + icount
      call fill_unique_array_1D_double( aq, DIM1+icount, base )
      base = base + DIM1 + icount
      call fill_unique_array_1D_double( amolecule, DIM1+icount, base )

      if ( myrank .EQ. 0 ) then
        write (method,'(A)') "mpi_pack_ddt"

        call MPI_Type_size( MPI_DOUBLE_PRECISION, typesize, ier )
        bytes = icount * 8 * typesize

        call timing_init( testname, method, bytes )
      endif

      do i=1,outer_loop

        call MPI_Type_size( MPI_DOUBLE_PRECISION, typesize, ier )
        bytes = 8 * icount * typesize 
        allocate( buffer(8*icount), stat=ier ) 

        allocate(index_displacement(icount), stat=ier)
        index_displacement(1:icount) = list(1:icount,i) - 1
        call MPI_Type_create_indexed_block( icount, 1, index_displacement, &
          MPI_DOUBLE_PRECISION, dtype_indexed1_t, ier )
      
        index_displacement(1:icount) = 3 * index_displacement(1:icount)
        call MPI_Type_create_indexed_block( icount, 3, index_displacement, &
          MPI_DOUBLE_PRECISION, dtype_indexed3_t, ier )

        call MPI_Get_address( ax(1,1), address_displacement(1), ier )
        call MPI_Get_address( atag(1), address_displacement(2), ier )
        call MPI_Get_address( atype(1), address_displacement(3), ier )
        call MPI_Get_address( amask(1), address_displacement(4), ier )
        call MPI_Get_address( amolecule(1), address_displacement(5), ier )
        call MPI_Get_address( aq(1), address_displacement(6), ier )

        oldtype(1) = dtype_indexed3_t
        oldtype(2:6) = dtype_indexed1_t
 
        blocklength = 1
 
        call MPI_Type_create_struct( 6, blocklength, address_displacement, &
          oldtype, dtype_send_t, ier )
        call MPI_Type_commit( dtype_send_t, ier )

        call MPI_Type_contiguous( icount, MPI_DOUBLE_PRECISION, &
          dtype_cont1_t, ier );
        call MPI_Type_contiguous( 3*icount, MPI_DOUBLE_PRECISION, &
          dtype_cont3_t, ier );

        call MPI_Get_address( ax(1,DIM1+1), address_displacement(1), ier )
        call MPI_Get_address( atag(DIM1+1), address_displacement(2), ier )
        call MPI_Get_address( atype(DIM1+1), address_displacement(3), ier )
        call MPI_Get_address( amask(DIM1+1), address_displacement(4), ier )
        call MPI_Get_address( aq(DIM1+1), address_displacement(5), ier )
        call MPI_Get_address( amolecule(DIM1+1), address_displacement(6), ier )

        oldtype(1) = dtype_cont3_t
        oldtype(2:6) = dtype_cont1_t

        call MPI_Type_create_struct( 6, blocklength, address_displacement, &
          oldtype, dtype_recv_t, ier )
        call MPI_Type_commit( dtype_recv_t, ier )

        call MPI_Type_free( dtype_indexed1_t, ier )
        call MPI_Type_free( dtype_indexed3_t, ier )
        call MPI_Type_free( dtype_cont1_t, ier )
        call MPI_Type_free( dtype_cont3_t, ier )

        deallocate(index_displacement, stat=ier )

        if ( myrank .EQ. 0 ) then
          call timing_record(1)
        endif

        do j=1,inner_loop

          if ( myrank .EQ. 0 ) then
            pos = 0
            call MPI_Pack( MPI_BOTTOM, 1, dtype_send_t, buffer, &
              bytes, pos, local_communicator, ier )
            call timing_record(2)
            call MPI_Send( buffer, pos, MPI_PACKED, 1, itag, &
              local_communicator, ier )
            call MPI_Recv( buffer, bytes, MPI_PACKED, 1, itag, &
              local_communicator, MPI_STATUS_IGNORE, ier )
            call timing_record(3)
            pos = 0
            call MPI_Unpack( buffer, bytes, pos, MPI_BOTTOM, 1, &
              dtype_recv_t, local_communicator, ier )
            call timing_record(4)
          else
            call MPI_Recv( buffer, bytes, MPI_PACKED, 0, itag, &
              local_communicator, MPI_STATUS_IGNORE, ier )
            pos = 0
            call MPI_Unpack( buffer, bytes, pos, MPI_BOTTOM, 1, &
              dtype_recv_t, local_communicator, ier )
            pos = 0
            call MPI_Pack( MPI_BOTTOM, 1, dtype_send_t, buffer, &
              bytes, pos, local_communicator, ier )           
            call MPI_Send( buffer, pos, MPI_PACKED, 0, itag, &
              local_communicator, ier )
          endif

        enddo ! inner loop

        deallocate( buffer, stat=ier )

        call MPI_Type_free( dtype_send_t, ier )
        call MPI_Type_free( dtype_recv_t, ier )

        if ( myrank .EQ. 0 ) then
          call timing_record(5)
        endif

      enddo ! outer loop

      if ( myrank .EQ. 0 ) then
        call timing_print( .true. )
      endif

      end subroutine timing_lammps_full_mpi_pack_ddt

      subroutine timing_lammps_atomic_ddt( DIM1, icount, list, &
        outer_loop, inner_loop, correct_flag, ptypesize, testname, &
        filehandle_debug, local_communicator)

!      use timing, only: timing_init, timing_record, timing_print
      use utilities, only: fill_unique_array_1D_double, &
        fill_unique_array_2D_double

      implicit none

      include 'mpif.h'

      integer, intent(in) :: DIM1
      integer, intent(in) :: icount

      integer, intent(in) :: outer_loop, inner_loop

      integer, intent(in), dimension(icount,outer_loop) :: list

      logical, intent(out) :: correct_flag
      integer, intent(out) :: ptypesize

      character(50), intent(in) :: testname

      integer, intent(in) :: filehandle_debug

      integer, intent(in) :: local_communicator

      double precision, dimension(DIM1+icount) :: atag, atype, amask
      double precision, dimension(3,DIM1+icount) :: ax

      integer :: myrank, ier
      integer, parameter :: itag = 0
      integer :: i, j
      integer :: typesize, bytes, base

      integer :: dtype_indexed1_t, dtype_indexed3_t, &
        dtype_send_t, dtype_cont1_t, dtype_cont3_t, dtype_recv_t
      integer(kind=MPI_ADDRESS_KIND), dimension(4) :: address_displacement
      integer, dimension(4) :: blocklength, oldtype
      integer, dimension(:), allocatable :: index_displacement

      character(50) :: method

!> just some statements to prevent compiler warnings of unused variables
!> those parameter are included for future features
      correct_flag = .false.
      ptypesize = 0
      typesize = filehandle_debug

      call MPI_Comm_rank( local_communicator, myrank, ier )

      base = myrank * (6*(DIM1+icount)) + 1
      call fill_unique_array_2D_double( ax, 3, DIM1+icount, base )
      base = base + 3*(DIM1+icount)
      call fill_unique_array_1D_double( atag, DIM1+icount, base )
      base = base + DIM1 + icount
      call fill_unique_array_1D_double( atype, DIM1+icount, base )
      base = base + DIM1 + icount
      call fill_unique_array_1D_double( amask, DIM1+icount, base )

      if ( myrank .EQ. 0 ) then
        write (method,'(A)') "mpi_ddt"

        call MPI_Type_size( MPI_DOUBLE_PRECISION, typesize, ier )
        bytes = icount * 6 * typesize

        call timing_init( testname, method, bytes )
      endif

      do i=1,outer_loop

        allocate(index_displacement(icount), stat=ier)

        index_displacement(1:icount) = list(1:icount,i) - 1
        call MPI_Type_create_indexed_block( icount, 1, index_displacement, &
          MPI_DOUBLE_PRECISION, dtype_indexed1_t, ier )
      
        index_displacement(1:icount) = 3 * index_displacement(1:icount)
        call MPI_Type_create_indexed_block( icount, 3, index_displacement, &
          MPI_DOUBLE_PRECISION, dtype_indexed3_t, ier )

        call MPI_Get_address( ax(1,1), address_displacement(1), ier )
        call MPI_Get_address( atag(1), address_displacement(2), ier )
        call MPI_Get_address( atype(1), address_displacement(3), ier )
        call MPI_Get_address( amask(1), address_displacement(4), ier )

        oldtype(1) = dtype_indexed3_t
        oldtype(2:4) = dtype_indexed1_t
 
        blocklength = 1
 
        call MPI_Type_create_struct( 4, blocklength, address_displacement, &
          oldtype, dtype_send_t, ier )
        call MPI_Type_commit( dtype_send_t, ier )

        call MPI_Type_contiguous( icount, MPI_DOUBLE_PRECISION, &
          dtype_cont1_t, ier )
        call MPI_Type_contiguous( 3*icount, MPI_DOUBLE_PRECISION, &
          dtype_cont3_t, ier )

        call MPI_Get_address( ax(1,DIM1+1), address_displacement(1), ier )
        call MPI_Get_address( atag(DIM1+1), address_displacement(2), ier )
        call MPI_Get_address( atype(DIM1+1), address_displacement(3), ier )
        call MPI_Get_address( amask(DIM1+1), address_displacement(4), ier )

        oldtype(1) = dtype_cont3_t
        oldtype(2:4) = dtype_cont1_t
 
        call MPI_Type_create_struct( 4, blocklength, address_displacement, &
          oldtype, dtype_recv_t, ier )
        call MPI_Type_commit( dtype_recv_t, ier )

        call MPI_Type_free( dtype_indexed1_t, ier )
        call MPI_Type_free( dtype_indexed3_t, ier )
        call MPI_Type_free( dtype_cont1_t, ier )
        call MPI_Type_free( dtype_cont3_t, ier )

        deallocate(index_displacement, stat=ier )

        if ( myrank .EQ. 0 ) then
          call timing_record(1)
        endif

        do j=1,inner_loop

          if ( myrank .EQ. 0 ) then
            call MPI_Send( MPI_BOTTOM, 1, dtype_send_t, 1, itag, &
              local_communicator, ier )
            call MPI_Recv( MPI_BOTTOM, 1, dtype_recv_t, 1, itag, &
              local_communicator, MPI_STATUS_IGNORE, ier )
            call timing_record(3)
          else
            call MPI_Recv( MPI_BOTTOM, 1, dtype_recv_t, 0, itag, &
              local_communicator, MPI_STATUS_IGNORE, ier )
            call MPI_Send( MPI_BOTTOM, 1, dtype_send_t, 0, itag, &
              local_communicator, ier )
          endif

        enddo ! inner loop

        call MPI_Type_free( dtype_send_t, ier )
        call MPI_Type_free( dtype_recv_t, ier )

        if ( myrank .EQ. 0 ) then
          call timing_record(5)
        endif

      enddo ! outer loop

      if ( myrank .EQ. 0 ) then
        call timing_print( .true. )
      endif

      end subroutine timing_lammps_atomic_ddt

      subroutine timing_lammps_atomic_manual( DIM1, icount, list, &
        outer_loop, inner_loop, correct_flag, ptypesize, testname, &
        filehandle_debug, local_communicator)

!      use timing, only: timing_init, timing_record, timing_print
      use utilities, only: fill_unique_array_1D_double, &
        fill_unique_array_2D_double

      implicit none

      include 'mpif.h'

      integer, intent(in) :: DIM1
      integer, intent(in) :: icount

      integer, intent(in) :: outer_loop, inner_loop

      integer, dimension(icount,outer_loop) :: list

      logical, intent(out) :: correct_flag
      integer, intent(out) :: ptypesize

      character(50), intent(in) :: testname

      integer, intent(in) :: filehandle_debug

      integer, intent(in) :: local_communicator

      double precision, dimension(DIM1+icount) :: atag, atype, amask
      double precision, dimension(3,DIM1+icount) :: ax
      double precision, dimension(:), allocatable :: buffer

      integer :: myrank, ier
      integer, parameter :: itag = 0
      integer :: i, j, k, l
      integer :: typesize, bytes, base, pos, isize

      character(50) :: method

!> just some statements to prevent compiler warnings of unused variables
!> those parameter are included for future features
      correct_flag = .false.
      ptypesize = 0
      typesize = filehandle_debug

      call MPI_Comm_rank( local_communicator, myrank, ier )

      base = myrank * (6*(DIM1+icount)) + 1
      call fill_unique_array_2D_double( ax, 3, DIM1+icount, base )
      base = base + 3*(DIM1+icount)
      call fill_unique_array_1D_double( atag, DIM1+icount, base )
      base = base + DIM1 + icount
      call fill_unique_array_1D_double( atype, DIM1+icount, base )
      base = base + DIM1 + icount
      call fill_unique_array_1D_double( amask, DIM1+icount, base )

      if ( myrank .EQ. 0 ) then
        write (method,'(A)') "manual"

        call MPI_Type_size( MPI_DOUBLE_PRECISION, typesize, ier )
        bytes = icount * 6 * typesize

        call timing_init( testname, method, bytes )
      endif

      do k=1,outer_loop

        isize = 6*icount
        call MPI_Type_size( MPI_DOUBLE_PRECISION, typesize, ier )
        bytes = isize * typesize
        allocate( buffer(isize), stat=ier ) 

        if ( myrank .EQ. 0 ) then
          call timing_record(1)
        endif

        do l=1,inner_loop

          if ( myrank .EQ. 0 ) then
            pos = 1
            do i=1,icount
              j = list(i,k)
              buffer(pos) = ax(1,j)
              buffer(pos+1) = ax(2,j)
              buffer(pos+2) = ax(3,j)
              buffer(pos+3) = atag(j)
              buffer(pos+4) = atype(j)
              buffer(pos+5) = amask(j)
              pos = pos + 6
            enddo
            call timing_record(2)
            call MPI_Send( buffer, isize, MPI_DOUBLE_PRECISION, 1, itag, &
              local_communicator, ier )
            call MPI_Recv( buffer, isize, MPI_DOUBLE_PRECISION, 1, itag, &
              local_communicator, MPI_STATUS_IGNORE, ier )
            call timing_record(3)
            pos = 1
            do i=1,icount
              j = DIM1+i
              ax(1,j) = buffer(pos)
              ax(2,j) = buffer(pos+1)
              ax(3,j) = buffer(pos+2)
              atag(j) = buffer(pos+3)
              atype(j) = buffer(pos+4)
              amask(j) = buffer(pos+5)
              pos = pos + 6
            enddo
            call timing_record(4)
          else
            call MPI_Recv( buffer, isize, MPI_DOUBLE_PRECISION, 0, itag, &
              local_communicator, MPI_STATUS_IGNORE, ier )
            pos = 1
            do i=1,icount
              j = DIM1+i
              ax(1,j) = buffer(pos)
              ax(2,j) = buffer(pos+1)
              ax(3,j) = buffer(pos+2)
              atag(j) = buffer(pos+3)
              atype(j) = buffer(pos+4)
              amask(j) = buffer(pos+5)
              pos = pos + 6
            enddo
            pos = 1
            do i=1,icount
              j = list(i,k)
              buffer(pos) = ax(1,j)
              buffer(pos+1) = ax(2,j)
              buffer(pos+2) = ax(3,j)
              buffer(pos+3) = atag(j)
              buffer(pos+4) = atype(j)
              buffer(pos+5) = amask(j)
              pos = pos + 6
            enddo
            call MPI_Send( buffer, isize, MPI_DOUBLE_PRECISION, 0, itag, &
              local_communicator, ier )
          endif

        enddo ! inner loop

        deallocate( buffer, stat=ier )

        if ( myrank .EQ. 0 ) then
          call timing_record(5)
        endif

      enddo ! outer loop

      if ( myrank .EQ. 0 ) then
        call timing_print( .true. )
      endif

      end subroutine timing_lammps_atomic_manual

      subroutine timing_lammps_atomic_mpi_pack_ddt( DIM1, icount, list,&
        outer_loop, inner_loop, correct_flag, ptypesize, testname, &
        filehandle_debug, local_communicator )

!      use timing, only: timing_init, timing_record, timing_print
      use utilities, only: fill_unique_array_1D_double, &
        fill_unique_array_2D_double

      implicit none

      include 'mpif.h'

      integer, intent(in) :: DIM1
      integer, intent(in) :: icount

      integer, intent(in) :: outer_loop, inner_loop

      integer, dimension(icount,outer_loop) :: list

      logical, intent(out) :: correct_flag
      integer, intent(out) :: ptypesize

      character(50), intent(in) :: testname

      integer, intent(in) :: filehandle_debug

      integer, intent(in) :: local_communicator

      double precision, dimension(DIM1+icount) :: atag, atype, amask
      double precision, dimension(3,DIM1+icount) :: ax
      double precision, dimension(:), allocatable :: buffer

      integer :: myrank, ier
      integer, parameter :: itag = 0
      integer :: i, j
      integer :: typesize, bytes, base, pos

      integer :: dtype_indexed1_t, dtype_indexed3_t, &
        dtype_send_t, dtype_cont1_t, dtype_cont3_t, dtype_recv_t
      integer(kind=MPI_ADDRESS_KIND), dimension(4) :: address_displacement
      integer, dimension(4) :: blocklength, oldtype
      integer, dimension(:), allocatable :: index_displacement

      character(50) :: method

!> just some statements to prevent compiler warnings of unused variables
!> those parameter are included for future features
      correct_flag = .false.
      ptypesize = 0
      typesize = filehandle_debug

      call MPI_Comm_rank( local_communicator, myrank, ier )

      base = myrank * (6*(DIM1+icount)) + 1
      call fill_unique_array_2D_double( ax, 3, DIM1+icount, base )
      base = base + 3*(DIM1+icount)
      call fill_unique_array_1D_double( atag, DIM1+icount, base )
      base = base + DIM1 + icount
      call fill_unique_array_1D_double( atype, DIM1+icount, base )
      base = base + DIM1 + icount
      call fill_unique_array_1D_double( amask, DIM1+icount, base )

      if ( myrank .EQ. 0 ) then
        write (method,'(A)') "mpi_pack_ddt"

        call MPI_Type_size( MPI_DOUBLE_PRECISION, typesize, ier )
        bytes = icount * 6 * typesize

        call timing_init( testname, method, bytes )
      endif

      do i=1,outer_loop

        call MPI_Type_size( MPI_DOUBLE_PRECISION, typesize, ier )
        bytes = 6 * icount * typesize 
        allocate( buffer(6*icount), stat=ier ) 

        allocate(index_displacement(icount), stat=ier)
        index_displacement(1:icount) = list(1:icount,i) - 1
        call MPI_Type_create_indexed_block( icount, 1, index_displacement, &
          MPI_DOUBLE_PRECISION, dtype_indexed1_t, ier )
      
        index_displacement(1:icount) = 3 * index_displacement(1:icount)
        call MPI_Type_create_indexed_block( icount, 3, index_displacement, &
          MPI_DOUBLE_PRECISION, dtype_indexed3_t, ier )

        call MPI_Get_address( ax(1,1), address_displacement(1), ier )
        call MPI_Get_address( atag(1), address_displacement(2), ier )
        call MPI_Get_address( atype(1), address_displacement(3), ier )
        call MPI_Get_address( amask(1), address_displacement(4), ier )

        oldtype(1) = dtype_indexed3_t
        oldtype(2:4) = dtype_indexed1_t
 
        blocklength = 1
 
        call MPI_Type_create_struct( 4, blocklength, address_displacement, &
          oldtype, dtype_send_t, ier )
        call MPI_Type_commit( dtype_send_t, ier )

        call MPI_Type_contiguous( icount, MPI_DOUBLE_PRECISION, &
          dtype_cont1_t, ier );
        call MPI_Type_contiguous( 3*icount, MPI_DOUBLE_PRECISION, &
          dtype_cont3_t, ier );

        call MPI_Get_address( ax(1,DIM1+1), address_displacement(1), ier )
        call MPI_Get_address( atag(DIM1+1), address_displacement(2), ier )
        call MPI_Get_address( atype(DIM1+1), address_displacement(3), ier )
        call MPI_Get_address( amask(DIM1+1), address_displacement(4), ier )

        oldtype(1) = dtype_cont3_t
        oldtype(2:4) = dtype_cont1_t

        call MPI_Type_create_struct( 4, blocklength, address_displacement, &
          oldtype, dtype_recv_t, ier )
        call MPI_Type_commit( dtype_recv_t, ier )

        call MPI_Type_free( dtype_indexed1_t, ier )
        call MPI_Type_free( dtype_indexed3_t, ier )
        call MPI_Type_free( dtype_cont1_t, ier )
        call MPI_Type_free( dtype_cont3_t, ier )

        deallocate(index_displacement, stat=ier )

        if ( myrank .EQ. 0 ) then
          call timing_record(1)
        endif

        do j=1,inner_loop

          if ( myrank .EQ. 0 ) then
            pos = 0
            call MPI_Pack( MPI_BOTTOM, 1, dtype_send_t, buffer, &
              bytes, pos, local_communicator, ier )
            call timing_record(2)
            call MPI_Send( buffer, pos, MPI_PACKED, 1, itag, &
              local_communicator, ier )
            call MPI_Recv( buffer, bytes, MPI_PACKED, 1, itag, &
              local_communicator, MPI_STATUS_IGNORE, ier )
            call timing_record(3)
            pos = 0
            call MPI_Unpack( buffer, bytes, pos, MPI_BOTTOM, 1, &
              dtype_recv_t, local_communicator, ier )
            call timing_record(4)
          else
            call MPI_Recv( buffer, bytes, MPI_PACKED, 0, itag, &
              local_communicator, MPI_STATUS_IGNORE, ier )
            pos = 0
            call MPI_Unpack( buffer, bytes, pos, MPI_BOTTOM, 1, &
              dtype_recv_t, local_communicator, ier )
            pos = 0
            call MPI_Pack( MPI_BOTTOM, 1, dtype_send_t, buffer, &
              bytes, pos, local_communicator, ier )           
            call MPI_Send( buffer, pos, MPI_PACKED, 0, itag, &
              local_communicator, ier )
          endif

        enddo ! inner loop

        deallocate( buffer, stat=ier )

        call MPI_Type_free( dtype_send_t, ier )
        call MPI_Type_free( dtype_recv_t, ier )

        if ( myrank .EQ. 0 ) then
          call timing_record(5)
        endif

      enddo ! outer loop

      if ( myrank .EQ. 0 ) then
        call timing_print( .true. )
      endif

      end subroutine timing_lammps_atomic_mpi_pack_ddt

      end module timing_lammps
