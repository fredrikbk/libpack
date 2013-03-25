      module timing_wrf

      private

      public :: timing_wrf_manual
      public :: timing_wrf_vec_ddt
      public :: timing_wrf_vec_mpi_pack_ddt
      public :: timing_wrf_sa_ddt
      public :: timing_wrf_sa_mpi_pack_ddt

      contains

      subroutine timing_wrf_manual ( &
        number_2D, number_3D, number_4D, &
        ims, ime, jms, jme, kms, kme, limit_4D_arrays, &
        is, ie, js, je, ks, ke, param_first_scalar, &
        outer_loop, inner_loop, correct_flag, ptypesize, testname, &
        filehandle_debug, local_communicator )

      use datatypes
!      use timing, only: timing_init, timing_record, timing_print
      use utilities, only: fill_unique_array_2D_real, &
        fill_unique_array_3D_real, fill_unique_array_4D_real

      implicit none

      include 'mpif.h'

      integer, intent(in) :: number_2D, number_3D, number_4D
      integer, intent(in) :: ims, ime, jms, jme, kms, kme
      integer, intent(in) :: is, ie, js, je, ks, ke

      integer, intent(in) :: param_first_scalar
      integer, dimension(number_4D), intent(in) :: limit_4D_arrays

      integer, intent(in) :: outer_loop, inner_loop

      logical, intent(out) :: correct_flag
      integer, intent(out) :: ptypesize

      character(50), intent(in) :: testname

      integer, intent(in) :: filehandle_debug
      integer, intent(in) :: local_communicator

      type (arrayPointer2D), allocatable :: array2Ds(:)
      type (arrayPointer3D), allocatable :: array3Ds(:)
      type (arrayPointer4D), allocatable :: array4Ds(:)

!> some local variables
      integer :: ier
      integer :: m, counter, ii, io, bytes, typesize, base
      integer :: i, j, k, l
      integer :: element_number
      integer :: myrank

      integer :: itag
      parameter (itag = 0)

      real, dimension(:), allocatable :: buffer

!> some variables for writing output
      character(50) :: method

!> just some statements to prevent compiler warnings of unused variables
!> those parameter are included for future features
      correct_flag = .false.
      ptypesize = 0
      typesize = filehandle_debug

      call MPI_Comm_rank( local_communicator, myrank, ier )

!> ================= initialize the arrays =================

!> allocate all needed arrays first
      allocate( array2Ds(number_2D), stat=ier )
      allocate( array3Ds(number_3D), stat=ier )
      allocate( array4Ds(number_4D), stat=ier )

!> allocate and initialize the arrays
!> compute the number of elements in the arrays
      counter = ( number_2D + number_3D * (kme-kms+1) ) * (ime-ims+1) &
        * (jme-jms+1)
      do m=1,number_4d
        counter = counter + limit_4D_arrays(m) * (ime-ims+1) * &
         (jme-jms+1) * (kme-kms+1)
      enddo
      base = myrank * counter + 1

      do m=1,number_2D
        allocate( array2Ds(m)%ptr(ims:ime,jms:jme), stat=ier )
        call fill_unique_array_2D_real( array2Ds(m)%ptr, (ime-ims+1), &
          (jme-jms+1), base )
        base = base + (ime-ims+1) * (jme-jms+1)
      enddo

      do m=1,number_3D
        allocate( array3Ds(m)%ptr(ims:ime,kms:kme,jms:jme), stat=ier )
        call fill_unique_array_3D_real( array3Ds(m)%ptr, (ime-ims+1), &
          (kme-kms+1), (jme-jms+1), base )
        base = base + (ime-ims+1) * (kme-kms+1) * (jme-jms+1)
      enddo

      do m=1,number_4D
        allocate (array4Ds(m)%ptr(ims:ime,kms:kme,jms:jme, &
          limit_4D_arrays(m)), stat=ier )
        call fill_unique_array_4D_real( array4Ds(m)%ptr, (ime-ims+1), &
          (kme-kms+1), (jme-jms+1), limit_4D_arrays(m), base )
        base = base + limit_4D_arrays(m) * (ime-ims+1) * (kme-kms+1) &
          * (jme-jms+1)
      enddo

      if ( myrank .EQ. 0 ) then
        write (method,'(A)') "manual"

!> compute the number of bytes to be communicated
!> first compute the number of elements in the subarrays

        counter = number_2D * (ie-is+1) * (je-js+1) + &
          number_3D * (ie-is+1) * (je-js+1) * (ke-ks+1)
        do m=1,number_4D
          if (limit_4D_arrays(m) .GE. param_first_scalar) then
            counter = counter+(limit_4D_arrays(m)-param_first_scalar+1) * &
              (ie-is+1) * (je-js+1) * (ke-ks+1)
          endif
        enddo
        call MPI_Type_size( MPI_REAL, typesize, ier )
        bytes = counter * typesize

        call timing_init( testname, method, bytes )
      endif

      do io=1,outer_loop

!> compute the number of elements in the subarray
        element_number = number_2D * (ie-is+1) * (je-js+1) + &
          number_3D * (ie-is+1) * (je-js+1) * (ke-ks+1)
        do m=1,number_4D
          if (limit_4D_arrays(m) .GE. param_first_scalar) then
            element_number = element_number + &
              (limit_4D_arrays(m)-param_first_scalar+1) * &
              (ie-is+1) * (je-js+1) * (ke-ks+1)
          endif
        enddo

        allocate( buffer(element_number), stat=ier )

        if ( myrank .EQ. 0 ) then
          call timing_record(1)
        endif

        do ii=1,inner_loop

!> =============== ping pong communication =================

!> send the data from rank 0 to rank 1 
          if ( myrank .EQ. 0 ) then
!> ==================== pack the data ======================
            counter = 1
            do m=1,number_2D
              do j=js,je
                do i=is,ie
                  buffer(counter) = array2Ds(m)%ptr(i,j)
                  counter = counter + 1
                enddo
              enddo
            enddo
            do m=1,number_3D
              do j=js,je
                do k=ks,ke
                  do i=is,ie
                    buffer(counter) = array3Ds(m)%ptr(i,k,j)
                    counter = counter + 1
                  enddo
                enddo
              enddo
            enddo
            do m=1,number_4D
              if (limit_4D_arrays(m) .GE. param_first_scalar) then
                do l=param_first_scalar,limit_4D_arrays(m)
                  do j=js,je
                    do k=ks,ke
                      do i=is,ie
                        buffer(counter) = array4Ds(m)%ptr(i,k,j,l)
                        counter = counter + 1
                      enddo
                    enddo
                  enddo
                enddo
              endif
            enddo
            call timing_record(2)
            call MPI_Send( buffer, element_number, MPI_REAL, &
              1, itag, local_communicator, &
              ier )
!> receive the data back from rank 1
            call MPI_Recv( buffer, element_number, MPI_REAL, &
              1, itag, local_communicator, &
              MPI_STATUS_IGNORE, ier )
            call timing_record(3)
!> =================== unpack the data =====================
            counter = 1
            do m=1,number_2D
              do j=js,je
                do i=is,ie
                  array2Ds(m)%ptr(i,j) = buffer(counter)
                  counter = counter + 1
                enddo
              enddo
            enddo
            do m=1,number_3D
              do j=js,je
                do k=ks,ke
                  do i=is,ie
                    array3Ds(m)%ptr(i,k,j) = buffer(counter)
                    counter = counter + 1
                  enddo
                enddo
              enddo
            enddo
            do m=1,number_4D
              if (limit_4D_arrays(m) .GE. param_first_scalar) then
                do l=param_first_scalar,limit_4D_arrays(m)
                  do j=js,je
                    do k=ks,ke
                      do i=is,ie
                        array4Ds(m)%ptr(i,k,j,l) = buffer(counter)
                        counter = counter + 1
                      enddo
                    enddo
                  enddo
                enddo
              endif
            enddo
            call timing_record(4)
!> now for rank 1
          else
!> receive from rank 0      
            call MPI_Recv( buffer, element_number, MPI_REAL, &
              0, itag, local_communicator, &
              MPI_STATUS_IGNORE, ier )
!> unpack the data 
            counter = 1
            do m=1,number_2D
              do j=js,je
                do i=is,ie
                  array2Ds(m)%ptr(i,j) = buffer(counter)
                  counter = counter + 1
                enddo
              enddo
            enddo

            do m=1,number_3D
              do j=js,je
                do k=ks,ke
                  do i=is,ie
                    array3Ds(m)%ptr(i,k,j) = buffer(counter)
                    counter = counter + 1
                  enddo
                enddo
              enddo
            enddo

            do m=1,number_4D
              if (limit_4D_arrays(m) .GE. param_first_scalar) then
                do l=param_first_scalar,limit_4D_arrays(m)
                  do j=js,je
                    do k=ks,ke
                      do i=is,ie
                        array4Ds(m)%ptr(i,k,j,l) = buffer(counter)
                        counter = counter + 1
                      enddo
                    enddo
                  enddo
                enddo
              endif
            enddo
!> pack the data
            counter = 1
            do m=1,number_2D
              do j=js,je
                do i=is,ie
                  buffer(counter) = array2Ds(m)%ptr(i,j)
                  counter = counter + 1
                enddo
              enddo
            enddo

            do m=1,number_3D
              do j=js,je
                do k=ks,ke
                  do i=is,ie
                    buffer(counter) = array3Ds(m)%ptr(i,k,j)
                    counter = counter + 1
                  enddo
                enddo
              enddo
            enddo

            do m=1,number_4D
              if (limit_4D_arrays(m) .GE. param_first_scalar) then
                do l=param_first_scalar,limit_4D_arrays(m)
                  do j=js,je
                    do k=ks,ke
                      do i=is,ie
                        buffer(counter) = array4Ds(m)%ptr(i,k,j,l)
                        counter = counter + 1
                      enddo
                    enddo
                  enddo
                enddo
              endif
            enddo

!> send to rank 0
            call MPI_Send( buffer, element_number, MPI_REAL, &
              0, itag, local_communicator, &
              ier )
          endif !> of myrank .EQ. 0?

        enddo !> inner_loop

!> ======================= clean up ========================

        deallocate( buffer, stat=ier )

        if ( myrank .EQ. 0 ) then 
          call timing_record(5)
        endif

      enddo !> outer_loop

      if ( myrank .EQ. 0 ) then
        call timing_print( .true. )
      endif

!> ======================= clean up ========================

      do m=1,number_2D
        deallocate( array2Ds(m)%ptr, stat=ier )
      enddo

      do m=1,number_3D
        deallocate( array3Ds(m)%ptr, stat=ier )
      enddo

      do m=1,number_4D
        deallocate (array4Ds(m)%ptr, stat=ier )
      enddo
     
      deallocate( array2Ds, stat=ier )
      deallocate( array3Ds, stat=ier )
      deallocate( array4Ds, stat=ier )

      end subroutine timing_wrf_manual

      subroutine timing_wrf_vec_ddt ( &
        number_2D, number_3D, number_4D, &
        ims, ime, jms, jme, kms, kme, limit_4D_arrays, &
        is, ie, js, je, ks, ke, param_first_scalar, &
        outer_loop, inner_loop, correct_flag, ptypesize, testname, &
        filehandle_debug, local_communicator )

      use datatypes
!      use timing, only: timing_init, timing_record, timing_print
      use utilities, only: fill_unique_array_2D_real, &
        fill_unique_array_3D_real, fill_unique_array_4D_real

      implicit none

      include 'mpif.h'

      integer, intent(in) :: number_2D, number_3D, number_4D
      integer, intent(in) :: ims, ime, jms, jme, kms, kme
      integer, intent(in) :: is, ie, js, je, ks, ke

      integer, intent(in) :: param_first_scalar
      integer, dimension(number_4D), intent(in) :: limit_4D_arrays

      integer, intent(in) :: outer_loop, inner_loop

      logical, intent(out) :: correct_flag
      integer, intent(out) :: ptypesize

      character(50), intent(in) :: testname

      integer, intent(in) :: filehandle_debug
      integer, intent(in) :: local_communicator

      type (arrayPointer2D), allocatable :: array2Ds(:)
      type (arrayPointer3D), allocatable :: array3Ds(:)
      type (arrayPointer4D), allocatable :: array4Ds(:)

!> some local variables
      integer :: ier
      integer :: m, counter, ii, io, bytes, base
      integer :: myrank

      integer :: itag
      parameter (itag = 0)

!> some variables for writing output
      character(50) :: method

!> variables for the MPI derived datatypes
      integer :: dtype_subarray_t, dtype_temp_2D_t
      integer :: dtype_temp_t, dtype_temp_3D_t
      integer, dimension(:), allocatable :: oldtype, blocklength
      integer(kind=MPI_ADDRESS_KIND), dimension(:), allocatable :: &
        displacement
      integer(kind=MPI_ADDRESS_KIND) :: stride
      integer :: typesize

!> just some statements to prevent compiler warnings of unused variables
!> those parameter are included for future features
      correct_flag = .false.
      ptypesize = 0
      typesize = filehandle_debug

      call MPI_Comm_rank( local_communicator, myrank, ier )

!> ================= initialize the arrays =================

!> allocate all needed arrays first
      allocate( array2Ds(number_2D), stat=ier )
      allocate( array3Ds(number_3D), stat=ier )
      allocate( array4Ds(number_4D), stat=ier )

!> allocate and initialize the arrays
!> compute the number of elements in the arrays
      counter = ( number_2D + number_3D * (kme-kms+1) ) * (ime-ims+1) &
        * (jme-jms+1)
      do m=1,number_4d
        counter = counter + limit_4D_arrays(m) * (ime-ims+1) * &
         (jme-jms+1) * (kme-kms+1)
      enddo
      base = myrank * counter + 1
 
      do m=1,number_2D
        allocate( array2Ds(m)%ptr(ims:ime,jms:jme), stat=ier )
        call fill_unique_array_2D_real( array2Ds(m)%ptr, (ime-ims+1), &
          (jme-jms+1), base )
        base = base + (ime-ims+1) * (jme-jms+1)
      enddo

      do m=1,number_3D
        allocate( array3Ds(m)%ptr(ims:ime,kms:kme,jms:jme), stat=ier )
        call fill_unique_array_3D_real( array3Ds(m)%ptr, (ime-ims+1), &
          (kme-kms+1), (jme-jms+1), base )
        base = base + (ime-ims+1) * (kme-kms+1) * (jme-jms+1)
      enddo

      do m=1,number_4D
        allocate (array4Ds(m)%ptr(ims:ime,kms:kme,jms:jme, &
          limit_4D_arrays(m)), stat=ier )
        call fill_unique_array_4D_real( array4Ds(m)%ptr, (ime-ims+1), &
          (kme-kms+1), (jme-jms+1), limit_4D_arrays(m), base )
        base = base + limit_4D_arrays(m) * (ime-ims+1) * (kme-kms+1) &
          * (jme-jms+1)
      enddo

      if ( myrank .EQ. 0 ) then
        write (method,'(A)') "mpi_ddt"

!> compute the number of bytes to be communicated
!> first compute the number of elements in the subarrays

        counter = number_2D * (ie-is+1) * (je-js+1) + &
          number_3D * (ie-is+1) * (je-js+1) * (ke-ks+1)
        do m=1,number_4D
          if (limit_4D_arrays(m) .GE. param_first_scalar) then
            counter = counter+(limit_4D_arrays(m)-param_first_scalar+1) * &
              (ie-is+1) * (je-js+1) * (ke-ks+1)
          endif
        enddo
        call MPI_Type_size( MPI_REAL, typesize, ier )
        bytes = counter * typesize

        call timing_init( testname, method, bytes )
      endif

      do io=1,outer_loop

!> ========== building the MPI derived datatype ============

!> compute the number of subarrays, so that the buffer for the datatype
!> creation are big enough

        counter = number_2D + number_3D
        do m=1,number_4D
          if ( limit_4D_arrays(m) .GE. param_first_scalar ) then
            counter = counter + 1
          endif
        enddo
 
        allocate( oldtype(counter), stat=ier )
        allocate( blocklength(counter), stat=ier )
        allocate( displacement(counter), stat=ier )

        blocklength(1:counter) = 1
  
        counter = 1
        
        call MPI_Type_vector( je-js+1, ie-is+1, ime-ims+1, MPI_REAL, &
          dtype_temp_2D_t, ier )

        do m=1,number_2D
          call MPI_Get_address( array2Ds(m)%ptr(is,js), &
            displacement(counter), ier )
          oldtype(counter) = dtype_temp_2D_t
          counter = counter + 1
        enddo

        call MPI_Type_size (MPI_REAL, typesize, ier )
        call MPI_Type_vector( ke-ks+1, ie-is+1, ime-ims+1, MPI_REAL, &
          dtype_temp_t, ier )
        stride = (ime-ims+1)*(kme-kms+1) * typesize
        call MPI_Type_create_hvector( je-js+1, 1, stride, &
          dtype_temp_t, dtype_temp_3D_t, ier )
        call MPI_Type_free( dtype_temp_t, ier )

        do m=1,number_3D
          call MPI_Get_address( array3Ds(m)%ptr(is,ks,js), &
            displacement(counter), ier )
          oldtype(counter) = dtype_temp_3D_t
          counter = counter + 1
        enddo

        stride = stride * (jme-jms+1)
        do m=1,number_4D
          if (limit_4D_arrays(m) .GE. param_first_scalar) then
            call MPI_Get_address( &
              array4Ds(m)%ptr(is,ks,js,param_first_scalar), &
              displacement(counter), ier )
            call MPI_Type_create_hvector( &
              limit_4D_arrays(m)-param_first_scalar+1, &
              1, stride, dtype_temp_3D_t, oldtype(counter), ier )
            counter = counter + 1
          endif 
        enddo

        counter = counter - 1 !> we need the actual number of arrays
        call MPI_Type_create_struct( counter, blocklength, displacement, &
          oldtype, dtype_subarray_t, ier )
        call MPI_Type_commit( dtype_subarray_t, ier )

        call MPI_Type_free( dtype_temp_2D_t, ier )
        call MPI_Type_free( dtype_temp_3D_t, ier )
        counter = number_2D + number_3D + 1
        do m=1,number_4D
          if (limit_4D_arrays(m) .GE. param_first_scalar) then
            call MPI_Type_free( oldtype(counter), ier )
            counter = counter + 1
          endif
        enddo

        deallocate( oldtype, stat=ier )
        deallocate( blocklength, stat=ier )
        deallocate( displacement, stat=ier )

        if ( myrank .EQ. 0 ) then
          call timing_record(1)
        endif

        do ii=1,inner_loop

!> =============== ping pong communication =================

!> send the data from rank 0 to rank 1 
          if ( myrank .EQ. 0 ) then
            call MPI_Send( MPI_BOTTOM, 1, dtype_subarray_t, &
              1, itag, local_communicator, ier )
!> receive the data back from rank 1
            call MPI_Recv( MPI_BOTTOM, 1, dtype_subarray_t, &
              1, itag, local_communicator, MPI_STATUS_IGNORE, ier )
            call timing_record(3)
!> now for rank 1
          else
!> receive from rank 0      
            call MPI_Recv( MPI_BOTTOM, 1, dtype_subarray_t, &
              0, itag, local_communicator, MPI_STATUS_IGNORE, ier )
!> send to rank 0
            call MPI_Send( MPI_BOTTOM, 1, dtype_subarray_t, &
              0, itag, local_communicator, ier )
          endif

        enddo !> inner_loop

!> ======================= clean up ========================

        call MPI_Type_free( dtype_subarray_t, ier )

        if ( myrank .EQ. 0 ) then 
          call timing_record(5)
        endif

      enddo !> outer_loop

      if ( myrank .EQ. 0 ) then
        call timing_print( .true. )
      endif

!> ======================= clean up ========================

      do m=1,number_2D
        deallocate( array2Ds(m)%ptr, stat=ier )
      enddo

      do m=1,number_3D
        deallocate( array3Ds(m)%ptr, stat=ier )
      enddo

      do m=1,number_4D
        deallocate (array4Ds(m)%ptr, stat=ier )
      enddo
     
      deallocate( array2Ds, stat=ier )
      deallocate( array3Ds, stat=ier )
      deallocate( array4Ds, stat=ier )

      end subroutine timing_wrf_vec_ddt

      subroutine timing_wrf_vec_mpi_pack_ddt ( &
        number_2D, number_3D, number_4D, &
        ims, ime, jms, jme, kms, kme, limit_4D_arrays, &
        is, ie, js, je, ks, ke, param_first_scalar, &
        outer_loop, inner_loop, correct_flag, ptypesize, testname, &
        filehandle_debug, local_communicator )

      use datatypes
!      use timing, only: timing_init, timing_record, timing_print
      use utilities, only: fill_unique_array_2D_real, &
        fill_unique_array_3D_real, fill_unique_array_4D_real

      implicit none

      include 'mpif.h'

      integer, intent(in) :: number_2D, number_3D, number_4D
      integer, intent(in) :: ims, ime, jms, jme, kms, kme
      integer, intent(in) :: is, ie, js, je, ks, ke

      integer, intent(in) :: param_first_scalar
      integer, dimension(number_4D), intent(in) :: limit_4D_arrays

      integer, intent(in) :: outer_loop, inner_loop

      logical, intent(out) :: correct_flag
      integer, intent(out) :: ptypesize

      character(50), intent(in) :: testname

      integer, intent(in) :: filehandle_debug
      integer, intent(in) :: local_communicator

      type (arrayPointer2D), allocatable :: array2Ds(:)
      type (arrayPointer3D), allocatable :: array3Ds(:)
      type (arrayPointer4D), allocatable :: array4Ds(:)

!> some local variables
      integer :: ier
      integer :: m, counter, ii, io, bytes, typesize, pos, base
      integer :: element_number
      integer :: myrank

      integer :: itag
      parameter (itag = 0)

      real, dimension(:), allocatable :: buffer

!> some variables for writing output
      character(50) :: method

!> variables for the MPI derived datatypes
      integer :: dtype_subarray_t, dtype_temp_2D_t
      integer :: dtype_temp_t, dtype_temp_3D_t
      integer, dimension(:), allocatable :: oldtype, blocklength
      integer(kind=MPI_ADDRESS_KIND), dimension(:), allocatable :: &
        displacement
      integer(kind=MPI_ADDRESS_KIND) :: stride

!> just some statements to prevent compiler warnings of unused variables
!> those parameter are included for future features
      correct_flag = .false.
      ptypesize = 0
      typesize = filehandle_debug

      call MPI_Comm_rank( local_communicator, myrank, ier )

!> ================= initialize the arrays =================

!> allocate all needed arrays first
      allocate( array2Ds(number_2D), stat=ier )
      allocate( array3Ds(number_3D), stat=ier )
      allocate( array4Ds(number_4D), stat=ier )

!> allocate and initialize the arrays
!> compute the number of elements in the arrays
      counter = ( number_2D + number_3D * (kme-kms+1) ) * (ime-ims+1) &
        * (jme-jms+1)
      do m=1,number_4d
        counter = counter + limit_4D_arrays(m) * (ime-ims+1) * &
         (jme-jms+1) * (kme-kms+1)
      enddo
      base = myrank * counter + 1

      do m=1,number_2D
        allocate( array2Ds(m)%ptr(ims:ime,jms:jme), stat=ier )
        call fill_unique_array_2D_real( array2Ds(m)%ptr, (ime-ims+1), &
          (jme-jms+1), base )
        base = base + (ime-ims+1) * (jme-jms+1)
      enddo

      do m=1,number_3D
        allocate( array3Ds(m)%ptr(ims:ime,kms:kme,jms:jme), stat=ier )
        call fill_unique_array_3D_real( array3Ds(m)%ptr, (ime-ims+1), &
          (kme-kms+1), (jme-jms+1), base )
        base = base + (ime-ims+1) * (kme-kms+1) * (jme-jms+1)
      enddo

      do m=1,number_4D
        allocate (array4Ds(m)%ptr(ims:ime,kms:kme,jms:jme, &
          limit_4D_arrays(m)), stat=ier )
        call fill_unique_array_4D_real( array4Ds(m)%ptr, (ime-ims+1), &
          (kme-kms+1), (jme-jms+1), limit_4D_arrays(m), base )
        base = base + limit_4D_arrays(m) * (ime-ims+1) * (kme-kms+1) &
          * (jme-jms+1)
      enddo

      if ( myrank .EQ. 0 ) then
        write (method,'(A)') "mpi_pack_ddt"

!> compute the number of bytes to be communicated
!> first compute the number of elements in the subarrays

        counter = number_2D * (ie-is+1) * (je-js+1) + &
          number_3D * (ie-is+1) * (je-js+1) * (ke-ks+1)
        do m=1,number_4D
          if (limit_4D_arrays(m) .GE. param_first_scalar) then
            counter = counter+(limit_4D_arrays(m)-param_first_scalar+1) * &
              (ie-is+1) * (je-js+1) * (ke-ks+1)
          endif
        enddo
        
        call MPI_Type_size( MPI_REAL, typesize, ier )
        bytes = counter * typesize

        call timing_init( testname, method, bytes )
      endif

      do io=1,outer_loop

!> compute the number of elements in the subarray
        element_number = number_2D * (ie-is+1) * (je-js+1) + &
          number_3D * (ie-is+1) * (je-js+1) * (ke-ks+1)
        do m=1,number_4D
          if (limit_4D_arrays(m) .GE. param_first_scalar) then
            element_number = element_number + &
              (limit_4D_arrays(m)-param_first_scalar+1) * &
              (ie-is+1) * (je-js+1) * (ke-ks+1)
          endif
        enddo

        allocate( buffer(element_number), stat=ier )
        call MPI_Type_size( MPI_REAL, typesize, ier )
        bytes = element_number * typesize

!> compute the number of arrays, so that the buffer for the datatype
!> creation are big enough
        counter = number_2D + number_3D
        do m=1,number_4D
          if ( limit_4D_arrays(m) .GE. param_first_scalar ) then
            counter = counter + 1
          endif
        enddo
 
        allocate( oldtype(counter), stat=ier )
        allocate( blocklength(counter), stat=ier )
        allocate( displacement(counter), stat=ier )

        blocklength(1:counter) = 1
  
        call MPI_Type_vector( je-js+1, ie-is+1, ime-ims+1, MPI_REAL, &
          dtype_temp_2D_t, ier )
 
        counter = 1
        do m=1,number_2D
          call MPI_Get_address( array2Ds(m)%ptr(is,js), &
            displacement(counter), ier )
          oldtype(counter) = dtype_temp_2D_t
          counter = counter + 1
        enddo

        call MPI_Type_size (MPI_REAL, typesize, ier )
        call MPI_Type_vector( ke-ks+1, ie-is+1, ime-ims+1, MPI_REAL, &
          dtype_temp_t, ier )
        stride = (ime-ims+1)*(kme-kms+1) * typesize
        call MPI_Type_create_hvector( je-js+1, 1, stride, &
          dtype_temp_t, dtype_temp_3D_t, ier )
        call MPI_Type_free( dtype_temp_t, ier )

        do m=1,number_3D
          call MPI_Get_address( array3Ds(m)%ptr(is,ks,js), &
            displacement(counter), ier )
          oldtype(counter) = dtype_temp_3D_t
          counter = counter + 1
        enddo

        stride = stride*(jme-jms+1)
        do m=1,number_4D
          if (limit_4D_arrays(m) .GE. param_first_scalar) then
            call MPI_Get_address( &
              array4Ds(m)%ptr(is,ks,js,param_first_scalar), &
              displacement(counter), ier )
            call MPI_Type_create_hvector( &
              limit_4D_arrays(m)-param_first_scalar+1, &
              1, stride, dtype_temp_3D_t, oldtype(counter), ier )
            counter = counter + 1
          endif 
        enddo

        counter = counter - 1 !> we need the actual number of arrays
        call MPI_Type_create_struct( counter, blocklength, displacement, &
          oldtype, dtype_subarray_t, ier )
        call MPI_Type_commit( dtype_subarray_t, ier )

        call MPI_Type_free( dtype_temp_2D_t, ier )
        call MPI_Type_free( dtype_temp_3D_t, ier )
        counter = number_2D + number_3D + 1
        do m=1,number_4D
          if (limit_4D_arrays(m) .GE. param_first_scalar) then
            call MPI_Type_free( oldtype(counter), ier )
            counter = counter + 1
          endif
        enddo

        deallocate( oldtype, stat=ier )
        deallocate( blocklength, stat=ier )
        deallocate( displacement, stat=ier )

        if ( myrank .EQ. 0 ) then
          call timing_record(1)
        endif

        do ii=1,inner_loop

!> =============== ping pong communication =================

!> send the data from rank 0 to rank 1 
          if ( myrank .EQ. 0 ) then
!> ==================== pack the data ======================
            pos = 0
            call MPI_Pack( MPI_BOTTOM, 1, dtype_subarray_t, &
              buffer, bytes, pos, local_communicator, ier )
            call timing_record(2)
            call MPI_Send( buffer, pos, MPI_PACKED, &
              1, itag, local_communicator, ier )
!> receive the data back from rank 1
            call MPI_Recv( buffer, bytes, MPI_PACKED, &
              1, itag, local_communicator, MPI_STATUS_IGNORE, ier )
            call timing_record(3)
!> =================== unpack the data =====================
            pos = 0
            call MPI_Unpack( buffer, bytes, pos, MPI_BOTTOM, 1, &
              dtype_subarray_t, local_communicator, ier )
            call timing_record(4)
!> now for rank 1
          else
!> receive from rank 0      
            call MPI_Recv( buffer, element_number, MPI_REAL, &
              0, itag, local_communicator, &
              MPI_STATUS_IGNORE, ier )
!> unpack the data 
            pos = 0
            call MPI_Unpack( buffer, bytes, pos, MPI_BOTTOM, 1, &
              dtype_subarray_t, local_communicator, ier )
!> pack the data
            pos = 0
            call MPI_Pack( MPI_BOTTOM, 1, dtype_subarray_t, &
              buffer, bytes, pos, local_communicator, ier )
!> send to rank 0
            call MPI_Send( buffer, element_number, MPI_REAL, &
              0, itag, local_communicator, &
              ier )
          endif !> of myrank .EQ. 0?

        enddo !> inner_loop

!> ======================= clean up ========================

        call MPI_Type_free( dtype_subarray_t, ier )

        deallocate( buffer, stat=ier )

        if ( myrank .EQ. 0 ) then 
          call timing_record(5)
        endif

      enddo !> outer_loop

      if ( myrank .EQ. 0 ) then
        call timing_print( .true. )
      endif

!> ======================= clean up ========================

      do m=1,number_2D
        deallocate( array2Ds(m)%ptr, stat=ier )
      enddo

      do m=1,number_3D
        deallocate( array3Ds(m)%ptr, stat=ier )
      enddo

      do m=1,number_4D
        deallocate (array4Ds(m)%ptr, stat=ier )
      enddo
     
      deallocate( array2Ds, stat=ier )
      deallocate( array3Ds, stat=ier )
      deallocate( array4Ds, stat=ier )

      end subroutine timing_wrf_vec_mpi_pack_ddt

      subroutine timing_wrf_sa_ddt ( &
        number_2D, number_3D, number_4D, &
        ims, ime, jms, jme, kms, kme, limit_4D_arrays, &
        is, ie, js, je, ks, ke, param_first_scalar, &
        outer_loop, inner_loop, correct_flag, ptypesize, testname, &
        filehandle_debug, local_communicator )

      use datatypes
!      use timing, only: timing_init, timing_record, timing_print
      use utilities, only: fill_unique_array_2D_real, &
        fill_unique_array_3D_real, fill_unique_array_4D_real

      implicit none

      include 'mpif.h'

      integer, intent(in) :: number_2D, number_3D, number_4D
      integer, intent(in) :: ims, ime, jms, jme, kms, kme
      integer, intent(in) :: is, ie, js, je, ks, ke

      integer, intent(in) :: param_first_scalar
      integer, dimension(number_4D), intent(in) :: limit_4D_arrays

      integer, intent(in) :: outer_loop, inner_loop

      logical, intent(out) :: correct_flag
      integer, intent(out) :: ptypesize

      character(50), intent(in) :: testname

      integer, intent(in) :: filehandle_debug
      integer, intent(in) :: local_communicator

      type (arrayPointer2D), allocatable :: array2Ds(:)
      type (arrayPointer3D), allocatable :: array3Ds(:)
      type (arrayPointer4D), allocatable :: array4Ds(:)

! some local variables
      integer :: ier
      integer :: m, counter, ii, io, bytes, base
      integer :: myrank

      integer :: itag
      parameter (itag = 0)

! some variables for writing output
      character(50) :: method

! variables for the MPI derived datatype
      integer :: dtype_subarray_t, dtype_temp_2D_t
      integer :: dtype_temp_3D_t
      integer, dimension(:), allocatable :: oldtype, blocklength
      integer(kind=MPI_ADDRESS_KIND), dimension(:), allocatable :: &
        displacement
      integer :: typesize
      integer, dimension(4) :: arraysize, subarraysize, subarraystart

!> just some statements to prevent compiler warnings of unused variables
!> those parameter are included for future features
      correct_flag = .false.
      ptypesize = 0
      typesize = filehandle_debug

      call MPI_Comm_rank( local_communicator, myrank, ier )

!> ================= initialize the arrays =================

!> allocate all needed arrays first
      allocate( array2Ds(number_2D), stat=ier )
      allocate( array3Ds(number_3D), stat=ier )
      allocate( array4Ds(number_4D), stat=ier )

!> allocate and initialize the arrays
!> compute the number of elements in the arrays
      counter = ( number_2D + number_3D * (kme-kms+1) ) * (ime-ims+1) &
        * (jme-jms+1)
      do m=1,number_4d
        counter = counter + limit_4D_arrays(m) * (ime-ims+1) * &
          (jme-jms+1) * (kme-kms+1)
      enddo
      base = myrank * counter + 1

      do m=1,number_2D
        allocate( array2Ds(m)%ptr(ims:ime,jms:jme), stat=ier )
        call fill_unique_array_2D_real( array2Ds(m)%ptr, (ime-ims+1), &
          (jme-jms+1), base )
        base = base + (ime-ims+1) * (jme-jms+1)
      enddo

      do m=1,number_3D
        allocate( array3Ds(m)%ptr(ims:ime,kms:kme,jms:jme), stat=ier )
        call fill_unique_array_3D_real( array3Ds(m)%ptr, (ime-ims+1), &
          (kme-kms+1), (jme-jms+1), base )
        base = base + (ime-ims+1) * (kme-kms+1) * (jme-jms+1)
      enddo

      do m=1,number_4D
        allocate (array4Ds(m)%ptr(ims:ime,kms:kme,jms:jme, &
          limit_4D_arrays(m)), stat=ier )
        call fill_unique_array_4D_real( array4Ds(m)%ptr, (ime-ims+1), &
          (kme-kms+1), (jme-jms+1), limit_4D_arrays(m), base )
        base = base + limit_4D_arrays(m) * (ime-ims+1) * (kme-kms+1) &
          * (jme-jms+1)
      enddo

      if ( myrank .EQ. 0 ) then
        write (method,'(A)') "mpi_ddt"

!> compute the number of bytes to be communicated
!> first compute the number of elements in the subarrays

        counter = number_2D * (ie-is+1) * (je-js+1) + &
          number_3D * (ie-is+1) * (je-js+1) * (ke-ks+1)
        do m=1,number_4D
          if (limit_4D_arrays(m) .GE. param_first_scalar) then
            counter = counter+(limit_4D_arrays(m)-param_first_scalar+1) &
              * (ie-is+1) * (je-js+1) * (ke-ks+1)
          endif
        enddo
        call MPI_Type_size( MPI_REAL, typesize, ier )
        bytes = counter * typesize

        call timing_init( testname, method, bytes )
      endif

      do io=1,outer_loop

!> ========== building the MPI derived datatype ============

!> compute the number of subarrays, so that the buffer for the datatype
!> creation are big enough
        counter = number_2D + number_3D
        do m=1,number_4D
          if ( limit_4D_arrays(m) .GE. param_first_scalar ) then
            counter = counter + 1
          endif
        enddo
 
        allocate( oldtype(counter), stat=ier )
        allocate( blocklength(counter), stat=ier )
        allocate( displacement(counter), stat=ier )

        blocklength(1:counter) = 1
  
        counter = 1
       
!> create the 2D subarray type 
        arraysize(1) = ime-ims+1
        arraysize(2) = jme-jms+1
        subarraysize(1) = ie-is+1
        subarraysize(2) = je-js+1
        subarraystart(1) = is-ims
        subarraystart(2) = js-jms
        call MPI_Type_create_subarray( 2, arraysize, subarraysize, &
          subarraystart, MPI_ORDER_FORTRAN, MPI_REAL, &
          dtype_temp_2D_t, ier )

!> create the parameters of the 2D arrays for the struct type
        do m=1,number_2D
          call MPI_Get_address( array2Ds(m)%ptr(ims,jms), &
            displacement(counter), ier )
          oldtype(counter) = dtype_temp_2D_t
          counter = counter + 1
        enddo

!> create the 3D subarray type
!>        arraysize(1) = ime-ims+1, already set
        arraysize(2) = kme-kms+1
        arraysize(3) = jme-jms+1
!>        subarraysize(1) = ie-is+1, already set
        subarraysize(2) = ke-ks+1
        subarraysize(3) = je-js+1
!>        subarraystart(1) = is-ims, already set
        subarraystart(2) = ks-kms
        subarraystart(3) = js-jms
        call MPI_Type_create_subarray( 3, arraysize, subarraysize, &
          subarraystart, MPI_ORDER_FORTRAN, MPI_REAL, &
          dtype_temp_3D_t, ier )

        do m=1,number_3D
          call MPI_Get_address( array3Ds(m)%ptr(ims,kms,jms), &
            displacement(counter), ier )
          oldtype(counter) = dtype_temp_3D_t
          counter = counter + 1
        enddo

!> create the 4D subarray types
        subarraystart(4) = param_first_scalar-1
        do m=1,number_4D
          if (limit_4D_arrays(m) .GE. param_first_scalar) then
            arraysize(4) = limit_4D_arrays(m)
            subarraysize(4) = limit_4D_arrays(m)-param_first_scalar+1
            call MPI_Type_create_subarray( 4, arraysize, subarraysize, &
              subarraystart, MPI_ORDER_FORTRAN, MPI_REAL, &
              oldtype(counter), ier )
            call MPI_Get_address( &
              array4Ds(m)%ptr(ims,kms,jms,1), &
              displacement(counter), ier )
            counter = counter + 1
          endif 
        enddo

!> create the all-embracing struct type
        counter = counter - 1 !> we need the actual number of arrays
        call MPI_Type_create_struct( counter, blocklength, displacement, &
          oldtype, dtype_subarray_t, ier )
        call MPI_Type_commit( dtype_subarray_t, ier )

        call MPI_Type_free( dtype_temp_2D_t, ier )
        call MPI_Type_free( dtype_temp_3D_t, ier )
        counter = number_2D + number_3D + 1
        do m=1,number_4D
          if (limit_4D_arrays(m) .GE. param_first_scalar) then
            call MPI_Type_free( oldtype(counter), ier )
            counter = counter + 1
          endif
        enddo

        deallocate( oldtype, stat=ier )
        deallocate( blocklength, stat=ier )
        deallocate( displacement, stat=ier )

        if ( myrank .EQ. 0 ) then
          call timing_record(1)
        endif

        do ii=1,inner_loop

!> =============== ping pong communication =================

!> send the data from rank 0 to rank 1 
          if ( myrank .EQ. 0 ) then
            call MPI_Send( MPI_BOTTOM, 1, dtype_subarray_t, &
              1, itag, local_communicator, &
              ier )
!> receive the data back from rank 1
            call MPI_Recv( MPI_BOTTOM, 1, dtype_subarray_t, &
              1, itag, local_communicator, &
              MPI_STATUS_IGNORE, ier )
            call timing_record(3)
!> now for rank 1
          else
!> receive from rank 0      
            call MPI_Recv( MPI_BOTTOM, 1, dtype_subarray_t, &
              0, itag, local_communicator, &
              MPI_STATUS_IGNORE, ier )
!> send to rank 0
            call MPI_Send( MPI_BOTTOM, 1, dtype_subarray_t, &
              0, itag, local_communicator, &
              ier )
          endif

        enddo !> inner_loop

!> ======================= clean up ========================

        call MPI_Type_free( dtype_subarray_t, ier )

        if ( myrank .EQ. 0 ) then 
          call timing_record(5)
        endif

      enddo !> outer_loop

      if ( myrank .EQ. 0 ) then
        call timing_print( .true. )
      endif

!> ======================= clean up ========================

      do m=1,number_2D
        deallocate( array2Ds(m)%ptr, stat=ier )
      enddo

      do m=1,number_3D
        deallocate( array3Ds(m)%ptr, stat=ier )
      enddo

      do m=1,number_4D
        deallocate (array4Ds(m)%ptr, stat=ier )
      enddo
     
      deallocate( array2Ds, stat=ier )
      deallocate( array3Ds, stat=ier )
      deallocate( array4Ds, stat=ier )

      end subroutine timing_wrf_sa_ddt

      subroutine timing_wrf_sa_mpi_pack_ddt ( &
        number_2D, number_3D, number_4D, &
        ims, ime, jms, jme, kms, kme, limit_4D_arrays, &
        is, ie, js, je, ks, ke, param_first_scalar, &
        outer_loop, inner_loop, correct_flag, ptypesize, testname, &
        filehandle_debug, local_communicator )

      use datatypes
!      use timing, only: timing_init, timing_record, timing_print
      use utilities, only: fill_unique_array_2D_real, &
        fill_unique_array_3D_real, fill_unique_array_4D_real

      implicit none

      include 'mpif.h'

      integer, intent(in) :: number_2D, number_3D, number_4D
      integer, intent(in) :: ims, ime, jms, jme, kms, kme
      integer, intent(in) :: is, ie, js, je, ks, ke

      integer, intent(in) :: param_first_scalar
      integer, dimension(number_4D), intent(in) :: limit_4D_arrays

      integer, intent(in) :: outer_loop, inner_loop

      logical, intent(out) :: correct_flag
      integer, intent(out) :: ptypesize

      character(50), intent(in) :: testname

      integer, intent(in) :: filehandle_debug
      integer, intent(in) :: local_communicator

      type (arrayPointer2D), allocatable :: array2Ds(:)
      type (arrayPointer3D), allocatable :: array3Ds(:)
      type (arrayPointer4D), allocatable :: array4Ds(:)

!> some local variables
      integer :: ier
      integer :: m, counter, ii, io, bytes, typesize, pos, base
      integer :: element_number
      integer :: myrank

      integer :: itag
      parameter (itag = 0)

      real, dimension(:), allocatable :: buffer

!> some variables for writing output
      character(50) :: method

!> variables for the MPI derived datatypes
      integer :: dtype_subarray_t, dtype_temp_2D_t, dtype_temp_3D_t
      integer, dimension(:), allocatable :: oldtype, blocklength
      integer(kind=MPI_ADDRESS_KIND), dimension(:), allocatable :: &
        displacement
      integer, dimension(4) :: arraysize, subarraysize, subarraystart

!> just some statements to prevent compiler warnings of unused variables
!> those parameter are included for future features
      correct_flag = .false.
      ptypesize = 0
      typesize = filehandle_debug

      call MPI_Comm_rank( local_communicator, myrank, ier )

!> ================= initialize the arrays =================

!> allocate all needed arrays first
      allocate( array2Ds(number_2D), stat=ier )
      allocate( array3Ds(number_3D), stat=ier )
      allocate( array4Ds(number_4D), stat=ier )
 
!> allocate and initialize the arrays
!> compute the number of elements in the arrays
      counter = ( number_2D + number_3D * (kme-kms+1) ) * (ime-ims+1) &
        * (jme-jms+1)
      do m=1,number_4d
        counter = counter + limit_4D_arrays(m) * (ime-ims+1) * &
          (jme-jms+1) * (kme-kms+1)
      enddo
      base = myrank * counter + 1
 
      do m=1,number_2D
        allocate( array2Ds(m)%ptr(ims:ime,jms:jme), stat=ier )
        call fill_unique_array_2D_real( array2Ds(m)%ptr, (ime-ims+1), &
          (jme-jms+1), base )
        base = base + (ime-ims+1) * (jme-jms+1) 
      enddo

      do m=1,number_3D
        allocate( array3Ds(m)%ptr(ims:ime,kms:kme,jms:jme), stat=ier )
        call fill_unique_array_3D_real( array3Ds(m)%ptr, (ime-ims+1), &
          (kme-kms+1), (jme-jms+1), base )
        base = base + (ime-ims+1) * (kme-kms+1) * (jme-jms+1)
      enddo

      do m=1,number_4D
        allocate (array4Ds(m)%ptr(ims:ime,kms:kme,jms:jme, &
          limit_4D_arrays(m)), stat=ier )
        call fill_unique_array_4D_real( array4Ds(m)%ptr, (ime-ims+1), &
          (kme-kms+1), (jme-jms+1), limit_4D_arrays(m), base )
        base = base + limit_4D_arrays(m) * (ime-ims+1) * (kme-kms+1) &
          * (jme-jms+1)
      enddo

      if ( myrank .EQ. 0 ) then
        write (method,'(A)') "mpi_pack_ddt"

!> compute the number of bytes to be communicated
!> first compute the number of elements in the subarrays

        counter = number_2D * (ie-is+1) * (je-js+1) + &
          number_3D * (ie-is+1) * (je-js+1) * (ke-ks+1)
        do m=1,number_4D
          if (limit_4D_arrays(m) .GE. param_first_scalar) then
            counter = counter+(limit_4D_arrays(m)-param_first_scalar+1) &
              * (ie-is+1) * (je-js+1) * (ke-ks+1)
          endif
        enddo
        call MPI_Type_size( MPI_REAL, typesize, ier )
        bytes = counter * typesize

        call timing_init( testname, method, bytes )
      endif

      do io=1,outer_loop

!> compute the number of elements in the subarray
        element_number = number_2D * (ie-is+1) * (je-js+1) + &
          number_3D * (ie-is+1) * (je-js+1) * (ke-ks+1)
        do m=1,number_4D
          if (limit_4D_arrays(m) .GE. param_first_scalar) then
            element_number = element_number + &
              (limit_4D_arrays(m)-param_first_scalar+1) * &
              (ie-is+1) * (je-js+1) * (ke-ks+1)
          endif
        enddo

        allocate( buffer(element_number), stat=ier )
        call MPI_Type_size( MPI_REAL, typesize, ier )
        bytes = element_number * typesize

!> compute the number of subarrays, so that buffer for the datatype
!> creation are big enough
        counter = number_2D + number_3D
        do m=1,number_4D
          if ( limit_4D_arrays(m) .GE. param_first_scalar ) then
            counter = counter + 1
          endif
        enddo
 
        allocate( oldtype(counter), stat=ier )
        allocate( blocklength(counter), stat=ier )
        allocate( displacement(counter), stat=ier )

        blocklength(1:counter) = 1

        counter = 1
 
!> create the subtype for the 2D case
        arraysize(1) = ime-ims+1
        arraysize(2) = jme-jms+1
        subarraysize(1) = ie-is+1
        subarraysize(2) = je-js+1
        subarraystart(1) = is-ims
        subarraystart(2) = js-jms
        call MPI_Type_create_subarray( 2, arraysize, subarraysize, &
          subarraystart, MPI_ORDER_FORTRAN, MPI_REAL, dtype_temp_2D_t, &
          ier )

        do m=1,number_2D
          call MPI_Get_address( array2Ds(m)%ptr(ims,jms), &
            displacement(counter), ier )
          oldtype(counter) = dtype_temp_2D_t
          counter = counter + 1
        enddo

!> create the subtype for the 3D case
!>        arraysize(1) = ime-ims+1, already set
        arraysize(2) = kme-kms+1
        arraysize(3) = jme-jms+1
!>        subarraysize(1) = ie-is+1, already set
        subarraysize(2) = ke-ks+1
        subarraysize(3) = je-js+1
!>        subarraystart(1) = is-ims, already set
        subarraystart(2) = ks-kms
        subarraystart(3) = js-jms
        call MPI_Type_create_subarray( 3, arraysize, subarraysize, &
          subarraystart, MPI_ORDER_FORTRAN, MPI_REAL, dtype_temp_3D_t, &
          ier )

        do m=1,number_3D
          call MPI_Get_address( array3Ds(m)%ptr(ims,kms,jms), &
            displacement(counter), ier )
          oldtype(counter) = dtype_temp_3D_t
          counter = counter + 1
        enddo

!> create the subtypes for the 4D cases
        subarraystart(4) = param_first_scalar-1
        do m=1,number_4D
          if (limit_4D_arrays(m) .GE. param_first_scalar) then
            arraysize(4) = limit_4D_arrays(m)
            subarraysize(4) = limit_4D_arrays(m)-param_first_scalar+1
            call MPI_Type_create_subarray( 4, arraysize, subarraysize, &
              subarraystart, MPI_ORDER_FORTRAN, MPI_REAL, &
              oldtype(counter), ier )
            call MPI_Get_address( &
              array4Ds(m)%ptr(ims,kms,jms,1), &
              displacement(counter), ier )
            counter = counter + 1
          endif 
        enddo
 
! create the all-embracing struct type
        counter = counter - 1 !> we need the actual number of arrays
        call MPI_Type_create_struct( counter, blocklength, displacement, &
          oldtype, dtype_subarray_t, ier )
        call MPI_Type_commit( dtype_subarray_t, ier )

        call MPI_Type_free( dtype_temp_2D_t, ier )
        call MPI_Type_free( dtype_temp_3D_t, ier )
        counter = number_2D + number_3D + 1
        do m=1,number_4D
          if (limit_4D_arrays(m) .GE. param_first_scalar) then
            call MPI_Type_free( oldtype(counter), ier )
            counter = counter + 1
          endif
        enddo

        deallocate( oldtype, stat=ier )
        deallocate( blocklength, stat=ier )
        deallocate( displacement, stat=ier )

        if ( myrank .EQ. 0 ) then
          call timing_record(1)
        endif

        do ii=1,inner_loop

!> =============== ping pong communication =================

!> send the data from rank 0 to rank 1 
          if ( myrank .EQ. 0 ) then
!> ==================== pack the data ======================
            pos = 0
            call MPI_Pack( MPI_BOTTOM, 1, dtype_subarray_t, &
              buffer, bytes, pos, local_communicator, ier )
            call timing_record(2)
            call MPI_Send( buffer, pos, MPI_PACKED, &
              1, itag, local_communicator, ier )
!> receive the data back from rank 1
            call MPI_Recv( buffer, bytes, MPI_PACKED, &
              1, itag, local_communicator, MPI_STATUS_IGNORE, ier )
            call timing_record(3)
!> =================== unpack the data =====================
            pos = 0
            call MPI_Unpack( buffer, bytes, pos, MPI_BOTTOM, 1, &
              dtype_subarray_t, local_communicator, ier )
            call timing_record(4)
!> now for rank 1
          else
!> receive from rank 0      
            call MPI_Recv( buffer, bytes, MPI_PACKED, &
              0, itag, local_communicator, &
              MPI_STATUS_IGNORE, ier )
!> unpack the data 
            pos = 0
            call MPI_Unpack( buffer, bytes, pos, MPI_BOTTOM, 1, &
              dtype_subarray_t, local_communicator, ier )
!> pack the data
            pos = 0
            call MPI_Pack( MPI_BOTTOM, 1, dtype_subarray_t, &
              buffer, bytes, pos, local_communicator, ier )
!> send to rank 0
            call MPI_Send( buffer, bytes, MPI_PACKED, &
              0, itag, local_communicator, &
              ier )
          endif !> of myrank .EQ. 0?

        enddo !> inner_loop

!> ======================= clean up ========================

        call MPI_Type_free( dtype_subarray_t, ier )

        deallocate( buffer, stat=ier )

        if ( myrank .EQ. 0 ) then 
          call timing_record(5)
        endif

      enddo !> outer_loop

      if ( myrank .EQ. 0 ) then
        call timing_print( .true. )
      endif

!> ======================= clean up ========================

      do m=1,number_2D
        deallocate( array2Ds(m)%ptr, stat=ier )
      enddo

      do m=1,number_3D
        deallocate( array3Ds(m)%ptr, stat=ier )
      enddo

      do m=1,number_4D
        deallocate (array4Ds(m)%ptr, stat=ier )
      enddo
     
      deallocate( array2Ds, stat=ier )
      deallocate( array3Ds, stat=ier )
      deallocate( array4Ds, stat=ier )

      end subroutine timing_wrf_sa_mpi_pack_ddt

      end module timing_wrf
