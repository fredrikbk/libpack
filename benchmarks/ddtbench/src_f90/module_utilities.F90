      module utilities

      contains

!> gives list_dim unique numbers in index list back (from the range of
!> 1:global_dim

      subroutine random_array_shuffle( index_list, list_dim, &
        global_dim )

      implicit none

      integer, intent(in) :: list_dim, global_dim
      integer, dimension(list_dim), intent(out) :: index_list
!> some local variables
      integer i, temp, irandom
      integer, dimension(global_dim) :: shuffle_array
      real random

      do i=1,global_dim
        shuffle_array(i) = i
      enddo

      do i=1,global_dim
        call random_number(random)
        irandom = int(random * global_dim) + 1

        temp = shuffle_array(i)
        shuffle_array(i) = shuffle_array(irandom)
        shuffle_array(irandom) = temp
      enddo

      index_list(1:list_dim) = shuffle_array(1:list_dim)

      end subroutine random_array_shuffle

!> write line into the file named by filehandle

      subroutine write_mpi( line, filehandle)

      implicit none

      include 'mpif.h'

      character(256), intent(in) :: line
      integer, intent(in) :: filehandle

      integer stringlength, ier

      stringlength = len(trim(line))

      CALL MPI_File_write( filehandle, line, stringlength, &
        MPI_CHARACTER, MPI_STATUS_IGNORE, ier )
      if ( ier .NE. MPI_SUCCESS ) then
        write (*,*) "Error at write with MPI I/O"
      endif

      end subroutine write_mpi

!> init the pseudo random number generator with the clocktime
!> copied from
!> http://gcc.gnu.org/onlinedocs/gfortran/RANDOM_005fSEED.html

      SUBROUTINE init_random_seed()

      implicit none
      
      INTEGER :: i, n, clock
      INTEGER, DIMENSION(:), ALLOCATABLE :: seed
          
      CALL RANDOM_SEED(size = n)
      ALLOCATE(seed(n))
          
      CALL SYSTEM_CLOCK(COUNT=clock)

      do i=1,n        
        seed = clock + 37 * (i - 1)
      enddo

      CALL RANDOM_SEED(PUT = seed)
          
      DEALLOCATE(seed)

      END SUBROUTINE init_random_seed

!> the following subroutines initialize the elements of an (1D/2D/3D/4D) array
!> with a unique number beginning from base
!> for each needed datatype there is a extra subroutine

      subroutine fill_unique_array_1D_real( array, DIM1, base ) 

      implicit none

      integer, intent(in) :: DIM1
      integer, intent(in) :: base
      real, dimension(DIM1), intent(out) :: array

!> local variables
      integer :: i

      do i=1,DIM1
        array(i) = base + i
      enddo

      end subroutine fill_unique_array_1D_real

      subroutine fill_unique_array_2D_real( array, DIM1, DIM2, base) 

      implicit none

      integer, intent(in) :: DIM1, DIM2
      integer, intent(in) :: base
      real, dimension(DIM1, DIM2), intent(out) :: array

!> local variables
      integer :: i, j
      integer :: counter

      counter = base
      do j=1,DIM2
        do i=1,DIM1
          counter = counter + 1
          array(i,j) = counter
        enddo
      enddo

      end subroutine fill_unique_array_2D_real

      subroutine fill_unique_array_3D_real( array, DIM1, DIM2, DIM3, &
        base) 

      implicit none

      integer, intent(in) :: DIM1, DIM2, DIM3
      integer, intent(in) :: base
      real, dimension(DIM1, DIM2, DIM3), intent(out) :: &
        array

!> local variables
      integer :: i, j, k
      integer :: counter

      counter = base
      do k=1,DIM3
        do j=1,DIM2
          do i=1,DIM1
            array(i,j,k) = counter
            counter = counter + 1
          enddo
        enddo
      enddo

      end subroutine fill_unique_array_3D_real

      subroutine fill_unique_array_4D_real( array, DIM1, DIM2, DIM3, &
        DIM4, base) 

      implicit none

      integer, intent(in) :: DIM1, DIM2, DIM3, DIM4
      integer, intent(in) :: base
      real, dimension(DIM1, DIM2, DIM3, DIM4), intent(out) :: &
        array

!> local variables
      integer :: i, j, k, l
      integer :: counter

      counter = base
      do l=1,DIM4
        do k=1,DIM3
          do j=1,DIM2
            do i=1,DIM1
              array(i,j,k,l) = counter
              counter = counter + 1
            enddo
          enddo
        enddo
      enddo

      end subroutine fill_unique_array_4D_real

      subroutine fill_unique_array_1D_double( array, DIM1, base ) 

      implicit none

      integer, intent(in) :: DIM1
      integer, intent(in) :: base
      double precision, dimension(DIM1), intent(out) :: array

!> local variables
      integer :: i

      do i=1,DIM1
        array(i) = base + i
      enddo

      end subroutine fill_unique_array_1D_double

      subroutine fill_unique_array_2D_double( array, DIM1, DIM2, base) 

      implicit none

      integer, intent(in) :: DIM1, DIM2
      integer, intent(in) :: base
      double precision, dimension(DIM1, DIM2), intent(out) :: array

!> local variables
      integer :: i, j
      integer :: counter

      counter = base
      do j=1,DIM2
        do i=1,DIM1
          counter = counter + 1
          array(i,j) = counter
        enddo
      enddo

      end subroutine fill_unique_array_2D_double

      subroutine fill_unique_array_3D_double( array, DIM1, DIM2, DIM3, &
        base) 

      implicit none

      integer, intent(in) :: DIM1, DIM2, DIM3
      integer, intent(in) :: base
      double precision, dimension(DIM1, DIM2, DIM3), intent(out) :: &
        array

!> local variables
      integer :: i, j, k
      integer :: counter

      counter = base
      do k=1,DIM3
        do j=1,DIM2
          do i=1,DIM1
            array(i,j,k) = counter
            counter = counter + 1
          enddo
        enddo
      enddo

      end subroutine fill_unique_array_3D_double

      subroutine fill_unique_array_2D_double_complex( array, DIM1, &
        DIM2, base) 

      implicit none

      integer, intent(in) :: DIM1, DIM2
      integer, intent(in) :: base
      double complex, dimension(DIM1, DIM2), intent(out) :: &
        array

!> local variables
      integer :: i, j
      integer :: counter

      counter = base
      do j=1,DIM2
        do i=1,DIM1
          array(i,j) = cmplx(counter, counter+1)
          counter = counter + 2
        enddo
      enddo

      end subroutine fill_unique_array_2D_double_complex

      subroutine fill_unique_array_5D_complex( array, DIM1, DIM2, &
        DIM3, DIM4, DIM5, base) 

      implicit none

      integer, intent(in) :: DIM1, DIM2, DIM3, DIM4, DIM5
      integer, intent(in) :: base
      complex, dimension(DIM1, DIM2, DIM3, DIM4, DIM5), intent(out) :: &
        array

!> local variables
      integer :: i, j, k, l, m
      integer :: counter

      counter = base
      do m=1,DIM5
        do l=1,DIM4
          do k=1,DIM3
            do j=1,DIM2
              do i=1,DIM1
                array(i,j,k,l,m) = cmplx(counter, counter+1)
                counter = counter + 2
              enddo
            enddo
          enddo
        enddo
      enddo

      end subroutine fill_unique_array_5D_complex

      end module utilities
