      subroutine MPI_INIT(ierr)
      integer ierr
      print*,"MPI_INIT"
      call interposer_init()
      call PMPI_INIT(ierr)
      end


      subroutine MPI_SEND(BUF, COUNT, DATATYPE, DEST, TAG, COMM, IERROR)
      INTEGER COUNT, DATATYPE, DEST, TAG, COMM, IERROR
      call interposer_pack(BUF)
      call PMPI_SEND(BUF, COUNT, DATATYPE, DEST, TAG, COMM, IERROR)
      end
