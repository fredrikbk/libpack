      module datatypes

      type arrayPointer2D
        real, pointer :: ptr(:,:)
      end type arrayPointer2D

      type arrayPointer3D
        real, pointer :: ptr(:,:,:)
      end type arrayPointer3D

      type arrayPointer4D
        real, pointer :: ptr(:,:,:,:)
      end type arrayPointer4D      

      end module datatypes
