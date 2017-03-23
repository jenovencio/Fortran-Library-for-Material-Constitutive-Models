      subroutine devvol(T,Tdev,Tvol)
! Calculate deviatoric and volumetric tensors of a matrix 3 by 3
! ----------------------------------------------------------------------------
! Parameters:
!   T: The input matrix 3 by 3
!   Tdev: Storage deviatoric tensor oT T
!   Tvol: Storage deviatoric tensor oT T
! ---------------------------------------------------------------------------- 
!     .. Arguments ..
      DOUBLE PRECISION T(3,3)
	  DOUBLE PRECISION Tdev(3,3)
	  DOUBLE PRECISION Tvol(3,3)
!    

!local variables
      DOUBLE PRECISION p
     
	 
	  p = T(1,1)+T(2,2)+T(3,3)
	  p = p/3.0d0
	  Tvol(1,1)=p
	  Tvol(1,2)=0.d0
	  Tvol(1,3)=0.d0
	  Tvol(2,1)=0.d0
	  Tvol(2,2)=p
	  Tvol(2,3)=0.d0
	  Tvol(3,1)=0.d0
	  Tvol(3,2)=0.d0
	  Tvol(3,3)=p
	
	  Tdev(1,1)=T(1,1)-p
	  Tdev(1,2)=T(1,2)
	  Tdev(1,3)=T(1,3)
	  Tdev(2,1)=T(2,1)
	  Tdev(2,2)=T(2,2)-p
	  Tdev(2,3)=T(2,3)
	  Tdev(3,1)=T(3,1)
	  Tdev(3,2)=T(3,2)
	  Tdev(3,3)=T(3,3)-p
	
      end subroutine
!end code  