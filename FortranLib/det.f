	  subroutine det(F,J)
! Calculate determinante of a matrix 3 by 3
! ----------------------------------------------------------------------------
! Parameters:
!   F: The input matrix 3 by 3
!   J: Storage determinante of F
! ---------------------------------------------------------------------------- 
!     .. Arguments ..
      DOUBLE PRECISION F(3,3)
	  DOUBLE PRECISION J,Jp,Jn
		Jp=F(1,1)*F(2,2)*F(3,3) 
		Jp=Jp + F(1,2)*F(2,3)*F(3,1) + F(1,3)*F(2,1)*F(3,2)
		Jn=F(1,3)*F(2,2)*F(3,1) 
		Jn=Jn + F(1,1)*F(2,3)*F(3,2) + F(2,1)*F(1,2)*F(3,3)
		J=Jp - Jn
	  end subroutine
!end code  
