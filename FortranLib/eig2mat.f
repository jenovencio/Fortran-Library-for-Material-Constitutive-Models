	  subroutine eig2mat(V,E,Vt,M)
! This subroutine evaluate M=V*E*V'
! Where
! V -> matrix (3 by 3)
! E -> matrix (3 by 3)
! V' -> matrix(3 by 3)
! M -> matrix (3 by 3)
! Global variable declaration
! Input
	  DOUBLE PRECISION V(3,3),E(3,3),Vt(3,3) 	
! Output 
	  DOUBLE PRECISION 	M(3,3)	
! local variable declaration
	  DOUBLE PRECISION EVt(3,3)
	  EVt=matmul(E,Vt)
	  M=matmul(V,EVt)
	  
	  end subroutine
!end code	  