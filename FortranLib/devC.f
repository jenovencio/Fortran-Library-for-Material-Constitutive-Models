	  subroutine DEVC(A,Y,dC)
! This subroutine evaluates deviatoric strain of C
! Where
! A-> 3 by 3 matrix
! Y-> 3 by 3 matrix
! devC-> 3 by 3 matrix
! This subroutine needs:
! matinv.f
! prodinv.f
	  
! Global variable declaration
	  DOUBLE PRECISION A(3,3), Y(3,3), dC(3,3)
	  
! local variable declaration
	  DOUBLE PRECISION  Yinv(3,3), c,t
	  t=1.0d0/3.0d0
	  call prodint(A,Y,c)
	  Yinv=Y
	  call matinv(Yinv,3)
	  dC=A-t*c*Yinv
	  
	  end subroutine
!end code	  