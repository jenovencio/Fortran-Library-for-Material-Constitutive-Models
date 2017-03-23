	  subroutine cauchyS(T,J,S)
! This subroutine evaluate Cauchy Stress based on:
!T - 3 by 3 matrix with Kirshhoff Stress
!J - determinant of F (deformation gradient) 
!S - vector 6 by 1 of Cauchy stress
! This subroutine needs of:
! mat2vec.f90
	  DOUBLE PRECISION T(3,3), J, S(6)
	  DOUBLE PRECISION CS(3,3)
	  CS=(1.0d0/J)*T
	  call mat2vec(CS,S)
	  end subroutine
!end code
