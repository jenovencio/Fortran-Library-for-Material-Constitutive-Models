       subroutine DSDE(J,D,S,A)
! This subroutine computes de derivative of Cauchu Stress (S) in 
! relation to Logarithm Strain (E)
! Input are
! D - derivative od Kirchhofff Stress in relation to logarithm Strain
! J - determinant of F
! S - Cauchy Stress

!  output
! A - DS/DE

! 	A = (1/J)*D - S(x)I
! Global variables declaration
      double precision J,D(6,6),S(3,3),A(6,6)
!local variable declaration
      double precision SI(6,6), I(3,3)

      call eye(3,3,I)
      call KroIsoProd(S,I,SI) 	  
	  
      A=(1.0d0/J)*D + SI
   
	   
       end subroutine
!end code	   