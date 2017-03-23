	  subroutine DX2DX(X,DX2)
! This subroutine evaluate Derivative of X^2
! X = symmetric second order tensor 3 by 3
! DX2 = symmetric matrix 6 by 6
! DI = Identity matrix
	  DOUBLE PRECISION X(3,3),DX2(6,6)
	  DOUBLE PRECISION A(3,3,3,3)
	  DOUBLE PRECISION DI(3,3)
	  integer i,j,k,l
	  call eye(3,3,DI)
! Inicialize matrix	
	  do i=1,3
		do j=1,3
			do k=1,3
				do l=1,3
					A(i,j,l,k)=0.0d0
				end do
			end do	
		end do
	  end do
	
	  do i=1,3
		do j=1,3
			do k=1,3
				do l=1,3
					A(i,j,k,l)=DI(i,k)*X(l,j) + DI(i,l)*X(k,j) +
     &              DI(j,l)*X(i,k) + DI(k,j)*X(i,l)					
				end do
			end do	
		end do
	  end do	
	  
	  A=0.5d0*A
	
	  DX2(1,1)=A(1,1,1,1)
	  DX2(1,2)=A(1,1,2,2)
	  DX2(1,3)=A(1,1,3,3)
	  DX2(1,4)=A(1,1,1,2)
	  DX2(1,5)=A(1,1,2,3)
	  DX2(1,6)=A(1,1,1,3)
	
	  DX2(2,1)=A(2,2,1,1)
	  DX2(2,2)=A(2,2,2,2)
	  DX2(2,3)=A(2,2,3,3)
	  DX2(2,4)=A(2,2,1,2)
	  DX2(2,5)=A(2,2,2,3)
	  DX2(2,6)=A(2,2,1,3)
	  
	  DX2(3,1)=A(3,3,1,1)
	  DX2(3,2)=A(3,3,2,2)
	  DX2(3,3)=A(3,3,3,3)
	  DX2(3,4)=A(3,3,1,2)
	  DX2(3,5)=A(3,3,2,3)
	  DX2(3,6)=A(3,3,1,3)	
	  
	  DX2(4,1)=A(1,2,1,1)
	  DX2(4,2)=A(1,2,2,2)
	  DX2(4,3)=A(1,2,3,3)
	  DX2(4,4)=A(1,2,1,2)
	  DX2(4,5)=A(1,2,2,3)
	  DX2(4,6)=A(1,2,1,3)
	  
	  DX2(5,1)=A(2,3,1,1)
	  DX2(5,2)=A(2,3,2,2)
	  DX2(5,3)=A(2,3,3,3)
	  DX2(5,4)=A(2,3,1,2)
	  DX2(5,5)=A(2,3,2,3)
	  DX2(5,6)=A(2,3,1,3)
	  
	  DX2(6,1)=A(1,3,1,1)
	  DX2(6,2)=A(1,3,2,2)
	  DX2(6,3)=A(1,3,3,3)
	  DX2(6,4)=A(1,3,1,2)
	  DX2(6,5)=A(1,3,2,3)
	  DX2(6,6)=A(1,3,1,3)
	
	
	
	
	  end subroutine
!end code  
