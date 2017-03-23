      subroutine zeros(n,m,A)
! this subroutine generates n by m zeros matrix
	  integer n,m,i,j
	  DOUBLE PRECISION A(n,n)
	
	  do i=1,n
		do j=1,m
			A(i,j)=0.0d0
		end do
	  end do	
      end subroutine 
!end code  
