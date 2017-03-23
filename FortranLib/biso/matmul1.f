      SUBROUTINE matmul1(F,A,C)
	  DOUBLE PRECISION F(3,3), A(3,3),C(3,3)
      integer i,j,k
		
      do i=1,3
         do j=1,3
                 C(i,j) = 0
         end do
      end do
	  
! *** C = F'*A
      do i=1,3
         do j=1,3
		      do k=1,3
                 C(i,j) = C(i,j)+F(k,i)*A(k,j)
		      end do
         end do
      end do
      END SUBROUTINE
 !end code  
