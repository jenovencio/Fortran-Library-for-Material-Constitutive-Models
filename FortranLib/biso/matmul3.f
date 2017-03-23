       SUBROUTINE matmul3(F,A,B)	   
	     DOUBLE PRECISION F(3,3), A(3,3), B(3,3)
        integer i,j,k
		
	  do i=1,3
         do j=1,3
                 B(i,j) = 0
         end do
      end do
		
! *** M = F*A
      do i=1,3
         do j=1,3
		      do k=1,3
                 B(i,j) = B(i,j)+F(i,k)*A(k,j)
		      end do
         end do
      end do
      END SUBROUTINE
!end code  
