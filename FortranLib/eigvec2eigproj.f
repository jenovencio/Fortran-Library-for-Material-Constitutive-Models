      subroutine eigvec2eigproj(e,M)
! This subroutine create eigen projections matrix using eigen vectors
! e - Matrix 3 by 3 of eigen vector 
! M(k) - 3 by 3 symmetric matrix eigen projection k, k=1,3
! M - 3 by 3 by 3 array 
	  DOUBLE PRECISION e(3,3),M(3,3,3),et(3,3)
	  integer i,j
	  do k=1,3
		  do i=1,3
			do j=1,3
				M(k,i,j)=e(i,k)*e(j,k)
			end do
		  end do
	  end do		
      end subroutine
!end code  
