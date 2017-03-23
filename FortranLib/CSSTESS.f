      subroutine CSSTESS(D,A)
	  double precision D(6,6),A(6,6)
	  integer i,j
	  
	  do i=3,6
	     do j=1,6
	  A(i,j)=2.0d0(D(i,j)
	     end do
	  end do	 
	  end subroutine
!end code	  