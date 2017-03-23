	  subroutine mattransp(A,n,At)
		integer i,j,n  
	    DOUBLE PRECISION At(n,n), A(n,n)

! ***	At=A'
		do i=1,n
		   do j=1,n
		      At(j,i)=A(i,j)
			end do
		end do	
	  end subroutine
!end code          
