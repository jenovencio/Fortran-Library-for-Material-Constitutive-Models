	  subroutine vec_print(A,m,n)
	    integer i,k,m
		DOUBLE PRECISION A(m)
	    CHARACTER(LEN=*) n
	    print *, n
		     do i=1,m
				print *, '|', A(i), '|'
			 end do
      end subroutine	
!end code        
