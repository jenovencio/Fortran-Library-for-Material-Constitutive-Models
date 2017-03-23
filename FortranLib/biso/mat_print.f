	  subroutine mat_print(A,n)
	   DOUBLE PRECISION A(3,3)
	   integer i,k
	   CHARACTER(LEN=*) n
	    print *, n
		do k=1,3
		     print *, '|', A(k,1), ',', A(k,2),',', A(k,3), '|'
		end do
	  return	
      end subroutine	
!end code        
