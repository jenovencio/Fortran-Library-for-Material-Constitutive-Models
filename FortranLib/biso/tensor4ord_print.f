	  subroutine tensor4ord_print(D,n)
	  DOUBLE PRECISION D(6,6)
	  integer i,k
	  CHARACTER(LEN=*) n
	    print *, n
		do k=1,6
		     print *, '|', D(1,k),D(2,k),D(3,k),D(4,k),D(5,k),D(6,k), '|'
		end do
      end subroutine
!end code        
