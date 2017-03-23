      subroutine hookeLaw(e,G,K,D,s)
! global variables
      double precision e(6),G,K,s(6)
	  
! local variables	  
      double precision D(6,6),lampda
      integer i,j
!	  
      lampda=K-(2.0d0/3.0d0)*G
!
      call zeros(6,6,D)
!
!	  
      do i=1,3
	     do j=1,3
	       D(i,j)=lampda
         end do
		    D(i+3,i+3)=G
            D(i,i)=D(i,i)+2.0d0*G
      end do		 
!
      s=matmul(D,e)	  
      return
      end subroutine
!end code