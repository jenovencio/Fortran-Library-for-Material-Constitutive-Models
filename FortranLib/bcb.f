      subroutine bcb(B,C,D)
! The aim of this routine is to compute B(i,m)*C4t(m,j,k,p)*B(p,l)      
	  double precision B(3,3),C(6,6),D(6,6)
      double precision C4t(3,3,3,3),D4t(3,3,3,3)
      integer i,j,k,l,m,p
      call matrix2tensor(C,C4t)
      call zeros(6,6,D)
      call matrix2tensor(D,D4t)	  

	
      do i=1,3
	     do j=1,3
		     do k=1,3
			    do l=1,3
				   do m=1,3
				      do p=1,3
					     D4t(i,j,k,l)=D4t(i,j,k,l)+
     &						 B(i,m)*C4t(m,j,k,p)*B(p,l)      
	                   end do
	                end do
	            end do
	         end do
	     end do
      end do	 
	  
      call tensor2matrix(D4t,D)
	  
      end subroutine
!end code	  