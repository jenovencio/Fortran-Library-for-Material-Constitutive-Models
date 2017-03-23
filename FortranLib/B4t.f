      subroutine B4t(B,B66)
      double precision B(3,3), B66(6,6),Ident(3,3)
      double precision Bt(3,3,3,3)
      integer i,j,k,l
	  
      call eye(3,3,Ident)
      do i=1,3
	     do j=1,3
		    do k=1,3
                do l=1,3
                    Bt(i,j,k,l)=Ident(i,k)*B(j,l)+Ident(j,k)*B(i,l)
                end do
            end do
          end do
      end do		  
	  
      call tensor2matrix(Bt,B66)
      end subroutine
!end code	  