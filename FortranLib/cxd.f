      subroutine cxd(S,SI)
      double precision S(3,3), SI(6,6),Ident(3,3)
      double precision St(3,3,3,3)
      integer i,j,k,l
	  
      call eye(3,3,Ident)
      do i=1,3
	     do j=1,3
		    do k=1,3
                do l=1,3
                    St(i,j,k,l)=S(i,l)*Ident(j,k)
                end do
            end do
          end do
      end do		  
	  
      call tensor2matrix(St,SI)
      end subroutine
!end code	  