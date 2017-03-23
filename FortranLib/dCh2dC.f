      subroutine dCh2dC(dCh2,Fpn1,dC)
! the aim of this routine is to convert derivative of S in relation to
! isochoric rigth Cauchy Strain (Ch) into derivative of S in relation to
! rigth Cauchy Strain C
!global variable 
      double precision dCh2(6,6),Fpn1(3,3),dC(6,6)
! local variables
      double precision f(3,3),dCh24t(3,3,3,3),dC4t(3,3,3,3)
      integer i,j,k,l,A,B,C,D
      call zeros(6,6,dC)
      call matrix2tensor(dCh2,dCh24t)
      call matrix2tensor(dC,dC4t)

! f= Fpn1^-1     
	  f=Fpn1
      call matinv(f,3)
	  
      do i=1,3
	     do j=1,3
		     do k=1,3
			    do l=1,3
				   do A=1,3
				      do B=1,3
					    do C=1,3
						   do D=1,3
					          dC4t(i,j,k,l)=dC4t(i,j,k,l) +
     &	f(i,A)*f(j,B)*f(k,C)*f(l,D)*dCh24t(A,B,C,D)	 
	                       end do 
						end do 
	                   end do
	                end do
	            end do
	         end do
	     end do
      end do	
	  
      end subroutine
!end code	  