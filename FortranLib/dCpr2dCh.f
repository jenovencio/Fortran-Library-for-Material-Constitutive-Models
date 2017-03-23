      subroutine dCpr2dCh(Cpr,dWdCpr,d2WdCpr,Fpn1,Ch,dWdCh,d2WdCh)
! the aim of this routine is to convert derivative of S in relation to
! isochoric preditec rigth Cauchy Strain (Cpr) into derivative of S in relation to
! rigth isochoric Cauchy Strain Ch
!global variable 
      double precision dWdCpr(3,3),d2WdCpr(6,6),Fpn1(3,3),dWdCh(3,3),
     & d2WdCh(6,6),Cpr(3,3),Ch(3,3)
! local variables
      double precision f(3,3),d2WdCpr4t(3,3,3,3),d2WdCh4t(3,3,3,3)
      integer i,j,k,l,A,B,C,D
      call zeros(6,6,d2WdCh)
      call matrix2tensor(d2WdCpr,d2WdCpr4t)
      call matrix2tensor(d2WdCh,d2WdCh4t)

! f= Fpn1^-1     
	  f=Fpn1
      call matinv(f,3)
!********************************************************************************
!               Ch = f(i,A)*Cpr(i,j)*f(j,B)
!********************************************************************************
          Ch=matmul((matmul(transpose(f),Cpr)),f)
!********************************************************************************	  
	  
!********************************************************************************
!                dWdCh(i,j) = f(i,A)*dWdCpr(i,j)*f(j,B)
!********************************************************************************
      dWdCh=matmul((matmul(transpose(f),dWdCpr)),f)
!********************************************************************************	  
	  
!********************************************************************************
!          d2WdCh = f(i,A)*f(j,B)*f(k,C)*f(l,D)*d2WdCpr(A,B,C,D)
!********************************************************************************	  
      do i=1,3
	     do j=1,3
		     do k=1,3
			    do l=1,3
				   do A=1,3
				      do B=1,3
					    do C=1,3
						   do D=1,3
					          d2WdCh4t(i,j,k,l)=d2WdCh4t(i,j,k,l) +
     &	f(i,A)*f(j,B)*f(k,C)*f(l,D)*d2WdCpr4t(A,B,C,D)	 
	                       end do 
						end do 
	                   end do
	                end do
	            end do
	         end do
	     end do
      end do	
	  
	  call tensor2matrix(d2WdCh4t,d2WdCh)
      end subroutine
!end code	  