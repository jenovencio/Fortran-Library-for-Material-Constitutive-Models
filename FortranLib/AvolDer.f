      subroutine AvolDer(cn1v,Jn1,dU,d2U,dUdc,d2Udc)
! the aim of this routine is to Assemble Volumetric Derivatives
! dJ/dc , d2J/dc(i)dc(j)
!the derivatives are based on chain rule
! dJ/dc(i)	= (dU/dJ)*(dJ/dc(i))
! d2U/dc(j)dc(i)=(d2U/dJ2)*(dJ/dc(j))*(dJ/dc(i))+dU/dJ*dJ/dc(i)dc(j)  
! global varibles      
      double precision cn1v(3),Jn1,dU,d2U,dUdc(3),d2Udc(3,3)
	  
	  
!local variables	  
      double precision d2Jdc(3,3),dJdc(3),l(3),dJdJ(3,3)
      integer i,j
      
      do i=1,3
          dJdc(i)=0.5d0*Jn1*(1.0d0/cn1v(i))
		  l(i)=sqrt(cn1v(i))
      end do	 
      
      
      d2Jdc(1,1)=-(l(2)*l(3))/(4.0d0*cn1v(1)**(3.0d0/2.0d0))			
      d2Jdc(1,1)=-(l(1)*l(3))/(4.0d0*cn1v(2)**(3.0d0/2.0d0))			
      d2Jdc(1,1)=-(l(2)*l(1))/(4.0d0*cn1v(3)**(3.0d0/2.0d0))	

      d2Jdc(1,2)=l(3)/(4.0d0*l(1)*l(2))
      d2Jdc(2,3)=l(1)/(4.0d0*l(2)*l(3))
      d2Jdc(1,3)=l(2)/(4.0d0*l(1)*l(3))	  
	  
      d2Jdc(2,1)=d2Jdc(1,2)
      d2Jdc(3,2)=d2Jdc(2,3)
      d2Jdc(3,1)=d2Jdc(1,3)
	  
      dUdc=dU*dJdc
      do i=1,3
	     do j=1,3
	       dJdJ(i,j)=dJdc(i)*dJdc(j)
         end do
      end do
	  
	  
      d2Udc=d2U*dJdJ+dU*d2Jdc
	 
      end subroutine
!end code	  