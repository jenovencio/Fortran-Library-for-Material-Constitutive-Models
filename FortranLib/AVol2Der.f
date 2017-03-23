		subroutine AVol2Der(d2U,dJ,dU,d2JdE2,dvol)
! The aim of this routine is to Assemble Volumetric 2nd Derivative
! Second derivative are assemble based of logarithm principal strain
! The formula to assemble volumetric part is based on chain rule
! dvol = d(dU)/dE*dJ + dU*d2JdE2

!Global Variable declaration
      double precision d2U(3),dJ(3),dU,d2JdE2(3,3),dvol(3,3)
	  
! local varible declaration 	  
      double precision d2UdJ(3,3)
      integer i,j
	  
      do i=1,3
	     do j=1,3
	         d2UdJ(i,j)=d2U(j)*dJ(3)
         end do
      end do		 
	  
	  dvol=d2UdJ+dU*d2JdE2		 
      end subroutine
!end code	  