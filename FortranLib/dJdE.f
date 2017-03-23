      subroutine dJdE(J,dJ,d2JdE2)
! This routine calulate de derivative of J in relation to Logarithm Strain
! J=exp(e1)*exp(e2)*exp(e3)
      double precision dJ(3),d2JdE2(3,3),J
      integer i,q
	  
      do i=1,3
	     dJ(i)=J
         do q=1,3
		     d2JdE2(i,q)=J
          end do
      end do
      return 
      end subroutine
!end code	  