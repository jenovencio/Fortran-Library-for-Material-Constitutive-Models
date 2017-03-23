      subroutine dKdee2dKdb(dtaldee,ee,dtaldb)
! The aim of this routine is to convert derivative of Kirshhoff Stress 
! in relation to logatithm principal strain (dtaldee /dKdee) to
! the derivative of Kirshhoff Stress in relation to de principal values of
! left Cauchy Strain (B=V*bi*V)
! Where B = F*F'   	  
! ee - is the principal logarithm strain ee=0.5*ln(bi)
! dtaldee - derivative of Kirshhoff Stress in relation to logatithm principal strain 
! dtaldb - derivative of Kirshhoff Stress in relation to de principal values of left Cauchy Strain
! This rotuine applies the chain rule as follow:
! dtaldb(i,j) = dtal(i)/dee(k)*dee(k)/db(i) 
!  Global variables
      double precision dtaldee(3,3),dtaldb(3,3),ee(3)
! Local Variables
      double precision deedb(3,3)
      integer i
      call zeros(3,3,deedb)

! derivative of ee(k) in relation to b(i) -> ee(i)=0.5*ln(b(i)) 	  
      do i=1,3
	     deedb(i,i)=0.50d0/exp(2.0d0*ee(i))
      end do
	  
      dtaldb=matmul(dtaldee,deedb)
	  
	  
      end subroutine
!end code	  