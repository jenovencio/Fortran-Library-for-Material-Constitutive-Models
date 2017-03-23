      subroutine ADev2Der(dE,d2PdEi,ddev)
! The aim of this routine is to Assemble Deviatoric 2nd Derivative
! Second derivative are assemble based of logarithm principal strain
! The formula to assemble deviatoric part is based on chain rule
! ddev = dEdev/dE*(d2P/dEdev2)*dEdev/dE + (dP/dEdev*d(dEdev/dE)/dEdev)*dEdev/dE
!However d(dEdev/dE)/dEdev)=0
!because dEdev/dE = K (constant matrix)

!Global Variable declaration
      double precision dE(3,3),d2PdEi(3,3),ddev(3,3)
	  
      call eig2mat(dE,d2PdEi,dE,ddev)
	  
      return
      end subroutine
!end code