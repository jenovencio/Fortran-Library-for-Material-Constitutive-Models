      subroutine j2calc(s,N,qtrial)
 !global variables     
      double precision s(6),qtrial,N(6)
!local variables 	  
      double precision pEl,sigDev(6),qEl,TWO
	  
      TWO=2.0d0
c *** hydrostatic pressure stress
      pEl = -(1.0d0/3.0d0) * (s(1) + s(2) + s(3))
c *** compute the deviatoric stress tensor
      sigDev(1) = s(1) + pEl
      sigDev(2) = s(2) + pEl
      sigDev(3) = s(3) + pEl
      sigDev(4) = s(4)
      sigDev(5) = s(5)
      sigDev(6) = s(6)
c *** compute von-mises stress
      qEl = 
     &  sigDev(1) * sigDev(1)+sigDev(2) * sigDev(2)+
     &  sigDev(3) * sigDev(3)+
     &  TWO*(sigDev(4) * sigDev(4)+ sigDev(5) * sigDev(5)+ 
     &  sigDev(6) * sigDev(6))
      qtrial = sqrt( 1.5d0 * qEl)
c ***      
c	  call vec_print(sigDev,6,"sigDev")
c ***	  

c	  print *,"normS=",normS
	  print *,"qtrial=",qtrial
c ***	  
c ***	  
      N=1.5d0*sigDev/qtrial
c ***
      return	  
      end subroutine
!end code	  