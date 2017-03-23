      subroutine dVolu(Jn1,dU)
! Global Variable Declaration
      DOUBLE PRECISION Jn1, dU
! Local Variable Declaration	  
      DOUBLE PRECISION k0,TP(4),MP(16)
      common /MatBl/ MP
      common /MatTP/ TP
	  k0=MP(15)
   
      if ( TP(2) .eq. 1.0D0) then
	      dU=k0*(log(Jn1)/Jn1)
      elseif ( TP(2) .eq. 2.0D0) then
   	      dU=k0*(Jn1-1.0d0)
      end if
      end subroutine
!end code	  