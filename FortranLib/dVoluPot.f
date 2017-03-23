      subroutine dVoluPot(Jn1,dU,d2U)
! Global Variable Declaration
      DOUBLE PRECISION Jn1, dU, d2U(3)
! Local Variable Declaration	  
      DOUBLE PRECISION k0, MP(16)
      integer tp, i
      common /MatBl/ MP
! This subroutine evaluate volumetric potential and volumetric derivatives (dU and d2U)	  
	  
	  k0=MP(15)
   
	  tp=2
      if ( tp .eq. 1) then
	      dU=k0*(log(Jn1)/Jn1)
      else
! Volumetric potential -> U=0.5*k0(J-1)^2	  
   	      dU=k0*(Jn1-1.0d0)
! d2U dependeds on the variable we are chosen to derive (
! in this case d(dU)/dE de derivative of dU in relation to E (logarithm Strain)
		  d2U(1)=k0
		  d2U(2)=k0
		  d2U(3)=k0
      end if
      end subroutine
!end code	  