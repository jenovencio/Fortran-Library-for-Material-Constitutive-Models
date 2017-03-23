      subroutine dMaxwellViscousPotential(deltaQ,dt,dPsi,ddPsi,Psi)
! Global variable declaration
      Double Precision deltaQ,dt,dPsi,ddPsi,Psi

! Local variable declaration
      Double Precision Y0,H,Qp0,m,S0, DV(50),TP(4)
      common /MatDV/ DV
      common /MatTP/ TP
! This routine needs:
! DP.f
!****************************************************************************************	  
!	  
	  Qp=deltaQ/dt
      if (TP(4) .eq. 1.0D0) then
		  Y0=DV(1)
		  Psi=Y0*Qp
		  dPsi=Y0
		  ddPsi=0.0d0
       elseif (TP(4) .eq. 2.0D0) then
	      Y0=DV(1)
	      Qp0=DV(2)
          m=DV(3)
		  Psi=(m/(m+1.0d0))*Qp0*Y0*(Qp/Qp0)**(m/(m+1.0d0))
		  dPsi=Qp0**(-1.0d0/m)*Y0*Qp**(1.0d0/m)
		  ddPsi=(Qp0**(-1.0d0/m))*(1.0d0/m)*(Y0*Qp**(-1+1/m))
      end if
      return	  
      end subroutine
!end code	  