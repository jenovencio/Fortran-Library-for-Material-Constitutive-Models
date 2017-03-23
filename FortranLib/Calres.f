	  subroutine Calres(x,Qn,Betab,Taub,epr,dt,r)
	! This subroutine calculates the residual of the minimization
	! Global Variable declaration
	  Double Precision Betab,epr(3),Taub,Qn, x(6), r(6), dt
	! Local variable declaration
	  Double Precision q(3), ee(3), dFip, ddFip, Fip
	  Double Precision dFie(3),ddFie(3), Fie(3), dPsi,ddPsi, Psi
	  
	! This subroutine needs:
    ! dMaxwellPlasticPotential(Qn,x(4),dFip,ddFip,Fip)  -> dMaxwellPlasticPotential.f
    ! dMaxwellElasticPotential(ee,Betab,dFie,ddFie,Fie) -> dMaxwellElasticPotential.f
    !	dMaxwellViscousPotential(x(4),dt,dPsi,ddPsi,Psi)  -> dMaxwellViscousPotential.f
!****************************************************************************************

	  q(1)=x(1)
	  q(2)=x(2)
	  q(3)=x(3)
	  ee=epr - q*x(4)
	  
      call dMaxwellElasticPotential(ee,Betab,dFie,ddFie,Fie)
	  call dMaxwellPlasticPotential(Qn,x(4),dFip,ddFip,Fip)
	  call dMaxwellViscousPotential(x(4),dt,dPsi,ddPsi,Psi)
	  
!	  print *,"dFip=",dFip 
!	  print *,"dPsi=",dPsi 
	  
!	  %Calculo do residuo
	  r(1) = -dFie(1)*x(4) + x(5) +2.0d0*x(6)*x(1)
	  r(2) = -dFie(2)*x(4) + x(5) +2.0d0*x(6)*x(2)
	  r(3) = -dFie(3)*x(4) + x(5) +2.0d0*x(6)*x(3)
	  r(4) = -dFie(1)*x(1)-dFie(2)*x(2)-dFie(3)*x(3)+dFip+dPsi
	  r(5)= x(1)+x(2)+x(3)
	  r(6)= x(1)**2 + x(2)**2 + x(3)**2 -3.0d0/2.0d0
	  
	  end subroutine
!end code	  