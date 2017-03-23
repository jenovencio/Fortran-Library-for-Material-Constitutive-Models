      subroutine Calcor(x,dt,r,Qn,Betab,Taub,epr,deltax)
! This subroutine calculates de correction\direction of minimization
! Declaration of global Variables
!input
      Double Precision Betab,epr(3),Taub,Qn, x(6), r(6), dt
! Output
      Double Precision deltax(6)
	  
! Declaration of local Variables
      Double Precision rr(6), q(3),ee(3), K(3), H(4)
      Double Precision dFie(3), ddFie(3), Fie(3)
      Double Precision dFip,ddFip,Fip
      Double Precision dPsi,ddPsi,Psi
      Double Precision Kzaana(6,6), deltax1(6)
      Double Precision dx(6),xf(6),xi(6),rf(6),ri(6),Kfi(6)
      Double Precision Kzao(6,6), Kzainv(6,6)
      integer c
	  
! Subroutine parameter
	  c=2
	  rr(1)=r(4)
	  rr(2)=r(1)
	  rr(3)=r(2)
	  rr(4)=r(3)
	  rr(5)=r(5)
	  rr(6)=r(6)
	  
	  
      if (c .eq. 2) then
		  q(1)=x(1)
		  q(2)=x(2)
		  q(3)=x(3)
	  
		  ee=epr-x(4)*q
	      call dMaxwellElasticPotential(ee,Betab,dFie,ddFie,Fie)
	      call dMaxwellPlasticPotential(Qn,x(4),dFip,ddFip,Fip)
	      call dMaxwellViscousPotential(x(4),dt,dPsi,ddPsi,Psi)
	      
          do i=1,3
			  K(i)=ddFie(i)*x(4)**2+2*x(6)
			  H(i) = ddFie(i)*x(4)*x(i)-dFie(i)
	      end do
		  
		  H(4)=ddFie(1)*x(1)**2+ddFie(2)*x(2)**2+ddFie(3)*x(3)**2 
     &    +ddFip+(1/dt)*ddPsi
	 
		  
		  
	      Kzaana(1,1)=H(4)
		  Kzaana(1,2)=H(1)
		  Kzaana(1,3)=H(2)
		  Kzaana(1,4)=H(3)
		  Kzaana(1,5)=0.0d0
		  Kzaana(1,6)=0.0d0
		  
		  Kzaana(2,1)=H(1)
		  Kzaana(2,2)=K(1)
		  Kzaana(2,3)=0.0d0
		  Kzaana(2,4)=0.0d0
		  Kzaana(2,5)=1.0d0
		  Kzaana(2,6)=2.0d0*q(1)
		  
		  Kzaana(3,1)=H(2)
		  Kzaana(3,2)=0.0d0
		  Kzaana(3,3)=K(2)
		  Kzaana(3,4)=0.0d0
		  Kzaana(3,5)=1.0d0
		  Kzaana(3,6)=2.0d0*q(2)
		  
		  Kzaana(4,1)=H(3)
		  Kzaana(4,2)=0.0d0
		  Kzaana(4,3)=0.0d0
		  Kzaana(4,4)=K(3)
		  Kzaana(4,5)=1.0d0
		  Kzaana(4,6)=2.0d0*q(3)
		  
		  Kzaana(5,1)=0.0d0
		  Kzaana(5,2)=1.0d0
		  Kzaana(5,3)=1.0d0
		  Kzaana(5,4)=1.0d0
		  Kzaana(5,5)=0.0d0
		  Kzaana(5,6)=0.0d0
		  
		  Kzaana(6,1)=0.0d0
		  Kzaana(6,2)=2.0d0*q(1)
		  Kzaana(6,3)=2.0d0*q(2)
		  Kzaana(6,4)=2.0d0*q(3)
		  Kzaana(6,5)=0.0d0
		  Kzaana(6,6)=0.0d0
		  
		  call qr_solve(6,6,Kzaana,-rr,deltax1)
!		  Kzainv=Kzaana
!		  call matinv(Kzainv,6)
!		  deltax1=matmul(Kzainv,-rr)

		  
		  deltax(1)=deltax1(2)
          deltax(2)=deltax1(3)
          deltax(3)=deltax1(4)
          deltax(4)=deltax1(1)
          deltax(5)=deltax1(5)
          deltax(6)=deltax1(6)
		  
	  else
! It is not implemented
	  print *, "It is not implemented"
	  end if
	  
	  return
	  end subroutine
!end code	  