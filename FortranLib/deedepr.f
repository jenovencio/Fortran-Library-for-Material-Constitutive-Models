      subroutine deedepr(q,deltaQ,bet,dt,Betab,epr,dee)
! This subroutine evaluates dee/depr = eye(i,j)+ (ddeltaQ/depr(j))*q(i)-deltaQ*dq(i)*depr(j)
! Inputs
! q(i) -> eigen values of D i=1:3 
! deltaQ -> plastic multiplier
! epr(i) -> plastic log deformation epr= ee - deltaq*q

! Output
! dee/depr(i,j)

! Variable declaration
! Declaration of global Variables
      Double Precision q(3),deltaQ,epr(3),dee(3,3),Betab,dt,bet
!Local Variable declaration
      Double Precision ee(3), K(3), H(4)
      Double Precision dFie(3), ddFie(3), Fie(3)
      Double Precision dFip,ddFip,Fip
      Double Precision dPsi,ddPsi,Psi
      Double Precision Kzaana(6,6), Ident(3,3)
      Double Precision Lam,y(6,3),yj(6),x(6,3),xj(6)  	  
      integer i,j
	  
      Lam = 1.0d-6

      call eye(3,3,Ident)
	  
      ee=epr-deltaQ*q
      call dMaxwellElasticPotential(ee,Betab,dFie,ddFie,Fie)
      call dMaxwellPlasticPotential(Qn,deltaQ,dFip,ddFip,Fip)
      call dMaxwellViscousPotential(deltaQ,dt,dPsi,ddPsi,Psi)

! Assemble K matrix      
	  	  do i=1,3
			  K(i)=ddFie(i)*deltaQ**2+2.0d0*bet
			  H(i) = ddFie(i)*deltaQ*q(i)-dFie(i)
	      end do
		  
		  H(4)=ddFie(1)*q(1)**2+ddFie(2)*q(2)**2+ddFie(3)*q(3)**2 
     &    +ddFip+(1.0d0/dt)*ddPsi
	 
		  
		  
	      Kzaana(1,1)=K(1)
		  Kzaana(1,2)=0.0d0
		  Kzaana(1,3)=0.0d0
		  Kzaana(1,4)=H(1)
		  Kzaana(1,5)=1.0d0
		  Kzaana(1,6)=2.0d0*q(1)
		  
		  Kzaana(2,1)=0.0d0
		  Kzaana(2,2)=K(2)
		  Kzaana(2,3)=0.0d0
		  Kzaana(2,4)=H(2)
		  Kzaana(2,5)=1.0d0
		  Kzaana(2,6)=2.0d0*q(2)
		  
		  Kzaana(3,1)=0.0d0
		  Kzaana(3,2)=0.0d0
		  Kzaana(3,3)=K(3)
		  Kzaana(3,4)=H(3)
		  Kzaana(3,5)=1.0d0
		  Kzaana(3,6)=2.0d0*q(3)
		  
		  Kzaana(4,1)=H(1)
		  Kzaana(4,2)=H(2)
		  Kzaana(4,3)=H(3)
		  Kzaana(4,4)=H(4)
		  Kzaana(4,5)=0.0d0
		  Kzaana(4,6)=0.0d0
		  
		  Kzaana(5,1)=1.0d0
		  Kzaana(5,2)=1.0d0
		  Kzaana(5,3)=1.0d0
		  Kzaana(5,4)=0.0d0
		  Kzaana(5,5)=0.0d0
		  Kzaana(5,6)=0.0d0
		  
		  Kzaana(6,1)=2.0d0*q(1)
		  Kzaana(6,2)=2.0d0*q(2)
		  Kzaana(6,3)=2.0d0*q(3)
		  Kzaana(6,4)=0.0d0
		  Kzaana(6,5)=0.0d0
		  Kzaana(6,6)=0.0d0
!--------------------------------------------------------------------
	      do j=1,3
		     y(1,j)=ddFie(1)*Ident(1,j)*deltaQ
	         y(2,j)=ddFie(2)*Ident(2,j)*deltaQ
		     y(3,j)=ddFie(3)*Ident(3,j)*deltaQ
		     y(4,j)=ddFie(j)*q(j)
			 y(5,j)=0.0d0
             y(6,j)=0.0d0
             yj=y(:,j)
             call qr_solve(6,6,Kzaana,yj,xj)
             x(:,j)=xj			 
           end do			 

       do i=1,3
        	do j=1,3
				dee(i,j)=Ident(i,j)-q(i)*(x(4,j))-deltaQ*(x(i,j))
        	end do
       end do
		
	  
      end subroutine
!end code	  