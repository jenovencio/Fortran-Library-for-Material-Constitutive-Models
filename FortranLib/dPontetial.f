	  subroutine dPontetial(q,deltaQ,dt,Betab,epr,cprv
     &	  ,dFiecpr,ddFiecpr)
! This subroutine evaluates the first (dFiedcpr) and second (ddFiedcpr)
! derivatives of the Elastic Potential Fie in relation to Cpr = (F^-T)*C*(F^-1)

! This subroutine needs:
! deedepr(q,deltaQ,dt,Betab,epr,dee)
! dMaxwellElasticPotential
!***********************************************************************
! Variable declaration
! Declaration of global Variables
	  DOUBLE PRECISION epr(3),Betab,cprv(3)
	  Double Precision q(3),deltaQ,dt
	  DOUBLE PRECISION dFiecpr(3),ddFiecpr(3,3) 
!     Local variable declaration
	  DOUBLE PRECISION ddFie(3),Fie(3), dFie(3), dee(3,3), ee(3), cc
     & cpr2	  
	  integer i, j

       ee=epr-deltaQ*q
	   call dMaxwellElasticPotential(ee,Betab,dFie,ddFie,Fie)
       call deedepr(q,deltaQ,dt,Betab,epr,dee)
	   
	   do j=1,3
		   dFiecpr(j)= dFie(j)*(1.0d0/(2.0d0*cprv(j)))
		   do i=1,3
		       cc=4.0d0*cprv(j)*cprv(i)
			   cpr2=2*(cprv(j))**2.0d0
		       ddFiecpr(j,i)=ddFie(j)*dee(j,i)*(1.0d0/cc)
     &         -dFie(j)*(1/cpr2)
		   end do
	   end do
	  
      return	  
	  end subroutine
!end code	  