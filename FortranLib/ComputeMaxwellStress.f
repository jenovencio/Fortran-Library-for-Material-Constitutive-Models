	  subroutine ComputeMaxwellStress(ee,Betab,cprv,MprProj,dFiedCpr)
!     Global Variable declaration
!	  Input variables	
	  DOUBLE PRECISION ee(3),Betab,cprv(3),MprProj(3,3,3)
!	  Output variables	  
	  DOUBLE PRECISION dFiedCpr(3,3)
!     Local variable declaration
	  DOUBLE PRECISION ddFie(3),Fie(3), dFie(3)
	  integer i
! This subroutine needs:
! zeros.f	  
! dMaxwellElasticPotential.f

	  call dMaxwellElasticPotential(ee,Betab,dFie,ddFie,Fie)
	  call zeros(3,3,dFiedCpr)
	  do i=1,3
	   dFiedCpr=dFiedCpr
     & + dFie(i)*(1.0d0/(2.0d0*cprv(i)))*MprProj(i,:,:)
	  end do
	  
	  return
	  end subroutine
!end code	  