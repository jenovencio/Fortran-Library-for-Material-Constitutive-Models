      subroutine A2DBP(edev,J,Betab,d)
! The aim of this routine is to Assemble Deviatoric + Volumetric  2nd Derivative (e1,e2,e3)	  
! where e1,e2,e3 are the principal values of logarithm strain
! edev(i) = e(i)-(1/3)*(e1+e2+e3) 
! d = ddev + dvol
! global varible declaration
      double precision edev(3),J,d(3,3),Betab

! local variable declaration	  
      double precision d2U(3),dJ(3),dU,d2JdE2(3,3),dvol(3,3)
      double precision dE(3,3),d2PdEi(3,3),ddev(3,3)
      double precision dFie(3),ddFie(3),Fie
      integer i 
	  
      call dVoluPot(J,dU,d2U)
      call dJdE(J,dJ,d2JdE2)
      call dEdevdE(dE)
      call dMaxwellElasticPotential(edev,Betab,dFie,ddFie,Fie)
      call zeros(3,3,d2PdEi)
      do i=1,3
	     d2PdEi(i,i)=ddFie(i)
      end do
	  
! printing to debug	  
!      print *, "J=",J
!      print *, "dU=",dU
!	  call vec_print(d2U,3,"d2U")
!	  call vec_print(dJ,3,"dJ")
!	  call mat_print(d2JdE2,"d2JdE2")
!	  call mat_print(dE,"dE")
!	  call mat_print(d2PdEi,"d2PdEi")
	  
!---------------------------------------------------------

      call AVol2Der(d2U,dJ,dU,d2JdE2,dvol)
      call ADev2Der(dE,d2PdEi,ddev)
	  
      d=ddev+dvol
!	  call mat_print(ddev,"ddev")
!	  call mat_print(dvol,"dvol")
!	  call mat_print(d,"d")
      end subroutine
!end code	  