      subroutine AKStress(dFie,J,ks)
! The aim of this routine is to Assemble Kirchhoff principal Stress (ks)
! and return Cauchy Stress

! global varible declaration
      double precision dFie(3),J,ks(3)

! local variable declaration	  
      double precision dE(3,3),dPhidE(3),dU
      double precision d2U(3),dJ(3),d2JdE2(3,3)
      integer i,k
	  
      call dVoluPot(J,dU,d2U)
      call dJdE(J,dJ,d2JdE2)
      call dEdevdE(dE)
      do i=1,3
	     dPhidE(i)=0.0d0
	     do k=1,3
		     dPhidE(i)=dPhidE(i)+dFie(k)*dE(k,i)
         end do
      end do	 
	  
	  ks=dPhidE+dU*dJ
	  
      return
      end subroutine
!end code	  