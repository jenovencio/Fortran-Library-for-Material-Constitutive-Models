      subroutine Ad2Phidcpr(cpr,Betab,q,deltaQ,bet,dt,dPhidcpr,
     &	  d2Phidcpr)
! The aim of this routine is to Assemble Second Derivative of phi in relation to cpr
! based on chain rule:
! d2Phi(i)/dcpr(i)dcpr(j) = (d2Phi/de(i)de(j))*(de/depr)*(1/4cpr(j),cpr(i)) - dPhi/de(j)*1/2*cpr^2


!global variables
      Double Precision q(3),deltaQ,Betab,dt,cpr(3),d2Phidcpr(3,3),
     & dPhidcpr(3) 
!local variables	  
      Double Precision dee(3,3),Ident(3,3),epr(3),ddFie(3),dFie(3),Fie,
     & s
      integer i,j
	  
      do i=1,3
	     epr(i)=0.50d0*log(cpr(i))
      end do
	  
      call eye(3,3,Ident)
      call dMaxwellElasticPotential(epr,Betab,dFie,ddFie,Fie)
      
	  s=q(1)**2+q(2)**2+q(2)**2+deltaQ**2+bet**2
      if (s .gt. 0.0d0) then
          call deedepr(q,deltaQ,bet,dt,Betab,epr,dee)   
      else
          dee=Ident
      end if 	   

      do i=1,3
        dPhidcpr(i)=0.5d0*dFie(i)/(cpr(i)) 
        do j=1,3
            d2Phidcpr(i,j)=ddFie(i)*dee(i,j)/(4.0d0*cpr(j)*cpr(i))
     & -0.5d0*Ident(i,j)*dFie(i)/(cpr(i)*cpr(i))
         end do
      end do
	  
      end subroutine
!end code	