      subroutine dMaxwellElasticPotential(ee,Betab,dFie,ddFie,Fie)
! This routine evaluate elastic potential derivatives
! Output
! Fie 
! dFie
! ddFie(i)
	  
! This subroutine needs:
! MPOgden4.f
	  
! global variable declaration
      Double Precision ee(3),Betab,dFie(3)
      Double Precision ddFie(3),Fie
	  
! local Variable declaration
      Double Precision MP(16), mu(3),alfm(3),TP(4) 
      integer i,j
      common /MatBl/ MP
      common /MatTP/ TP
!********************************************************************	  


	  mu(1) = MP(1)*Betab
	  mu(2) = MP(2)*Betab
	  mu(3) = MP(3)*Betab
	 
	  alfm(1) = MP(4)
	  alfm(2) = MP(5)
	  alfm(3) = MP(6)	  

      
      if (TP(1) .eq. 1.0d0) then
        do i=1,3
			Fie=mu(1)*ee(i)*ee(i)
			dFie(i)=2*mu(1)*ee(i)
			ddFie(i)=2*mu(1)
        end do
      elseif (TP(1) .eq. 2.0d0) then
		do i=1,3
			Fie=0.0d0
			dFie(i)=0.0d0
			ddFie(i)=0.0d0
            do j=1,3
			 Fie= Fie + (mu(j)/alfm(j))*((exp(ee(i)))**alfm(j)-1) 
			 dFie(i)= dFie(i) + mu(j)*(exp(ee(i)))**alfm(j)
			 ddFie(i)= ddFie(i) + mu(j)*(exp(ee(i)))**alfm(j)*alfm(j) 
            end do
         end do
      end if

      end subroutine
!end code	  