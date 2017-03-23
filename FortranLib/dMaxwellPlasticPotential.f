      subroutine dMaxwellPlasticPotential(Qn,x,dFip,ddFip,Fip)
! Global variable declaration
      Double Precision Qn,x,dFip,ddFip,Fip
! Local variable declaration
      Double Precision Y0,H,Qp0,m,S0, DP(5),deltaQ,TP(4)
      common /MatDP/ DP
      common /MatTP/ TP
! This routine needs:
! DP.f
!****************************************************************************************	  
	  Y0=DP(1)
	  H=DP(2)
	  Qp0=DP(3)
	  m=DP(4)
	  S0=DP(5)
	  
	  deltaQ=x
!	  call vec_print(DP,5,"DP")
	  
	  Fip=S0*(Qn+deltaQ) + 0.5d0*H*(Qn+deltaQ)*(Qn+deltaQ)  
	  dFip=S0 + H*(Qn+deltaQ)    
	  ddFip=H
      end subroutine
!end code	  