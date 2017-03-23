	 subroutine DP(Y0,H,Qp0,m,S0)
! Global Variable Declararion
	 Double Precision Y0,H,Qp0,m,S0
	 
!Constante Material
		Y0=0.0d0

!Modulo de encruamento [N/mm2]
		H=30.0d0  

!Constante do material
		Qp0=0.1d0

!Const Material
		m=0.8d0

!Tensao de escoamento [N/mm2]
		S0=200.0d0 
	 end subroutine
!end code	 