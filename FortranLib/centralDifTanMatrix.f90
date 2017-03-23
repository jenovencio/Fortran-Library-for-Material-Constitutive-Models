program teste

    DOUBLE PRECISION  prop(10),stress(6)
    DOUBLE PRECISION  dsdePl(6,6)
	DOUBLE PRECISION  defGrad_t(3,3)
    DOUBLE PRECISION  B(3,3), V(3,3), l(3), t(3)   
    DOUBLE PRECISION  tB(3,3), KS(3,3), J3, D(6,6)
	DOUBLE PRECISION  lampda, etrial(3,3),alfan(10)
	DOUBLE PRECISION  alfan1(10),CS(6),e(3), poisson
	DOUBLE PRECISION  M1(3,3),M2(3,3), YM , deform
	DOUBLE PRECISION  Strain(6),t11, deform12, DI(3,3)
	DOUBLE PRECISION  Strial(6), JF2, JF, sl(3), c1
	DOUBLE PRECISION  CDM(6,6),bi(6),bcp(6),bcn(6),p
	DOUBLE PRECISION  Bp(3,3),lp(3),Vp(3,3),KSvp(6)
	DOUBLE PRECISION  Bn(3,3),ln(3),Vn(3,3),KSvn(6),tn(3), KSn(3,3)
	DOUBLE PRECISION  DP(6),Da(6,6), VL(3,3),VLVt(3,3),LM(3,3),LV(6)
	DOUBLE PRECISION  tv(6),tm(3,3),vtV(3,3)
	call eye(3,3,DI)
	call eye4sym(CDM)
	call zeros(6,6,D)

		t11=10.0d0
	poisson=0.22d0
	YM=1.0d1
	deform=t11/YM
	deform12=-poisson*t11/YM
	
	prop(1)=YM/(2.0d0*(1.0d0 + poisson))
	print *, "G=", prop(1)
	prop(2)=YM/(3.0d0*(1.0d0 - 2.0d0*poisson))
	print *, "K=", prop(2)
	
	Strain(1)=deform
	Strain(2)=deform12
	Strain(3)=deform12
	Strain(4)=0.0d0
	Strain(5)=0.0d0
	Strain(6)=0.0d0
	call vec2mat(Strain,etrial)
	! relative perturbation
	p=1.0d-9
	c1=1.0d0/(2*p)
	B=DI+etrial
	
	!call tensor4ord_print(CDM,"Central Difference Matrix ")
	
	call mat2vec(B,bi)
	do i=1,6
		print *, "iteration ", i
		bcp = bi + p*CDM(i,:)
		bcn = bi - p*CDM(i,:)
    
		call vec2mat(bcp,Bp)
		call vec2mat(bcn,Bn)
	   !call mat_print(Bp,"Bp")
	! Forward difference
		call mat_print(Bp,"Bp")
		call eigsym33(Bp,Vp,lp)
		!call mat_print(Vp,"V")
		!call vec_print(lp,3,"bp")
		LV(1)=lp(1)
		LV(2)=lp(2)
		LV(3)=lp(3)
		LV(4)=0.0D0
		LV(5)=0.0D0
		LV(6)=0.0D0
		call vec2mat(LV,LM)
		!call mat_print(LM,"LM")
		call matmul3(Vp,LM,VL)
		call matmul2(VL,Vp,VLVt)
		!call mat_print(VLVt,"VLVt")
		call odgen(lp,1,1.0d0,2.0d0,1.0d0,t,tB)
	! Kirchhoff Stress
		!call vec_print(t,3,"t")
		tv(1)=t(1)*lp(1)*2.0d0
		tv(2)=t(2)*lp(2)*2.0d0
		tv(3)=t(3)*lp(3)*2.0d0
		tv(4)=0.0d0
		tv(5)=0.0d0
		tv(6)=0.0d0
		call vec2mat(tv,tm)
		!call mat_print(tm,"tm")
		call matmul3(Vp,tm,VL)
		call matmul2(VL,Vp,VTV)
		!call mat_print(VTV,"VTV")
		call kirchhoff(t,lp,Vp,KS) 
		!call mat_print(KS,"KSp")
		call mat2vec(KS,KSvp)
		!call vec_print(KSvp,6,"KSvp")
	  
	  ! Backward difference
		call mat_print(Bn,"Bn")
		call eigsym33(Bn,Vn,ln)
		!call mat_print(Vn,"V")
		!call vec_print(ln,3,"bn")
		call odgen(ln,1,1.0d0,2.0d0,1.0d0,tn,tB)
	! Kirchhoff Stress
		!call vec_print(tn,3,"tn")
		call kirchhoff(tn,ln,Vn,KSn) 
		!call mat_print(KSn,"KSn")
		call mat2vec(KSn,KSvn)
		!call vec_print(KSvn,6,"KSvn")
	  
		DP=c1*(KSvp-KSvn)
	  
		!call vec_print(DP,6,"DP")
	  
		D(i,:)=DP
	  end do
	  
	  call tensor4ord_print(D,"Jacobian Tangent Matrix")

	  ! Derivative of General Isotropic Tensor Function
	  call mat_print(B,"B")
	  call eigsym33(B,V,l)
	  call odgen(l,1,1.0d0,2.0d0,1.0d0,t,tB)
	  call mat_print(tB,"tB")
	  call vec_print(t,3,"t")
	  call kirchhoff(t,l,V,KS) 
	  
	  call mat_print(KS,"KS")
	  call DGISOT(B,l,t,tB,V,Da)
	  call tensor4ord_print(Da,"D")
	  
    ! Consistent Tanget Matrix
	 !call CST(D,JF,B,stress,dsdePl)
	 !call tensor4ord_print(dsdePl,"Jacobian Tangent Matrix")

end program