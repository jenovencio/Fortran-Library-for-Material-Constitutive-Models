	  subroutine CDTModgen(B,N,m,a,K,CDTM)
! Central Difference Tangent Matria for Odgen
! this routine generates Odgen Tangent Matrix based on Central Differences
	  integer i,q,N
	  DOUBLE PRECISION B(3,3),l(3),m(N),a(N),t(3),
     &	  tB(3,3),K,CDTM(6,6) 
	  
! Local Variables	  
	  DOUBLE PRECISION DI(3,3),CDM(6,6),
     & bi(6),bcp(6),bcn(6),p,c1
      DOUBLE PRECISION  Bp(3,3),lp(3),Vp(3,3),KSvp(6)
	  DOUBLE PRECISION  Bn(3,3),ln(3),Vn(3,3),KSvn(6)
	  DOUBLE PRECISION  tn(3), tp(3), KSn(3,3), KSp(3,3)	 
	  DOUBLE PRECISION  DP(6),Da(6,6), VL(3,3),VLVt(3,3),LM(3,3),LV(6)
	  DOUBLE PRECISION  tv(6),tm(3,3),vtV(3,3)
	  
! Inicializing Central Difference calculation	  
	  ! p = central difference interanal perturbation 
	  p=1.0d-8  
	  c1=1.0d0/(2*p)
	  
	  call eye(3,3,DI)
	  call eye4sym(CDM)
	  call zeros(6,6,CDTM)	  
	  call mat2vec(B,bi)
	  
	  do i=1,6
		bcp = bi + p*CDM(:,i)
		bcn = bi - p*CDM(:,i)
	    call vec2mat(bcp,Bp)
		call vec2mat(bcn,Bn)
!-----------------------------------------------------------------------    
	! Forward difference
		call eigsym33(Bp,Vp,lp)	
		call odgen(lp,N,m,a,K,tp,tB)
	! Kirchhoff Stress
		call kirchhoff(tp,lp,Vp,KSp) 
		call mat2vec(KSp,KSvp)
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------	  
	  ! Backward difference
		call eigsym33(Bn,Vn,ln)
		call odgen(ln,N,m,a,K,tn,tB)
	! Kirchhoff Stress
		call kirchhoff(tn,ln,Vn,KSn) 
		call mat2vec(KSn,KSvn)

!-----------------------------------------------------------------------	  
		
		DP=c1*(KSvp-KSvn)
		CDTM(i,:)=DP
	  end do

	  
	  end subroutine
!end code