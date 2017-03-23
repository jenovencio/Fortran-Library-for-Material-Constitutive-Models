      subroutine ComputeTanMat(Fn1,Fn,Fvn,Fpn,dt,bmax,Qn,MP,
     & Cauchy,Epsi,D)
c***************************************************************************
! This subroutine evaluate Material Jacobian Matrix (Tangent Matrix) based on
! central differences derivatives
! This subroutine evaluate variation in Cauchy Stress (dCS) and variation on logarithm 
! Strain based on Forward differences perturbation
!  Inputs
! Fn1,Fn,Fvn,Fpn,dt,bmax,Qn,MP
!
! Output
! dCSdEpsi -> derivative of Cachy Stress (CS) in relation to logarithm Strain (Epsi)
! Variable declaration  
c***************************************************************************
! Global Variable Declaration
! Input Variables
      DOUBLE PRECISION Fn1(3,3), Fn(3,3), Fvn(3,3), Fpn(3,3)
	  DOUBLE PRECISION dt, Qn, MP(16)
	  integer bmax 
! Output 
      DOUBLE PRECISION D(6,6)
  
	  
c***************************************************************************
	 !lobal Variable Declaration	
      DOUBLE PRECISION vM,Epsi(3,3),Piola(3,3), Cauchy(3,3),
     &  Qn1, Fpn1(3,3), Fvn1(3,3), deltaQ, nor,CSv(6),Ev(6)		 
      DOUBLE PRECISION Fhn1(3,3),Fvol(3,3), Chn1(3,3), Cn1(3,3),Cn(3,3)
	  DOUBLE PRECISION Mn1(3,3),c1(3),e1(3),En1(3,3),M(3,3,3) 
	  DOUBLE PRECISION dFiChmv(3,3), dFiChmp(3,3),q(3)
	  DOUBLE PRECISION Fpr(3,3),FpnInv(3,3),small
	  DOUBLE PRECISION Cpr(3,3),Mpr(3,3),cprv(3), epr(3), Betab,Taub
	  DOUBLE PRECISION q0(3), deltaQc, ee_p(3), deltaQqj(3)
	  DOUBLE PRECISION deltaQq(3,3), MprProj(3,3,3), dFiedCpr(3,3)
	  DOUBLE PRECISION dFiedChp(3,3), FpnInvt(3,3), dU, Svol(3,3)
	  DOUBLE PRECISION Cn1Inv(3,3), dC(3,3), Sdevp(3,3),Sp(3,3)
	  DOUBLE PRECISION Piolap(3,3), Caudevp(3,3),PFn1t(3,3),Ident(3,3)
	  DOUBLE PRECISION Cau2, eigVec(3,3), eigVal(3), Jn1, incn
	  DOUBLE PRECISION dCSp(6,6),dEpsip(6,6), dCSn(6,6), dEpsin(6,6)
	  DOUBLE PRECISION dCS(6,6),dEpsi(6,6),dEpinv(6,6), inc, Fn1s(3,3)
	  integer i,j, qtbmax, verifdQ, a, t
      DOUBLE PRECISION DP(5)
	  common /MatDP/ DP
! This subroutine needs:
! verifdeltaQzero
! MPOgden4.f
! isovol.f
! eye.f
! vec2mat.f
! mat_print.f
! det.f
! matmul1.f
! etc ...
! See f.bat for all subroutines needed
      inc=1.0d-6
      small=1.0d-10
      call zeros(6,6,D)
! Store initial Fn1      
	  Fn1s=Fn1 
      call vardCdEpsiFD(Fn1,Fn,Fvn,Fpn,dt,bmax,Qn,MP,inc,
     & dCSp,dEpsip)
      
      incn=-inc
	  Fn1=Fn1s
      call vardCdEpsiFD(Fn1,Fn,Fvn,Fpn,dt,bmax,Qn,MP,incn,
     & dCSn,dEpsin)	 
      

 

	  call mat2vec(Cauchy,CSv)
!	  call vec_print(CSv,6,"CSv")
	  call mat2vec(Epsi,Ev)
      call zeros(6,6,dCS)
	  call zeros(6,6,dEpsi)
! for t=1 -> Forward differeces, t=2 -> Central diferences
      t=2
	  
      if (t .eq. 1) then
         do i=1,6
		    dCS(i,:)=dCSp(i,:)-CSv
		    dEpsi(i,:)=dEpsip(i,:)-Ev
		    D(i,:)=dCS(i,:)/inc
         end do
      else
         do i=1,6
		    dCS(i,:)=dCSp(i,:)-dCSn(i,:)
		    dEpsi(i,:)=dEpsip(i,:)-dEpsin(i,:)
		    D(i,:)=dCS(i,:)/(2.0d0*inc)
         end do  
      end if	
! Check small values	  
      do i=1,6
         do j=1,6
            if (abs(D(i,j)) .le. small) then
			    D(i,j)=0.0d0
            end if
         end do
      end do		 
	  
	    
	   
!  checking local variables	 
   
!      print *, "D - first"
!      call tensor4ord_print(D,6)
!      print *,"dCS"
!      call tensor4ord_print(dCS,6)
!      print *,"dEpsi"
!      call tensor4ord_print(dEpsi,6)
	
      return  	   
      end subroutine
!end code	  
