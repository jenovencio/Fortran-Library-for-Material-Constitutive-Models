      subroutine ComputeStress(Fn1,Fn,Fvn,Fpn,dt,bmax,Qn,MP,
     & Cpr,cprv,Mpr,q,bet,vM,Piola,Cauchy,Epsi,Fvn1,Fpn1,Qn1,deltaQ,
     & norma)
c
c
c***************** User defined part *************************************	  
! User material 
! Variable declaration  
c***************************************************************************
! Global Variable Declaration
! Input Variables
      DOUBLE PRECISION Fn1(3,3), Fn(3,3), Fvn(3,3), Fpn(3,3)
      DOUBLE PRECISION dt, Qn, MP(16)
      integer bmax 
! Output 
      DOUBLE PRECISION vM,Epsi(3,3),Piola(3,3), Cauchy(3,3),
     &  Qn1, Fpn1(3,3), Fvn1(3,3), deltaQ, norma, bet	  
c	  
c***************************************************************************
	 !lobal Variable Declaration	  
      DOUBLE PRECISION Fhn1(3,3),Fvol(3,3), Chn1(3,3), Cn1(3,3),Cn(3,3)
      DOUBLE PRECISION Mn1(3,3),c1(3),e1(3),En1(3,3),M(3,3,3) 
      DOUBLE PRECISION dFiChmv(3,3), dFiChmp(3,3),q(3)
      DOUBLE PRECISION Fpr(3,3),FpnInv(3,3)
      DOUBLE PRECISION Cpr(3,3),Mpr(3,3),cprv(3), epr(3), Betab,Taub
      DOUBLE PRECISION q0(3), deltaQc, ee_p(3), deltaQqj(3)
      DOUBLE PRECISION deltaQq(3,3), MprProj(3,3,3), dFiedCpr(3,3)
      DOUBLE PRECISION dFiedChp(3,3), FpnInvt(3,3), dU, Svol(3,3)
      DOUBLE PRECISION Cn1Inv(3,3), dC(3,3), Sdevp(3,3),Sp(3,3)
      DOUBLE PRECISION Piolap(3,3), Caudevp(3,3),PFn1t(3,3),Ident(3,3)
      DOUBLE PRECISION Cau2, eigVec(3,3), eigVal(3), Jn1
      integer i,j, qtbmax, verifdQ, a
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
!**************************************************************************************
!                         Printing inputs
!**************************************************************************************
!      call mat_print(Fn1,"Fn1")
!      call mat_print(Fn,"Fn")
!      call mat_print(Fpn,"Fpn")
!      print *, "dt=", dt
!      print *, "Qn=", Qn
!**************************************************************************************
!    Compute Stress	  
      call det(Fn1,Jn1)
c
! Decompositon of F in isochoric (Fhn1) and volumetric (Fvol) 	  
      call isovol(Fn1,Fhn1,Fvol)
c	  
! Generating rigth Cauchy-Green tensor Chn1, Cn1, Cn
      call matmul1(Fhn1,Fhn1,Chn1)
      call matmul1(Fn1,Fn1,Cn1)
      call matmul1(Fn,Fn,Cn)
c	  
! Evaluating Eigen values and Eigen Vectors of Chn1
      call eigsym33(Chn1,Mn1,c1)
	  e1=0.5d0*log(c1)
      call zeros(3,3,En1)
      call eigproj(Chn1,c1,M)
c
! Calculating Logarithm Strain Tensor	  
      do i=1,3
	  		En1=En1+e1(i)*M(i,:,:)
      end do       
c	  
! Deviatoric stress. contribution of maxwell spring.	
! Evaluating on it each maxwell component
      call zeros(3,3,dFiChmv)
      call zeros(3,3,dFiChmp)
      qtbmax=1
! caso nao entre nos dois laços
	  deltaQ=0.0d0
	  q(1)=0.0d0
	  q(2)=0.0d0
	  q(3)=0.0d0
	  Qn1=0.0d0
      call eye(3,3,Fpn1)
      call eye(3,3,Fvn1)
c
c******************************************************************
c            CALCULATING PLASTIC CONTRIBUTION
c******************************************************************
	  FpnInv=Fpn
      call matinv(FpnInv,3)
      call zeros(3,3,Fpr)
      call matmul3(Fhn1,FpnInv,Fpr)
      call zeros(3,3,Cpr)
      call matmul1(Fpr,Fpr,Cpr)
c
c	  
      call eigsym33(Cpr,Mpr,cprv)
      call eigproj(Cpr,cprv,MprProj)	  
      epr=0.5d0*log(cprv)
c
c	  
	  Betab=MP(5+2*qtbmax)
      Taub=MP(6+2*qtbmax)
c	  
c
      call verifdeltaQzero(epr,Betab,Taub,dt,Qn,q,verifdQ)
c
c	  
      if (verifdQ .eq. 1) then
        q0=q
        call ComputeMinimization(epr,Betab,Taub,dt,Qn,q0,q,deltaQ,norma
     &        ,bet)
        deltaQc=deltaQ
      else
          deltaQ=0.0d0
          q(1)=0.0d0
          q(2)=0.0d0
          q(3)=0.0d0
		  bet=0.0d0
          norma=0.0d0  
      end if
c     
c     
      if (deltaQ .lt. 0.0d0) then
		  deltaQ=0.0d0
      end if
c
      ee_p=epr-deltaQ*q
      Qn1=Qn+deltaQ
 !    deltaQqj = exp(dp);   dp = deltaQ*q - Where dp are the eigen values of dt*Dp - Plastic streching :  L = D + W ; W=0 
	  deltaQqj=exp(deltaQ*q)
c	  
      call zeros(3,3,deltaQq)
c     
c
!       deltaQq = dt*Dp
      do a=1,3
        deltaQq=deltaQq+deltaQqj(a)*MprProj(a,:,:)
      end do 
c
c	  
!	  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c    
!    %Gradient of viscous deformation n+1
!     Fpn+1=(dt*Dp)*Fpn 
      Fpn1=matmul(deltaQq,Fpn)
c
c
c	  
      call ComputeMaxwellStress(ee_p,Betab,cprv,MprProj,dFiedCpr)
!      Transposing (Fpn^-1)'
      call mattransp(FpnInv,3,FpnInvt)
c
!      calculo do dFie/dCh 	  = (Fpn^-1)*(dFiedCpr)*(Fpn^-1)'
      call eig2mat(FpnInv,dFiedCpr,FpnInvt,dFiedChp)
	  dFiedChp=2.0d0*dFiedChp
c
      dFiChmv=dFiChmv+dFiedChp
      qtbmax=qtbmax + 1
c	  
!***************************************************************************
!******************End of CALCULATING PLASTIC CONTRIBUTION******************	  
!***************************************************************************	  
      call dVolu(Jn1,dU)
	  Cn1Inv=Cn1
      call matinv(Cn1Inv,3)
!     Pvol = (dU/dJ)*J*(C^-1)
	  Svol=dU*Jn1*Cn1Inv
c
c
! Tensao nos bracos visco 	 
! Falta implementar
c	  
! Tensao nos bracos plast
      call DEVC(dFiedChp,Cn1,dC)
c
	  Sdevp=(Jn1**(-2.0d0/3.0d0))*dC
	  Sp=Sdevp + Svol
!	  Sp=Sdevp
	  Piolap=matmul(Fn1,Sp)
      call matmul2(Piolap,Fn1,PFn1t)
	  Caudevp=(1.0d0/Jn1)*PFn1t
      call eye(3,3,Ident)
!	  Caudevp=Caudevp-Caudevp(3,3)*Ident
c	  
!S=S_p+S_v;
!Piola=Piolap + Piolav
!Cauchy=Cauchyp + Cauchyv	  
	  Piola=Piolap	  
	  Cauchy=Caudevp	  
      call prodint(Cauchy,Cauchy,Cau2)
	  vM=(1.5d0*Cau2)**0.5d0
	  Epsi=En1
!	  call eigsym33(Cauchy,eigVec,eigVal)
c     
!***************************************************************************
!*************************End of Subroutine*********************************	  
!***************************************************************************	
!***************************************************************************
!                         Printing outputs
!***************************************************************************
!      call mat_print(Cauchy,"Cauchy")
!      call mat_print(Epsi,"Epsi")
! 	  print *, "norma=",norma
      return  	   
      end subroutine
!end code	  
