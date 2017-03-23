      subroutine ComNumTanMat(Fn1,Fn,Fvn,Fpn,dt,bmax,Qn,MP,dFiedCpr0,C)
! Global Variable Declaration
      DOUBLE PRECISION Fn1(3,3), Fn(3,3), Fvn(3,3), Fpn(3,3)
      DOUBLE PRECISION dt, Qn, MP(16),C(6,6),dFiedCpr0(3,3)
	  
	  ! Output 
      DOUBLE PRECISION vM,Epsi(3,3),Piola(3,3), Cauchy(3,3),
     &  Qn1, Fpn1(3,3), Fvn1(3,3), deltaQ, norma
	  
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
      DOUBLE PRECISION Cau2, eigVec(3,3), eigVal(3), Jn1,inc,
     & P(6,6),chv(6),chvinc(6),ChS(3,3),dF0(6),cv(6)
      integer i,j, qtbmax, verifdQ, a
      DOUBLE PRECISION DP(5)
      common /MatDP/ DP	  
!****************************************************************
!      PARAMETER FOR FORWARD DIFFERENCES DERIVATIVE
!**************************************************************** 
	  inc=1.0d-3
      call eye(6,6,P)
      call zeros(6,6,C)
	  call mat2vec(dFiedCpr0,dF0)
!**************************************************************** 
! compute det de Fn1 
      call det(Fn1,Jn1)

! Decompositon of F in isochoric (Fhn1) and volumetric (Fvol) 	  
      call isovol(Fn1,Fhn1,Fvol)
	  
! Generating rigth Cauchy-Green tensor Chn1, Cn1, Cn
      call matmul1(Fhn1,Fhn1,Chn1)
      call matmul1(Fn1,Fn1,Cn1)
      call matmul1(Fn,Fn,Cn)

!****************************************************************	  	  
!	               Compute Stress
!****************************************************************	  
!	   Evaluating Eigen values and Eigen Vectors of Chn1
      call eigsym33(Chn1,Mn1,c1)
	  e1=0.5d0*log(c1)
      call zeros(3,3,En1)
      call eigproj(Chn1,c1,M)

! Calculating Logarithm Strain Tensor	  
      do i=1,3
	  		En1=En1+e1(i)*M(i,:,:)
      end do       
	  
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

c******************************************************************
c            CALCULATING PLASTIC CONTRIBUTION
c******************************************************************
	  FpnInv=Fpn
      call matinv(FpnInv,3)
      call zeros(3,3,Fpr)
      call matmul3(Fhn1,FpnInv,Fpr)
      call zeros(3,3,Cpr)
      call matmul1(Fpr,Fpr,Cpr)
      ChS=Cpr

! Loop start	  
      do j=1,6	  

!	  print *,"j=",j    
	  call mat2vec(ChS,chv)
	  chvinc=chv+inc*P(1:6,j)
	  call vec2mat(chvinc,Cpr)
!     call mat_print(Cpr,"Cpr")
	  
      call eigsym33(Cpr,Mpr,cprv)
      call eigproj(Cpr,cprv,MprProj)	  
      epr=0.5d0*log(cprv)

	  
	  Betab=MP(5+2*qtbmax)
      Taub=MP(6+2*qtbmax)
	  

      call verifdeltaQzero(epr,Betab,Taub,dt,Qn,q,verifdQ)

	  
      if (verifdQ .eq. 1) then
        q0=q
         call ComputeMinimization(epr,Betab,Taub,dt,Qn,q0,q,deltaQ,norma
     &        ,bet)
        deltaQc=deltaQ
      else
          deltaQ=0.0d0
          q(1)=0.0d0
          q(2)=0.0d0
		  norma=0.0d0
          q(3)=0.0d0
            
      end if
      
      
      if (deltaQ .lt. 0.0d0) then
		  deltaQ=0.0d0
      end if

      ee_p=epr-deltaQ*q
      Qn1=Qn+deltaQ
 !    deltaQqj = exp(dp);   dp = deltaQ*q - Where dp are the eigen values of dt*Dp - Plastic streching :  L = D + W ; W=0 
	  deltaQqj=exp(deltaQ*q)
	  
      call zeros(3,3,deltaQq)
      

!       deltaQq = dt*Dp
      do a=1,3
        deltaQq=deltaQq+deltaQqj(a)*MprProj(a,:,:)
      end do 

	  
!	  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
!    %Gradient of viscous deformation n+1
!     Fpn+1=(dt*Dp)*Fpn 
      Fpn1=matmul(deltaQq,Fpn)

	  
      call ComputeMaxwellStress(ee_p,Betab,cprv,MprProj,dFiedCpr)
      call mat2vec(dFiedCpr,cv)
!      Transposing (Fpn^-1)'
!      call mattransp(FpnInv,3,FpnInvt)

!      calculo do dFie/dCh 	  = (Fpn^-1)*(dFiedCpr)*(Fpn^-1)'
!      call eig2mat(FpnInv,dFiedCpr,FpnInvt,dFiedChp)


!****************************************************************	 
         do i=1,6
	         C(j,i)=(1.0d0/inc)*(cv(i)-dF0(i))
         end do
	 
 !     call tensor4ord_print(C,"C")
      end do
      return
      end subroutine
!end code	  