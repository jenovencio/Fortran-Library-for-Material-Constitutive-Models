      subroutine CTMOP(Fn1,Cauchy,Epsi,MP,D)
c***************************************************************************
! Compute Tangent Matrix for Odgen Potential - CTMOP
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
      DOUBLE PRECISION Fn1(3,3),Cauchy(3,3)
	  DOUBLE PRECISION MP(16), Epsi(3,3)
	  integer bmax 
! Output 
      DOUBLE PRECISION D(6,6)
  
	  
c***************************************************************************
	 !lobal Variable Declaration	
      DOUBLE PRECISION inc,CSv(6),Ev(6), incn
      DOUBLE PRECISION dCSp(6,6), dEpsip(6,6), small
      DOUBLE PRECISION dCSn(6,6), dEpsin(6,6),dEinv(6,6)	  
      DOUBLE PRECISION dCS(6,6), dEpsi(6,6), Fn1s(3,3)
      integer t,i,j
	  

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
      inc=1.0d-5
	  small=1.0d-6

      call zeros(6,6,D)
 ! Store initial Fn1     
	  Fn1s=Fn1
      call varOdgen(Fn1,MP,inc,dCSp,dEpsip)
!      call tensor4ord_print(dCSp,"dCSp")
	  incn=-inc
      Fn1=Fn1s
	  call varOdgen(Fn1,MP,incn,dCSn,dEpsin)
!      call tensor4ord_print(dCSn,"dCSn")
	  call mat2vec(Cauchy,CSv)
!	  call vec_print(CSv,6,"CSv")
      call mat2vec(Epsi,Ev)
      
      call zeros(6,6,dCS)
      call zeros(6,6,dEpsi)
      ! for t=1 -> Forward differeces, t=2 -> Central diferences
	  t=1
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
	  
	  
!      dEinv=dEpsi
!      call matinv(dEinv,6)
	    
	   
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
