      subroutine CSBOP(Fdt,MP,Epsi,CS)
! Compute Stress Based on Odgen Potential
! This subroutine computes Cauchy Stress CS given on Fn+1=Fdt
! Global Varibles 
      double precision Fn1(3,3), CS(3,3),Epsi(3,3)
! LOCAL VARIBLES	  
	  DOUBLE PRECISION B(3,3), V(3,3), l(3), t(3),Fdt(3,3)   
      DOUBLE PRECISION tB(3,3), KS(3,3), J3, D(6,6), s(6)
      Double Precision Dt(6,6), Ds(6,6), Vt(3,3), MP(4),lnc(3,3)
      integer N,i
c ********* Begin of subroutine *******************************88     
      call matmul2(Fdt,Fdt,B)
      call eigsym33(B,V,l)
      call zeros(3,3,lnc)
      do i=1,3
	    lnc(i,i)=0.5d0*log(l(i))
      end do
      call mattransp(V,3,Vt) 
      call eig2mat(V,lnc,Vt,Epsi)
      N=int(MP(1))
      call odgen(l,N,MP(2),MP(3),MP(4),t,tB)	  
      call kirchhoff(t,l,V,KS) 
      call det(Fdt,J3)
      call cauchyS(KS,J3,s)
      call vec2mat(s,CS)
	  
      end subroutine
!end code