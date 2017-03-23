      subroutine Fn12Etrial(Fn1,Fn,etrial)
      double precision Fn1(3,3),Fn(3,3),etrial(6)
	  
!local Variables	  
      double precision Fninv(3,3),Fdelta(3,3),Bn(3,3),en1(3,3),
     & Emat(3,3),Btrial(3,3),bi(3),M(3,3),lnbi(3,3)
      integer i,j
c
	  Fninv=Fn 
      call matinv(Fninv,3)
      Fdelta=Fn1*Fninv
c	  
c   Compute elastic trial
      Bn=matmul(Fn,transpose(Fn))	  
	  Btrial=matmul(Fdelta,(matmul(Bn,transpose(Fdelta))))
      call eigsym33(Btrial,M,bi) 
      call zeros(3,3,lnbi)
      do i=1,3
		    lnbi(i,i)=log(bi(i)) 
      end do	 
	  
      Emat=0.5d0*matmul(M,(matmul(lnbi,transpose(M))))
	  
	  etrial(1)=Emat(1,1)
	  etrial(2)=Emat(2,2)
	  etrial(3)=Emat(3,3)
	  etrial(4)=2.0d0*Emat(1,2)
	  etrial(5)=2.0d0*Emat(2,3)
	  etrial(6)=2.0d0*Emat(1,3)
	  
      end subroutine
!end code	  