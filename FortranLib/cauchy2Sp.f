      subroutine cauchy2Sp(Cauchy,Fn1,Sp)
      double precision Cauchy(3,3),Fn1(3,3),
     & Sp(3,3)
      double precision Jn1,Fn1inv(3,3),Ks(3,3)	 
c	  
c	  
c     Compute Kirshhoff Stress	  
      Fn1inv=Fn1
      call matinv(Fn1inv,3)
      call det(Fn1,Jn1)
	  Ks=Jn1*Cauchy
c	  
      Sp=matmul(Fn1inv,(matmul(Ks,transpose(Fn1inv))))
      return     
      end subroutine
!end code	  