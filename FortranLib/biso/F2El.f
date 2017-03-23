      subroutine F2El(F,El)
! This routine get F and transform in Logarithm Strain (El)
! Global variable declaration
      double precision F(3,3), El(3,3)
!Local Variable declaration      
      double precision B(3,3), Bi(3), V(3,3),Vt(3,3)
      double precision e(3,3)
      integer i
	  
      Call matmul2(F,F,B)
      Call eigsym33(B,V,Bi)
      Call zeros(3,3,e)
      do i=1,3
	     e(i,i)=0.5D0*log(Bi(i))
      end do
      Call mattransp(V,3,Vt)
      Call eig2mat(V,e,Vt,El)
	  
	  
      end subroutine
!end Bode	  