      subroutine plasticLaw(e,G,K,H,dpleq,qtrial,N,Dep,strial)
!global variables      
      double precision e(6),G,dpleq,qtrial,strial(6),Dep(6,6)
	  
! local variables	  
      double precision Idev(6,6),Is(6,6),II(6,6),De(6,6),D(6,6),
     & lampda,con1,con2,N(6),NN(6,6),twoG,threeG,TWOTHIRD,JM(6,6)
c
      integer i,j	  
	  
      twoG=2.0d0*G 
      threeG=3.0d0*G
	  TWOTHIRD=(2.0d0/3.0d0)
      call zeros(6,6,Is)
      Is(1,1)=1.0d0
      Is(2,2)=1.0d0  
      Is(3,3)=1.0d0
      Is(4,4)=0.50d0
      Is(5,5)=0.50d0
      Is(6,6)=0.50d0  
	  
      call zeros(6,6,II)
      do i=1,3
	     do j=1,3
            II(i,j)=1.0d0
         end do
      end do
      
      Idev=Is-(1.0d0/3.0d0)*II
	  
	  lampda=K-(2.0d0/3.0d0)*G
      call zeros(6,6,De)

! elastic Jacobian Matrix	
 	   
      do i=1,3
	     do j=1,3
	       De(i,j)=lampda
		   NN(i,j)=N(i)*N(j)
         end do
		    De(i+3,i+3)=G
            De(i,i)=De(i,i)+twoG
      end do		 
	  
      con1 = threeG * dpleq / qtrial
      con2 = threeG/(threeG+H) - con1
      con2 = TWOTHIRD * con2

	  D=De-con1*Idev

	  strial=matmul(D,e)
!      Dep = De - twoG(con1*Idev + con2*NN)
      DO i=1,ncomp
         DO j=1,ncomp
            JM(j,i) = ZERO
         END DO
      END DO
      DO i=1,nDirect
         DO j=1,nDirect
            JM(i,j) = -THIRD
         END DO
         JM(i,i) = JM(i,i) + ONE
      END DO
      DO i=nDirect + 1,ncomp
         JM(i,i) = HALF
      END DO
      DO i=1,ncomp
         DO j=1,ncomp
            Dep(i,j) =    De(i,j) - twoG
     &           * (  con2 * N(i) * N(j) + con1 * JM(i,j) )
         END DO
      END DO	 
 	  
      return
      end subroutine
!end code	  
!end code	  