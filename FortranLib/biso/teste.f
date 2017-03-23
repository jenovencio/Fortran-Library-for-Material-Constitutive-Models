      program teste
	  double precision s(6),qtrial,N(6),JM(6,6),
     & ZERO,THIRD,ONE,HALF
	  double precision Idev(6,6),Is(6,6),II(6,6),De(6,6),lampda,con
	  INTEGER i,j,ncomp,nDirect
	  s(1)=0.d0
	  s(2)=0.d0
	  s(3)=0.d0
	  s(4)=10.d0
	  s(5)=10.d0
	  s(6)=10.d0
	
	  ZERO=0.0d0
	  THIRD=1.0d0/3.0d0
	  ONE=1.0d0
	  HALF=0.5d0
	  
	  ncomp=6
	  nDirect=3
	  
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

	

	  call tensor4ord_print(JM,"JM")
	  call tensor4ord_print(Idev,"Idev")
	  end program