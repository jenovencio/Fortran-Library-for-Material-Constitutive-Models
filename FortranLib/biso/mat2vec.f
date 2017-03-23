      subroutine mat2vec(M,V)
        DOUBLE PRECISION M(3,3), V(6)
		DOUBLE PRECISION dif12,dif23,dif13,small
		small=1.0d-2
		dif12=M(1,2) - M(2,1)
		dif23=M(2,3) - M(3,2)
		dif13=M(1,3) - M(3,1)
		if (dif12 .gt. small) then
		    print *,"Matrix is not symmetryc"
		end if	
		if (dif23 .gt. small) then
		    print *,"Matrix is not symmetryc"
		end if	
        if (dif13 .gt. small) then
		    print *,"Matrix is not symmetryc"
         end if			
		
		 V(1)=M(1,1)
		 V(2)=M(2,2)
		 V(3)=M(3,3)
		 V(4)=(M(1,2)+M(2,1))*0.5d0
		 V(5)=(M(2,3)+M(3,2))*0.5d0
		 V(6)=(M(1,3)+M(3,1))*0.5d0
		
		
      end subroutine
!end code  
 
