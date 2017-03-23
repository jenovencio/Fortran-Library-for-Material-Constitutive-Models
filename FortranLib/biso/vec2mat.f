	  subroutine vec2mat(V,M)
         DOUBLE PRECISION M(3,3), V(6)
		 M(1,1)=V(1)
		 M(2,2)=V(2)
		 M(3,3)=V(3)
		 M(1,2)=V(4)
		 M(2,1)=V(4)
		 M(2,3)=V(5)
		 M(3,2)=V(5)
		 M(1,3)=V(6)
		 M(3,1)=V(6)		
	  end subroutine
!end code   
 
 
