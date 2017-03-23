      subroutine eigproj(X,xi,M)
! This subroutine create eigen projections matrix using eigen values xi, and the matrix X
! X - Matrix 3 by 3 of eigen vector 
! M(k) - 3 by 3 symmetric matrix eigen projection k, k=1,3
! M - 3 by 3 by 3 array 
      DOUBLE PRECISION X(3,3),M(3,3,3),xi(3),c,I1,I3,Ident(3,3),
     & dif12,dif13,dif23, small  
      integer i,j,k
      call eye(3,3,Ident)
	  small=1.0d-8
	  dif12=abs(xi(1)-xi(2))
	  dif13=abs(xi(1)-xi(3))
	  dif23=abs(xi(3)-xi(2))  
	  
        if ( (dif12 .gt. small) .and. (dif13 .gt. small) .and. 
     & 	 (dif23 .gt. small)) then 
		    goto 100
         else if ( (dif12 .lt. small) .and. (dif13 .lt. small) .and. 
     & 	 (dif23 .lt. small)) then 
            goto 200
         else
		    goto 300
         end if	
	  
	  

!***************************************************************************
100   continue	  
!	  x(1)~=x(2)  and x(2)~=x(3)   and x(3)~=x(1)
	  I1=xi(1)+xi(2)+xi(3)
	  I3=xi(1)*xi(2)*xi(3)

      do i=1,3
	     c=xi(i)/(2.0d0*xi(i)**3-I1*xi(i)*xi(i)+I3)
         M(i,1:3,1:3)=c*(matmul(X,X)-(I1-xi(i))*X+(I3/xi(i))*Ident)
      end do		
      goto 400
!***************************************************************************

!***************************************************************************	  
200   continue
!      x(1)=x(2)=x(3)  
       M(1,:,:)=Ident
       M(2,1:3,1:3)=0.0d0
       M(3,1:3,1:3)=0.0d0
      goto 400
!***************************************************************************

!***************************************************************************	 
300   continue
!	  x(1)=x(2)  or x(2)=x(3)   or x(3)=x(1)
      if (dif12 .lt. small) then
         ! x(1)=x(2) , x(j)=x(k)              
          i=3
          j=1
          k=2
      else if (dif23 .lt. small) then	  
         ! x(2)=x(3) , x(j)=x(k)              
          i=1
          j=2
          k=3
      else	  
	      ! x(1)=x(3) , x(j)=x(k)              
          i=2
          j=1
          k=3
      end if
		  
	   c=xi(i)/(2.0d0*xi(i)**3-I1*xi(i)*xi(i)+I3)
       M(i,1:3,1:3)=c*(matmul(X,X)-(I1-xi(i))*X+(I3/xi(i))*Ident)
	   M(j,1:3,1:3)=Ident(1:3,1:3) - M(i,1:3,1:3)
	   M(k,1:3,1:3)=0.0d0
	   
400   continue	  
      end subroutine
!end code  
