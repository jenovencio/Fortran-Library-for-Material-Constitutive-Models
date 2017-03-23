      subroutine DGISOT(Xm,x,y,dy,V,D)
! This subroutine computes the derivative of a general isotropic tensor function
! D(Y(X))/DX
! This subroutine needs of others subroutines:
! eigvec2eigproj.f90
! eye4sym.f90
! DX2DX.F90
! KroIsoProd.f90
! SCAMULT.f90
! zeros.f90
! Inputs
! Xm - matrix 3 by 3
! x - vector of eigen values of X
! y - vector of eigen values of Y
! dy - derivatives dyi/dxj - matrix 3 by 3
! V - matrix (3 by 3) with X eigen vectors
!Output
! D - derivative of the tensor function Y(X) - 6 by 6 matrix 
! Local Variable
! Is - Fourth ordem simmetric identity
! DX2 - derivative of X^2  - D(X^2)/DX
	  DOUBLE PRECISION y(3),x(3),dy(3,3),V(3,3),D(6,6)
	  DOUBLE PRECISION Xm(3,3)
	  DOUBLE PRECISION E(3,3,3)
	  DOUBLE PRECISION Is(6,6),K(3,3)
	  DOUBLE PRECISION DX2(6,6)
	  integer i,j
	  integer a(3),b(3),c(3)
	  DOUBLE PRECISION D1(6,6),D2(6,6),D3(6,6),D4(6,6)
	  DOUBLE PRECISION XX(6,6),XI(6,6),IX(6,6),II(6,6)
	  DOUBLE PRECISION Ea(3,3),Eb(3,3),Ec(3,3),Ei(3,3),Ej(3,3)
	  DOUBLE PRECISION ya,xa,xb,xc,c1,c2,c3,c4
	  DOUBLE PRECISION s1, s2, s3, s4, s5, s6
	  DOUBLE PRECISION small,dif12, dif13, dif23
	  
	  call DX2DX(Xm,DX2)
	  call eye4sym(Is)
	  call zeros(6,6,D)
!	  call eigvec2eigproj(V,E)
      call eigproj(Xm,x,E)
      call eye(3,3,K)
	  small=1.0d-5
	  dif12=abs(x(1)-x(2))
	  dif13=abs(x(1)-x(3))
	  dif23=abs(x(2)-x(3))
	  
! Check eigen values	    
		if ( (dif12 .gt. small) .and. (dif13 .gt. small) .and. 
     & 	 (dif23 .gt. small)) then 
		goto 500
		else if ( (dif12 .lt. small) .and. (dif13 .lt. small) .and. 
     & 	 (dif23 .lt. small)) then 
		goto 700
		else
		goto 600
		end if	
	
500   continue	
	  
! Permutation vectors
	  a(1)=1
	  a(2)=2
	  a(3)=3
	  b(1)=3
	  b(2)=1
	  b(3)=2	
	  c(1)=2
	  c(2)=3
	  c(3)=1
	  
	  do i=1,3
		ya=y(a(i))
		xa=x(a(i))
		xb=x(b(i))
		xc=x(c(i))
		c1=ya/((xa-xb)*(xa-xc))
		c2=xb+xc
		c3=(xa-xb)+(xa-xc)
		c4=(xb-xc)
		Ea=E(a(i),:,:)
		Eb=E(b(i),:,:)
		Ec=E(c(i),:,:)
		call KroIsoProd(Ea,Ea,D1) 
		call KroIsoProd(Eb,Eb,D2)
		call KroIsoProd(Ec,Ec,D3)
		D=D+c1*(DX2 - c2*Is - c3*D1 - c4*(D2 - D3))
	  end do	
!		call tensor4ord_print(DX2,"DX2")
!		 call tensor4ord_print(Is,"Is")
!		 call tensor4ord_print(D1,"D1")
!		 call tensor4ord_print(D2,"D2")
!		 call tensor4ord_print(D3,"D3")

	  do i=1,3
		do j=1,3
			Ei=E(i,:,:)
			Ej=E(j,:,:)
			call KroIsoProd(Ei,Ej,D4) 
			D=D+dy(i,j)*D4
		end do
	  end do		
	  goto 1000
	
600   continue
 
		    if (dif23 .lt. small) then ! xa~=xb=xc
                  !print *, "x1~=x2=x3"
				  ! Permutation vectors
				    a(1)=1
					b(1)=2
					c(1)=3
			else if (dif12 .lt. small) then ! xa=xb~=xc		
				 !print *, "x3~=x2=x1"
				! Permutation vectors	
					a(1)=3
					c(1)=1
					b(1)=2			
			else  ! xa=xc~=xb
			     !print *, "x2~=x3=x1"
				! Permutation vectors	
					a(1)=2
					c(1)=3
					b(1)=1				
			end if
			
     	    i=1			
			ya=y(a(i))
			yb=y(b(i))
			yc=y(c(i))
			xa=x(a(i))
			xb=x(b(i))
			xc=x(c(i))
			s1=(ya -yc)/(xa-xc)**2 +
     &      1.0d0/(xa-xc)*(dy(c(i),b(i))-dy(c(i),c(i)))
            s2=2.0d0*xc*(ya -yc)/(xa-xc)**2 +
     &      (xa+xc)/(xa-xc)*(dy(c(i),b(i))-dy(c(i),c(i)))
			s3=2.0d0*(ya -yc)/(xa-xc)**3 +
     &		1.0d0/(xa-xc)**2*(dy(a(i),c(i)) + dy(c(i),a(i)) -
     &       dy(a(i),a(i)) - dy(c(i),c(i)))	
			s4=2.0d0*xc*(ya -yc)/(xa-xc)**3 +
     &      1.0d0/(xa-xc)*(dy(a(i),c(i)) - dy(c(i),b(i))) +	
     &      xc/(xa-xc)**2*(dy(a(i),c(i)) + dy(c(i),a(i)) -
     &	    dy(a(i),a(i)) - dy(c(i),c(i)))
			s5=2.0d0*xc*(ya -yc)/(xa-xc)**3 +
     &      1.0d0/(xa-xc)*(dy(c(i),a(i)) - dy(c(i),b(i))) +	
     &      xc/(xa-xc)**2*(dy(a(i),c(i)) + dy(c(i),a(i)) -
     &	    dy(a(i),a(i)) - dy(c(i),c(i)))
			s6=2.0d0*xc*xc*(ya -yc)/(xa-xc)**3 +
     &      xa*xc/(xa-xc)**2*(dy(a(i),c(i)) + dy(c(i),a(i))) -
     &      xc*xc/(xa-xc)**2*(dy(a(i),a(i)) + dy(c(i),c(i))) -
     &      (xa+xc)/(xa-xc)*dy(c(i),b(i))	

	     call KroIsoProd(Xm,Xm,XX)
		 call KroIsoProd(Xm,K,XI)
		 call KroIsoProd(K,Xm,IX)
		 call KroIsoProd(K,K,II)
		 !call mat_print(K,"I")
		 !call mat_print(Xm,"B")
!         print *, "s1=",s1
!		  print *, "s2=",s2
!		  print *, "s3=",s3
!		  print *, "s4=",s4
!		  print *, "s5=",s5
!		  print *, "s6=",s6
!		 call tensor4ord_print(DX2,"DX2")
!		 call tensor4ord_print(Is,"Is")
!		 call tensor4ord_print(XX,"XX")
!		 call tensor4ord_print(XI,"XI")
!		 call tensor4ord_print(IX,"IX")
!		 call tensor4ord_print(II,"II")
         !print *, "s1=",s1
		 D=D+s1*DX2 - s2*Is - s3*XX + s4*XI 
     &	 + s5*IX - s6*II  	 

      goto 1000
700   continue
      !print *, "x3=x2=x1"
	  call KroIsoProd(K,K,II)
	  s1=dy(1,1)-dy(1,2)
	  s2=dy(1,2)
	  D =s1*Is + s2*II 
	  
	  
	  goto 1000
1000   continue

	
      end subroutine
!end code
