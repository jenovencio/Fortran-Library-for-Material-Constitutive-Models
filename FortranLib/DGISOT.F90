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
	  DOUBLE PRECISION Is(6,6)
	  DOUBLE PRECISION DX2(6,6)
	  integer i,j
	  integer a(3),b(3),c(3)
	  DOUBLE PRECISION D1(6,6),D2(6,6),D3(6,6),D4(6,6)
	  DOUBLE PRECISION Ea(3,3),Eb(3,3),Ec(3,3),Ei(3,3),Ej(3,3)
	  DOUBLE PRECISION ya,xa,xb,xc,c1,c2,c3,c4
	  DOUBLE PRECISION s1, s2, s3, s4, s5, s6
	  DOUBLE PRECISION small,dif12, dif13, dif23
	  call DX2DX(Xm,DX2)
	  call eye4sym(Is)
	  call zeros(6,6,D)
	  call eigvec2eigproj(V,E)
      small=1.0d-10
	  dif12=abs(x(1)-x(2))
	  dif13=abs(x(1)-x(3))
	  dif23=abs(x(2)-x(3))
! Check eigen values	    
		if ( (dif12 .lt. small) .and. (dif13 .lt. small) .and. 
     & 	 (dif23 .lt. small)) then 
		goto 500
		else if (a .eq. 2) then 
		goto 600
		else if (a .eq. 3) then
		goto 700
		end if
	
	
	
	
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
	
	
500   continue	
      print *, "entrou no 500"
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
      print *, "entrou no 600"
      goto 1000
700   continue
      print *, "entrou no 700"

	  goto 1000
1000   continue

	
      end subroutine
!end code
