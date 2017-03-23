	  subroutine CalcorZero(x,r,Betab,epr,deltax)
! This routine calculates Newton incremet deltax
	  
! Global Variable declaration
	  Double Precision x(5),Betab,epr(3),r(5),deltax(5)
	  
! Local Variable declaration
	  Double Precision K(5,5), Kinv(5,5), detK, dxp(5)
	  integer i, j
	  
! This subroutine needs:
! matinv.f90
!************************************************************************************8	  
	  do i=1,5
		 do j=1,5
			K(i,j)=0.0d0
		 end do
	  end do

      K(1,1)=2.0d0*x(5)
      K(2,2)=2.0d0*x(5)
      K(3,3)=2.0d0*x(5)
	  K(1,4)=1.0d0
	  K(1,5)=2.0d0*x(1)
      K(2,4)=1.0d0
	  K(2,5)=2.0d0*x(2)
	  K(3,4)=1.0d0
	  K(3,5)=2.0d0*x(3)
	  K(4,1)=1.0d0
	  K(4,2)=1.0d0
	  K(4,3)=1.0d0
	  K(5,1)=2.0d0*x(1)
	  K(5,2)=2.0d0*x(2)
	  K(5,3)=2.0d0*x(3)
	  
!	  call det(K,detK)
!	  print *, "detK=",detK
!	  Kinv=K
!	  if (abs(detK) .lt. 0.1d0) then
	     call qr_solve(5,5,K,-r,deltax)
!	  else	 
!         call matinv(Kinv,5)
!		  deltax=matmul(Kinv,-r)
!	  end if	  
!	  dxp=deltax
!	  call vec_print(dxp,5,"deltax")
	  
	  end subroutine
!end code	  