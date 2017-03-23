      subroutine cspa2cjau(cspa,Cauchy,cjau)
! the aim of this routine is to convert Spacial Tangent Matrix (cspa)
! in to Jauman Tangent Matrix cjau
! This routine is based on the following equation:
! cjau(i,j,k,l)=cspa(i,j,k,l)+0.5(I(i,k)*C(j,l)+I(j,l)*C(i,k)+I(i,l)*C(j,k)+I(j,k)*C(i,l))
! global varibles      
      double precision cspa(6,6),Cauchy(3,3),cjau(6,6)
      double precision SI(6,6)
      integer i,j

! SI = 	0.5(I(i,k)*C(j,l)+I(j,l)*C(i,k)+I(i,l)*C(j,k)+I(j,k)*C(i,l))
       call zeros(6,6,SI)
	   SI(1,1)=2.0D0*Cauchy(1,1) 
	   SI(2,2)=2.0D0*Cauchy(2,2) 
	   SI(3,3)=2.0D0*Cauchy(3,3)
	   SI(4,4)=(Cauchy(1,1)+Cauchy(2,2))/2.0D0
	   SI(5,5)=(Cauchy(2,2)+Cauchy(3,3))/2.0D0
	   SI(6,6)=(Cauchy(1,1)+Cauchy(3,3))/2.0D0
	   
	   SI(1,4)=Cauchy(1,2)
	   SI(1,6)=Cauchy(1,3)
	   
	   SI(2,4)=Cauchy(1,2)
       SI(2,5)=Cauchy(2,3)
	   
	   SI(3,5)=Cauchy(2,3)
	   SI(3,6)=Cauchy(1,3)
	   
	   SI(4,5)=Cauchy(1,3)/2.0D0
	   SI(4,6)=Cauchy(2,3)/2.0D0
	   
	   SI(5,6)=Cauchy(1,2)/2.0D0
	   
      do i=1,6
	      do j=i,6
		     SI(j,i)=SI(i,j)
          end do
      end do		  
	   
	  cjau= cspa + SI

      end subroutine
!end code	  