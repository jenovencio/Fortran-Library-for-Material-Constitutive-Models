      subroutine AdSdJ(Cpr,Jn1,Fpn1,dPhidcpr,dSdJ)
! The aim of this routine is to Assemble Derivative of S in relation to J(dSdC)
! based on chain rule:
! dSdJ=(-4/3)*(1/Jn1)^(5/3)*[dW/dCh - 1/3*dW/dCh:Ch*Chinv]
!global variables
      double precision Cpr(3,3),Jn1,Fpn1(3,3),dPhidcpr(3),dSdJ(3,3)

!local variables 
      double precision Ch(3,3),Fpinv(3,3),dWdC(3,3),dW(3,3),
     & k1,k2,Chinv(3,3)
      integer i,j,u,v
      k1=-(4.0d0/3.0d0)*((1.0d0/Jn1)**(5.0d0/3.0d0))	 
      k2=1.0d0/3.0d0
      Fpinv=Fpn1
      call matinv(Fpinv,3)
	  
      Ch=matmul((matmul(transpose(Fpinv),Cpr)),Fpinv)
      call mat_print(Ch,"Ch")
	  Chinv=Ch
      call matinv(Chinv,3)
	  call mat_print(Chinv,"Chinv")
	  
      call zeros(3,3,dW)
      do i=1,3
	     dW(i,i)=dPhidcpr(i)
      end do	 
	  call mat_print(dW,"dW")
	  
      dWdC=matmul((matmul(transpose(Fpinv),dW)),Fpinv)
	  call mat_print(dWdC,"dWdC")
      call zeros(3,3,dSdJ)
	  do i=1,3
	     do j=1,3
		    do u=1,3
			   do v=1,3
		        dSdJ(i,j)=dSdJ(i,j)+dWdC(i,j)
     &				    -k2*(dWdC(u,v)*Ch(u,v)*Chinv(i,j))
	           end do
	        end do
	     end do
      end do		 
      print *,"k1=",k1
      dSdJ=k1*dSdJ
      end subroutine
!end code	  