      subroutine AdSdC(Cn1,Ch,Jn1,dWdCh,d2WdCh,dSdC)
! The aim of this routine is to Assemble Material Tagente Matrix (dSdC)
! based on: 	  
!dSdCh - > derivative of Second Piola Stress in relation to C
!C -> right Cauchy Strain 
!Jn1 -> determinant of F
! dPhidCh -> rerivative of potential in relation to isochoric left Cauchy Strain
! dSdC = 2*dSdCh*[H] + Jn1*dS/dJ*C^-1
! H= (1/Jn1)^(2/3)*[0.5(I(m,k)*I(n,l) + I(m,l)*I(n,k) )- (1/3)*Ch(m,n)*Ch(k,l)^-1] 
! global variables      
      double precision Cn1(3,3),Cpr(3,3),Jn1,dWdCh(3,3),
     &  d2WdCh(6,6),Fpn1(3,3),dSdC(6,6)
	  
! local variables
      double precision dSdCh(3,3,3,3), dSdC4t(3,3,3,3),H(3,3,3,3),
     & Ident(3,3),Ch(3,3),Cinv(3,3),Chinv(3,3),
     & Fpinv(3,3),T1(3,3,3,3),T2(3,3,3,3),T3(3,3,3,3),dW,Zt(6,6),
     & d2WdCh4t(3,3,3,3),dSdJ(3,3)
      double precision k1,k2,k3,k4
      double precision HM(6,6),dSdChM(6,6)
      integer i,j,m,n,k,l,u,v
! Basic Variables      
      k1=(1.0d0/Jn1)**(2.0d0/3.0d0)
      k2=1.0d0/3.0d0
      k3=(-4.0d0/3.0d0)*((1.0d0/Jn1)**(5.0d0/3.0d0))

      Cinv=Cn1
      Chinv=Ch
      call matinv(Cinv,3)
      call matinv(Chinv,3)
      call eye(3,3,Ident)
      call zeros(6,6,Zt)
      call matrix2tensor(Zt,dSdC4t)
	  call matrix2tensor(Zt,T1)
	  call matrix2tensor(Zt,T2)
	  call matrix2tensor(Zt,T3)
	  call matrix2tensor(d2WdCh,d2WdCh4t)
	  
!	  dWdCh=0.50d0*dWdCh
!********************************************************************************
!                           dSdCh Calculation
!    dSdCh = 2*k1(T1-1/3*T2) +  1/3*k1*dW*T3
!     T1(i,j,m,n) = d2WdCh -1/3*d2WdCh*Ch*C^-1
!     T2(i,j,m,n) = dWdCh*Ch^-1
!     dw = dWdCh:Ch^-1
!     T3 = Ch^-1(i,m)Ch^-1(n,j)+Ch^-1(i,n)Ch^-1(m,j)
!********************************************************************************
      dW=0.0d0
      do u=1,3
	     do v=1,3
	        dW=dW+dWdCh(u,v)*Ch(u,v)
         end do
      end do

      do i=1,3
	     do j=1,3
		    do m=1,3
			   do n=1,3
			       do u=1,3
				       do v=1,3
					     T1(i,j,m,n)=T1(i,j,m,n)
     &                   -d2WdCh4t(u,v,m,n)*Ch(u,v)*Chinv(i,j)
                       end do
                    end do
					T1(i,j,m,n)=k2*T1(i,j,m,n) + d2WdCh4t(i,j,m,n)
                    T2(i,j,m,n)=dWdCh(m,n)*Chinv(i,j)
                    T3(i,j,m,n)=Chinv(i,m)*Chinv(n,j)
     &					         +Chinv(i,n)*Chinv(m,j)
	               dSdCh(i,j,m,n)=2.0d0*k1*T1(i,j,m,n)
     &                           -2.0d0*k1*k2*T2(i,j,m,n)
     &				              +k1*k2*dW*T3(i,j,m,n)	  
                end do
            end do
         end do
      end do		 
	  

  
!********************************************************************************
!                  H(i,j,k,l) Calculation
!H=(1/Jn1)^(2/3)*[0.5(I(m,k)*I(n,l) + I(m,l)*I(n,k) )- (1/3)*Ch(m,n)*Ch(k,l)^-1] 
!******************************************************************************** 	  
      do m=1,3
	     do n=1,3
		    do k=1,3
			   do l=1,3
			       H(m,n,k,l)=0.5d0*(Ident(m,k)*Ident(n,l) 
     &			+ Ident(m,l)*Ident(n,k))-(k2)*Ch(m,n)*Chinv(k,l)
	            end do
            end do	
        end do
      end do
      H=k1*H


!********************************************************************************

!********************************************************************************
!                             dS/dJ Calculation
!             dSdJ=(-4/3)*(1/Jn1)^(5/3)*[dW/dCh - 1/3*dW/dCh:Ch*Chinv]
!********************************************************************************
      call zeros(3,3,dSdJ)
      k4=0.0d0
      do u=1,3
	     do v=1,3
		        k4=k4+dWdCh(u,v)*Ch(u,v)
	     end do
      end do
	  
      do i=1,3
	     do j=1,3
			dSdJ(i,j)=k3*(dWdCh(i,j)-k2*k4*Chinv(i,j))
	     end do
      end do	

!********************************************************************************




!********************************************************************************
!                             dSdC Calculation
!                    dSdC = 2*dSdCh*[H] + Jn1*dS/dJ*C^-1
!********************************************************************************

	  do i=1,3
	     do j=1,3
		    do k=1,3
			   do l=1,3
			      do m=1,3
				     do n=1,3
			            dSdC4t(i,j,k,l)=dSdC4t(i,j,k,l)+
     & 	dSdCh(i,j,m,n)*H(m,n,k,l) 
	                 end do
				  end do
                 dSdC4t(i,j,k,l)=2.0d0*dSdC4t(i,j,k,l)
     &                        +Jn1*dSdJ(i,j)*Cinv(k,l)		
	            end do
            end do	
        end do
      end do

!********************************************************************************

! Convert 4ord tensor into Voight Notation	  
	  dSdC4t=(1.0d0/4.0d0)*dSdC4t
	  call tensor2matrix(dSdC4t,dSdC)
      end subroutine
!end code