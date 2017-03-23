      subroutine Acmatvol(Cn1,Jn1,Cmatvol)
! This routine contains the Volumetric part of Spatial Tagent Matrix  	  
      double precision Jn1,Cmatvol(6,6),Mvol(6,6),dU,
     & Cn1(3,3),Cn1inv(3,3),dJCdC(3,3,3,3)	  
      integer i,j,k,l
      common /MatTP/ TP
	  	  
      call dVolu(Jn1,dU)
	  Cn1inv=Cn1
      call matinv(Cn1inv,3)	  
      do i=1,3
	     do j=1,3
		    do k=1,3
			   do l=1,3
			      dJCdC(i,j,k,l)=Jn1*(Cn1inv(i,j)*Cn1inv(k,l)
     &                         -2.0d0*Cn1inv(i,k)*Cn1inv(j,l))
               end do
            end do
         end do			
      end do		 

      call tensor2matrix(dJCdC,Mvol)	  
      Cmatvol=dU*Mvol
      end subroutine
!end code	  