      subroutine Adsde(dtaldee,ks,Jn1,dsidej)
! The aim of this routine is to Assemble Deviatoric + Volumetric  2nd Derivative (e1,e2,e3)	  
! where e1,e2,e3 are the principal values of logarithm strain
! edev(i) = e(i)-(1/3)*(e1+e2+e3) 
! dsde(i,j) = dJ^(-1)/de(j)*ks(i)+1/J*dtaldee

! global varible declaration
      double precision dtaldee(3,3),Jn1,dsidej(3,3)

! local variable declaration	  
      double precision dJinv(3),ks(3),dJinvks(3,3)
      integer i,j

	  

      do i=1,3
         do j=1,3
		    dJinv(j)=-1.0d0/Jn1
	        dJinvks(i,j)=dJinv(j)*ks(i)
			dsidej(j,i)=dJinvks(i,j)+(1.d0/Jn1)*dtaldee(i,j)
         end do
      end do
	  
! printing to debug	  
!      print *, "J=",J
!      print *, "dU=",dU
!	  call vec_print(d2U,3,"d2U")
!	  call vec_print(dJ,3,"dJ")
!	  call mat_print(dJinvks,"dJinvks")
!	  call mat_print(dE,"dE")
!	  call mat_print(d2PdEi,"d2PdEi")
	  
!---------------------------------------------------------


      end subroutine
!end code	