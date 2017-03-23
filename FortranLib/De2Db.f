      subroutine De2Db(b,dtaldee,ksp,dPhidb,d2Phidb)
! the aim of this routine is to calculate d2Phi/dbjdbi based on d2Phi/dejdei:
!where 
! inputs 
! b -> principal values of Left or Right Cauchy Strain
! dtaldee -> derivative of kirshhoff in relation to logarithm strain (second derivative of potential)
! ksp -> kirshhoff principal stress (derivative of potential in relation to logarithm strain)
!
! Output
! d2Phidb -> second derivative of Phi in relation to b
! d2Phi/db(i)db(j) = -1/b(j)^2*dPhide(j)*I + 1/(4bi)*d2Phide(k,i)*dedb(k,j) 

      double precision dtaldee(3,3),ksp(3),b(3),d2Phidb(3,3),dPhidb(3)
      double precision de2db2(3,3,3),dedb(3,3),p1(3,3),p2(3,3)
      integer i,j,k
! derivative of logarithm principal deviatoric strain in relation to bi 	   
      call zeros(3,3,p1)
      call zeros(3,3,p2)
      call zeros(3,3,dedb)
	  
      do i=1,3
	     do j=1,3
		    do k=1,3
			   if ((i .eq. j) .and. (j .eq. k)) then
	              de2db2(k,j,i)=-0.5d0/(b(i)*b(i))
			   else
			       de2db2(k,j,i)=0.0d0
			   end if
               p2(i,j)=p2(i,j)+ksp(k)*de2db2(k,j,i)			   
            end do
         end do			
      end do
	        
	  
      do i=1,3
	     dedb(i,i)=0.5d0/(b(i))
      end do
	  
    
      p1=matmul((matmul(dedb,dtaldee)),dedb)
       
      dPhidb=matmul(ksp,dedb) 
      d2Phidb=p1+p2
	  
      end subroutine
!end code	  