      subroutine dEdevdE(M)
! this subroutine create a matrix with derivatives of Eiso to E
!  eiso(i)=	e(i)-(1/3)*trace(e)
! where e are the principal logarithm strain
! and eiso are the principal isochoric logarithm strain
!  M(i,j)=deiso(i)/de(j)      
      double precision M(3,3)
      
	  M(1,1)=2.0d0/3.0d0
	  M(1,2)=-1.0d0/3.0d0
	  M(1,3)=-1.0d0/3.0d0
	  
	  M(2,1)=-1.0d0/3.0d0
	  M(2,2)=2.0d0/3.0d0
	  M(2,3)=-1.0d0/3.0d0
	  
	  M(3,1)=-1.0d0/3.0d0
	  M(3,2)=-1.0d0/3.0d0
	  M(3,3)=2.0d0/3.0d0
	  
	  
	  
      end subroutine
!end code	  