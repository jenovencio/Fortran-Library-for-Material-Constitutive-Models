      subroutine CST(D,J3,B,S,A)
! This routine computes the Consistent Spatial Tangent modulus based on:
! inputs
! D(6,6) - fourth order tensor with derivative of kirshhoff stress DT/DB
! J3 - determinant of F
! B(3,3) - left Cauchy-Green strain tensor
! S(6) - Cauchy Stress
! Output
! A - Consistent Spatial Tangent modulus
! This routine needs:
! eye.f90
! vec2mat.f90
      integer indexix(2,6),indexjy(2,6)
      integer i,j,k,m,l,ix,jy, mid(6,3)
	  DOUBLE PRECISION D(6,6),J3,B(3,3),S(6),Ieye(3,3)
	  DOUBLE PRECISION A(6,6)
	  DOUBLE PRECISION c
	  DOUBLE PRECISION DTB
	  DOUBLE PRECISION CS(3,3)
	
	  call eye(3,3,Ieye)
	  call vec2mat(S,CS)
	  c=2.0d0/J3
! ix = i,j
	  indexix(1,1)=1
	  indexix(1,2)=2
	  indexix(1,3)=3
	  indexix(1,4)=1
	  indexix(1,5)=2
	  indexix(1,6)=1
	  indexix(2,1)=1
	  indexix(2,2)=2
	  indexix(2,3)=3
	  indexix(2,4)=2
	  indexix(2,5)=3
	  indexix(2,6)=3
!jy = k,l
	  indexjy(1,1)=1
	  indexjy(1,2)=2
	  indexjy(1,3)=3
	  indexjy(1,4)=1
	  indexjy(1,5)=2
	  indexjy(1,6)=1
	  indexjy(2,1)=1
	  indexjy(2,2)=2
	  indexjy(2,3)=3
	  indexjy(2,4)=2
	  indexjy(2,5)=3
	  indexjy(2,6)=3

! This is the mapping matrix for multiplication Dijkm*Bml -> summation (m)	  
	  mid(1,1)=1
	  mid(1,2)=4
	  mid(1,3)=6
	  
	  mid(2,1)=4
	  mid(2,2)=2
	  mid(2,3)=5
	  
	  mid(3,1)=6
	  mid(3,2)=5
	  mid(3,3)=3
	  
	  mid(4,1)=1
	  mid(4,2)=4
	  mid(4,3)=6
	  
	  mid(5,1)=4
	  mid(5,2)=2
	  mid(5,3)=5
	  
	  mid(6,1)=1
	  mid(6,2)=4
	  mid(6,3)=6
	  
	  do ix=1,6
		do jy=1,6
		  A(ix,jy)=0.0d0
		end do
	  end do		
	  
	  do ix=1,6
	     do jy=1,6
			!print *, "iteracao i j", ix, jy
				i=indexix(1,ix)
				j=indexix(2,ix)
				k=indexjy(1,jy)
				l=indexjy(2,jy)	
			    DTB=0.0d0
			do m=1,3
				!print *, ix,"DTB part=", D(ix,mid(ix,m))		
				!print *, "B(m,l)=", B(m,l)
				DTB= DTB + D(ix,mid(jy,m))*B(m,l)
			end do
        	    !print *, ix,"DTB=", DTB	
				!print *, jy, "c=", c
				!print *, jy,"DTB=", DTB
				!print *, jy,"CS part =", CS(i,l)
				!print *, jy,"I part =", Ieye(j,k)
				A(ix,jy)=c*DTB - CS(i,l)*Ieye(j,k)
			end do
      end do
	  end subroutine
!end code
