      SUBROUTINE matinv(a,n)
	  INTEGER n,np,NMAX
	  double precision a(n,n)
	  PARAMETER (NMAX=50)
	  INTEGER i,icol,irow,j,k,l,ll,indxc(NMAX),indxr(NMAX), ipiv(NMAX)

	  do j=1,n
		  ipiv(j)=0
	  end do

	  do i=1,n
		  big=0
		  do j=1,n 
			  if (ipiv(j) .ne. 1)then
				  do k=1,n
					  if (ipiv(k) .eq. 0) then
						  if (abs(a(j,k)) .ge. big)then
							  big=abs(a(j,k))
							  irow=j
							  icol=k
						  end if
					  else if (ipiv(k) .gt. 1) then
						      goto 500
					  end if
				  end do
			  end if
		  end do

	    ipiv(icol)=ipiv(icol)+1

		if (irow .ne. icol) then
			do  l=1,n
				dum=a(irow,l)
				a(irow,l)=a(icol,l)
				a(icol,l)=dum
			end do 
		end if

		indxr(i)=irow 
		indxc(i)=icol 
		if (a(icol,icol) .eq. 0.0d0) goto 600
			pivinv=1.0d0/a(icol,icol)
			a(icol,icol)=1.0d0
		do l=1,n
			a(icol,l)=a(icol,l)*pivinv
		end do
		
		do ll=1,n 
			if(ll .ne. icol)then 
				dum=a(ll,icol)
				a(ll,icol)=0.0d0
				do l=1,n
					a(ll,l)=a(ll,l)-a(icol,l)*dum
				end do
			end if
		end do
	  end do  
	 
	  do l=n,1,-1 
		if(indxr(l).ne.indxc(l))then
			do k=1,n
				dum=a(k,indxr(l))
				a(k,indxr(l))=a(k,indxc(l))
				a(k,indxc(l))=dum
			end do
		end if
	  end do
      return 
     
500   continue
   	   print *, 'singular matrix in gaussj - error 1'
	   goto 700
	 
600   continue
   	 print *, 'singular matrix in gaussj - error 2'
	 
700   continue 
	  end SUBROUTINE
!end code  