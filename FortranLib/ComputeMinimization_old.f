	  subroutine ComputeMinimization(epr,Betab,Taub,dt,Qn,q0,q,
     &	  deltaQ,nor)
	  DOUBLE PRECISION q0(3), nor, deltaQ,dt
	  DOUBLE PRECISION epr(3),Betab,Taub,Qn,q(3)
	  ! Local variables
	  DOUBLE PRECISION Lam, bet,deltaQmin, x(6),r(6), d(6) 
	  DOUBLE PRECISION step, xt(6),maxnor
	  integer it,conv, pau, maxint
! This subrotine needs:
! Calres.f
! Calcor.f
	  
!     Subroutine parameters
	  maxint=60
	  maxnor=1.0d-6
	  conv = 0
	  it = 1
	  Lam = 1.0d-6
	  bet= 1.0d-6
	  deltaQmin = 1.0d-8
	  deltaQ = deltaQmin
	  nor =1.0d0
	  
	  !---------------------------------------------------
	  x(1)=q(1) 
	  x(2)=q(2) 
	  x(3)=q(3) 
	  x(4)=deltaQ 
	  x(5)=Lam 
	  x(6)=bet
	  
	  
	  do it=1,maxint
	    if (nor .gt. maxnor)  then
		  call Calres(x,Qn,Betab,Taub,epr,dt,r)
		  nor = dnrm2(6,r,1)
	      call Calcor(x,dt,r,Qn,Betab,Taub,epr,d)

		  step=1.0d0
		  do i=1,100
			  xt=x+step*d
!             call vec_print(xt,6,"xt")
			  if (xt(4) .gt. deltaQmin) then
				  goto 100
			  else
				  step = step*8.0d-1
			  end if
		  end do	
100       if (i .eq. 100) then
			pau = 1
		  end if
		  x=xt
		  
		  call Calres(x,Qn,Betab,Taub,epr,dt,r)
		  nor = dnrm2(6,r,1)
		  
	    else
            goto 200
		end if		
	  end do
200   continue

! Updating output variables
	  q(1)=x(1)
	  q(2)=x(2)
	  q(3)=x(3)
      deltaQ=x(4)
      Lam=x(5)
      bet=x(6)
	  
	  return
	  end subroutine
!end code	  