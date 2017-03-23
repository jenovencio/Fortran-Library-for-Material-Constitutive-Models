	  SUBROUTINE eigsym33(M,v,d)
	  INTEGER n,np,nrot,NMAX
      PARAMETER (NMAX=500,n=3,np=3)	 
	  double precision M(np,np),d(np),v(np,np)
	 
* Computes all eigenvalues and eigenvectors of a real symmetric matrix a, which is of size n
* by n, stored in a physical np by np array. On output, elements of a above the diagonal are
*destroyed. d returns the eigenvalues of a in its 1rst n elements. v is a matrix with the same
*logical and physical dimensions as a, whose columns contain, on output, the normalized
*eigenvectors of a. nrot returns the number of Jacobi rotations that were required.	  

! ----------------------------------------------------------------------------
! Parameters:
!   a: The symmetric input matrix
!   v: Storage buffer for eigenvectors
!   d: Storage buffer for eigenvalues
! ----------------------------------------------------------------------------

	  INTEGER i,ip,iq,j
	  double precision c,g,h,s,sm,t,tau,theta,tresh
	  double precision b(NMAX),z(NMAX),a(np,np)
c Store matrix M in matrix a
      a=M	  
	  
c Initialize to the identity matrix.
	  do 12 ip=1,n 
		do 11 iq=1,n
			v(ip,iq)=0.d0
11		continue
        v(ip,ip)=1.d0
12    continue

c Initialize b and d to the diagonal of a.
	  do 13 ip=1,n
		b(ip)=a(ip,ip) 
	    d(ip)=b(ip)
		z(ip)=0.d0
13    continue
	  nrot=0
	  
	  
	  do 24 i=1,50
			sm=0.0d0
		do 15 ip=1,n-1
			do 14 iq=ip+1,n
				sm=sm+abs(a(ip,iq))
14			continue
15      continue 

		if (sm .eq. 0.0d0) return
		if (i.lt.4) then
			tresh=0.2d0*sm/n**2.0d0 
		else
			tresh=0.0d0 
		endif
		
		do 22 ip=1,n-1
			do 21 iq=ip+1,n
				g=100.d0*abs(a(ip,iq))
				if((i.gt.4).and.(abs(d(ip))+g.eq.abs(d(ip)))
     & .and.(abs(d(iq))+g.eq.abs(d(iq))))then
					a(ip,iq)=0.d0
				else if(abs(a(ip,iq)).gt.tresh)then
					h=d(iq)-d(ip)
					if(abs(h)+g.eq.abs(h))then	
					t=a(ip,iq)/h
				else
					theta=0.50d0*h/a(ip,iq)
					t=1.d0/(abs(theta)+sqrt(1.0d0+theta**2.0d0))
					if(theta.lt.0.0d0) t=-t
                endif
				c=1./sqrt(1+t**2)
				s=t*c
				tau=s/(1.+c)
				h=t*a(ip,iq)
				z(ip)=z(ip)-h
				z(iq)=z(iq)+h
				d(ip)=d(ip)-h
				d(iq)=d(iq)+h
				a(ip,iq)=0.0d0
				do 16 j=1,ip-1 
					g=a(j,ip)
					h=a(j,iq)
					a(j,ip)=g-s*(h+g*tau)
					a(j,iq)=h+s*(g-h*tau)
16				continue
				do 17 j=ip+1,iq-1 
					g=a(ip,j)
					h=a(j,iq)
					a(ip,j)=g-s*(h+g*tau)
					a(j,iq)=h+s*(g-h*tau)
17				continue
				do 18 j=iq+1,n 
					g=a(ip,j)
					h=a(iq,j)
					a(ip,j)=g-s*(h+g*tau)
					a(iq,j)=h+s*(g-h*tau)
18              continue
				do 19 j=1,n
					g=v(j,ip)
					h=v(j,iq)
					v(j,ip)=g-s*(h+g*tau)
					v(j,iq)=h+s*(g-h*tau)
19              continue
				nrot=nrot+1
			endif
21      continue
22      continue
		do 23 ip=1,n
			b(ip)=b(ip)+z(ip)
			d(ip)=b(ip) 
			z(ip)=0.0d0
23      continue
24    continue
	  print *, 'too many iterations in jacobi'
	  return
	  END SUBROUTINE
!end code	  		  