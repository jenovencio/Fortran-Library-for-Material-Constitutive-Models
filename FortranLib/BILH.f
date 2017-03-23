	  subroutine BILH(ep,prop,sy)
! This routine evaluate BILinear Hardening function
! input  
! ep -> cumulative plastic
! prop -> vector of material properties
	  double precision ep,prop(10),sy
	  sy = prop(1) + prop(2)*ep
	  
	  end subroutine
!end code	  

 