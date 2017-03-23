*deck,usermat      USERDISTRIB  parallel                                gal
      subroutine usermat(
     &                   matId, elemId,kDomIntPt, kLayer, kSectPt,
     &                   ldstep,isubst,keycut,
     &                   nDirect,nShear,ncomp,nStatev,nProp,
     &                   Time,dTime,Temp,dTemp,
     &                   stress,ustatev,dsdePl,sedEl,sedPl,epseq,
     &                   Strain,dStrain, epsPl, prop, coords, 
     &                   var0, defGrad_t, defGrad,
     &                   tsstif, epsZZ,
     &                   var1, var2, var3, var4, var5,
     &                   var6, var7, var8)
c*************************************************************************
c     *** primary function ***
c
c           user defined material constitutive model
c
c      Attention:
c           User must define material constitutive law properly
c           according to the stress state such as 3D, plane strain
c           and axisymmetry, plane stress and 3D/1D beam.
c
c           A 3D material constitutive model can be used for
c           plane strain and axisymmetry cases.
c
c           When using shell elements, a plane stress algorithm
c           must be used.
c
c                                             gal July, 1999
c
c       The following demonstrates a USERMAT subroutine for 
c       a plasticity model, which is the same as TB, BISO,
c       for different stress states. 
c       See "ANSYS user material subroutine USERMAT" for detailed
c       description of how to write a USERMAT routine.
c
c       This routine calls four routines,
c       usermat3d.F, usermatps.F usermatbm.F and usermat1d.F, w.r.t.
c       the corresponding stress states.
c       Each routine can be also a usermat routine for the specific 
c       element.
c
c*************************************************************************
c
c     input arguments
c     ===============
c      matId     (int,sc,i)               material #
c      elemId    (int,sc,i)               element #
c      kDomIntPt (int,sc,i)               "k"th domain integration point
c      kLayer    (int,sc,i)               "k"th layer
c      kSectPt   (int,sc,i)               "k"th Section point
c      ldstep    (int,sc,i)               load step number
c      isubst    (int,sc,i)               substep number
c      nDirect   (int,sc,in)              # of direct components
c      nShear    (int,sc,in)              # of shear components
c      ncomp     (int,sc,in)              nDirect + nShear
c      nstatev   (int,sc,i)               Number of state variables
c      nProp     (int,sc,i)               Number of material constants
c
c      Temp      (dp,sc,in)               temperature at beginning of
c                                         time increment
c      dTemp     (dp,sc,in)               temperature increment 
c      Time      (dp,sc,in)               time at beginning of increment (t)
c      dTime     (dp,sc,in)               current time increment (dt)
c
c      Strain   (dp,ar(ncomp),i)          Strain at beginning of time increment
c      dStrain  (dp,ar(ncomp),i)          Strain increment
c      prop     (dp,ar(nprop),i)          Material constants defined by TB,USER
c      coords   (dp,ar(3),i)              current coordinates
c      defGrad_t(dp,ar(3,3),i)            Deformation gradient at time t
c      defGrad  (dp,ar(3,3),i)            Deformation gradient at time t+dt
c
c     input output arguments              
c     ======================             
c      stress   (dp,ar(ncomp),io)         stress
c      ustatev   (dp,ar(nstatev),io)      user state variables
c      sedEl    (dp,sc,io)                elastic work
c      sedPl    (dp,sc,io)                plastic work
c      epseq    (dp,sc,io)                equivalent plastic strain
c      epsPl   (dp,ar(ncomp),io)          plastic strain
c      var?     (dp,sc,io)                not used, they are reserved arguments 
c                                         for further development
c
c     output arguments
c     ================
c      keycut   (int,sc,o)                loading bisect/cut control
c                                         0 - no bisect/cut
c                                         1 - bisect/cut 
c                                         (factor will be determined by solution control)
c      dsdePl   (dp,ar(ncomp,ncomp),o)    material jacobian matrix
c      tsstif   (dp,ar(2),o)              transverse shear stiffness
c                                         tsstif(1) - Gxz
c                                         tsstif(2) - Gyz
c                                         tsstif(1) is also used to calculate hourglass
c                                         stiffness, this value must be defined when low
c                                         order element, such as 181, 182, 185 with uniform 
c                                         integration is used.
c      epsZZ    (dp,sc,o)                 strain epsZZ for plane stress,
c                                         define it when accounting for thickness change
c                                         in shell and plane stress states
c
c*************************************************************************
c
c      ncomp   6   for 3D  (nshear=3)
c      ncomp   4   for plane strain or axisymmetric (nShear = 1)
c      ncomp   3   for plane stress (nShear = 1)
c      ncomp   3   for 3d beam      (nShear = 2)
c      ncomp   1   for 1D (nShear = 0)
c
c      stresses and strains, plastic strain vectors
c          11, 22, 33, 12, 23, 13    for 3D
c          11, 22, 33, 12            for plane strain or axisymmetry
c          11, 22, 12                for plane stress
c          11, 13, 12                for 3d beam
c          11                        for 1D
c
c      material jacobian matrix
c        3D
c           dsdePl    |  1111   1122   1133   1112   1123   1113 |
c           dsdePl    |  2211   2222   2233   2212   2223   2213 |
c           dsdePl    |  3311   3322   3333   3312   3323   3313 |
c           dsdePl    |  1211   1222   1233   1212   1223   1213 |
c           dsdePl    |  2311   2322   2333   2312   2323   2313 |
c           dsdePl    |  1311   1322   1333   1312   1323   1313 |
c        plane strain or axisymmetric (11, 22, 33, 12)
c           dsdePl    |  1111   1122   1133   1112 |
c           dsdePl    |  2211   2222   2233   2212 |
c           dsdePl    |  3311   3322   3333   3312 |
c           dsdePl    |  1211   1222   1233   1212 |
c        plane stress (11, 22, 12)
c           dsdePl    |  1111   1122   1112 |
c           dsdePl    |  2211   2222   2212 |
c           dsdePl    |  1211   1222   1212 |
c        3d beam (11, 13, 12)
c           dsdePl    |  1111   1113   1112 |
c           dsdePl    |  1311   1313   1312 |
c           dsdePl    |  1211   1213   1212 |
c        1d
c           dsdePl    |  1111 |
c
c*************************************************************************
#include "impcom.inc"
c
      INTEGER          
     &                 matId, elemId,
     &                 kDomIntPt, kLayer, kSectPt,
     &                 ldstep,isubst,keycut,
     &                 nDirect,nShear,ncomp,nStatev,nProp
      DOUBLE PRECISION 
     &                 Time,    dTime,   Temp,    dTemp,
     &                 sedEl,   sedPl,   epseq,   epsZZ
      DOUBLE PRECISION 
     &                 stress  (ncomp  ), ustatev (nStatev),
     &                 dsdePl  (ncomp,ncomp),
     &                 Strain  (ncomp  ), dStrain (ncomp  ), 
     &                 epsPl   (ncomp  ), prop    (nProp  ), 
     &                 coords  (3),       
     &                 defGrad (3,3),     defGrad_t(3,3),
     &                 tsstif  (2)
c
      EXTERNAL         usermat3d, usermatps, usermatbm, usermat1d

      DOUBLE PRECISION var0, var1, var2, var3, var4, var5,
     &                 var6, var7, var8
c
c*************************************************************************
c
      IF(ncomp .GE. 4) THEN
c ***    3d, plane strain and axisymmetric example
         call usermat3d (
     &                   matId, elemId,kDomIntPt, kLayer, kSectPt,
     &                   ldstep,isubst,keycut,
     &                   nDirect,nShear,ncomp,nStatev,nProp,
     &                   Time,dTime,Temp,dTemp,
     &                   stress,ustatev,dsdePl,sedEl,sedPl,epseq,
     &                   Strain,dStrain, epsPl, prop, coords,
     &                   var0, defGrad_t, defGrad,
     &                   tsstif, epsZZ,
     &                   var1, var2, var3, var4, var5,
     &                   var6, var7, var8)

      ELSE IF(nDirect.eq. 2 .and. ncomp .EQ. 3) THEN
c ***    plane stress example
         call usermatps (
     &                   matId, elemId,kDomIntPt, kLayer, kSectPt,
     &                   ldstep,isubst,keycut,
     &                   nDirect,nShear,ncomp,nStatev,nProp,
     &                   Time,dTime,Temp,dTemp,
     &                   stress,ustatev,dsdePl,sedEl,sedPl,epseq,
     &                   Strain,dStrain, epsPl, prop, coords,
     &                   var0, defGrad_t, defGrad,
     &                   tsstif, epsZZ,
     &                   var1, var2, var3, var4, var5,
     &                   var6, var7, var8)

      ELSE IF(ncomp .EQ. 3) THEN
c ***    3d beam example
         call usermatbm (
     &                   matId, elemId,kDomIntPt, kLayer, kSectPt,
     &                   ldstep,isubst,keycut,
     &                   nDirect,nShear,ncomp,nStatev,nProp,
     &                   Time,dTime,Temp,dTemp,
     &                   stress,ustatev,dsdePl,sedEl,sedPl,epseq,
     &                   Strain,dStrain, epsPl, prop, coords,
     &                   var0, defGrad_t, defGrad,
     &                   tsstif, epsZZ,
     &                   var1, var2, var3, var4, var5,
     &                   var6, var7, var8)

      ELSE IF(ncomp .EQ. 1) THEN
c ***    1d beam example
         call usermat1d (
     &                   matId, elemId,kDomIntPt, kLayer, kSectPt,
     &                   ldstep,isubst,keycut,
     &                   nDirect,nShear,ncomp,nStatev,nProp,
     &                   Time,dTime,Temp,dTemp,
     &                   stress,ustatev,dsdePl,sedEl,sedPl,epseq,
     &                   Strain,dStrain, epsPl, prop, coords,
     &                   var0, defGrad_t, defGrad,
     &                   tsstif, epsZZ,
     &                   var1, var2, var3, var4, var5,
     &                   var6, var7, var8)

      END IF
      return
      end
*deck,usermat1d    USERDISTRIB  parallel                                gal
      subroutine usermat1d(
     &                   matId, elemId,kDomIntPt, kLayer, kSectPt,
     &                   ldstep,isubst,keycut,
     &                   nDirect,nShear,ncomp,nStatev,nProp,
     &                   Time,dTime,Temp,dTemp,
     &                   stress,ustatev,dsdePl,sedEl,sedPl,epseq,
     &                   Strain,dStrain, epsPl, prop, coords, 
     &                   var0, defGrad_t, defGrad,
     &                   tsstif, epsZZ,
     &                   var1, var2, var3, var4, var5,
     &                   var6, var7, var8)
c*************************************************************************
c     *** primary function ***
c
c           user defined material constitutive model
c
c      Attention:
c           User must define material constitutive law properly
c           according to the stress state such as 3D, plane strain
c           and axisymmetry, plane stress and beam.
c
c           a 3D material constitutive model can use for
c           plane strain and axisymmetry cases.
c
c           When using shell elements, a plane stress algorithm
c           must be use.
c
c                                             gal July, 1999
c
c       The following demonstrates a USERMAT subroutine for
c       a plasticity model of 1D truss element (LINK180). 
c       The plasticity model is the same as TB, BISO.
c
c       See "ANSYS user material subroutine USERMAT" for detailed
c       description of how to write a USERMAT routine.
c
c*************************************************************************
c
c     input arguments
c     ===============
c      matId     (int,sc,i)               material #
c      elemId    (int,sc,i)               element #
c      kDomIntPt (int,sc,i)               "k"th domain integration point
c      kLayer    (int,sc,i)               "k"th layer
c      kSectPt   (int,sc,i)               "k"th Section point
c      ldstep    (int,sc,i)               load step number
c      isubst    (int,sc,i)               substep number
c      nDirect   (int,sc,in)              # of direct components
c      nShear    (int,sc,in)              # of shear components
c      ncomp     (int,sc,in)              nDirect + nShear
c      nstatev   (int,sc,l)               Number of state variables
c      nProp     (int,sc,l)               Number of material ocnstants
c
c      Temp      (dp,sc,in)               temperature at beginning of
c                                         time increment
c      dTemp     (dp,sc,in)               temperature increment 
c      Time      (dp,sc,in)               time at beginning of increment (t)
c      dTime     (dp,sc,in)               current time increment (dt)
c
c      Strain   (dp,ar(ncomp),i)          Strain at beginning of time increment
c      dStrain  (dp,ar(ncomp),i)          Strain increment
c      prop     (dp,ar(nprop),i)          Material constants defined by TB,USER
c      coords   (dp,ar(3),i)              current coordinates
c      defGrad_t(dp,ar(3,3),i)            Deformation gradient at time t
c      defGrad  (dp,ar(3,3),i)            Deformation gradient at time t+dt
c
c     input output arguments              
c     ======================             
c      stress   (dp,ar(nTesn),io)         stress
c      ustatev   (dp,ar(nstatev),io)      user state variables
c            ustatev(1)                     - equivalent plastic strain
c            ustatev(2) - ustatev(1+ncomp)   - plastic strain vector
c            ustatev(nStatev)               - von-Mises stress
c      sedEl    (dp,sc,io)                elastic work
c      sedPl    (dp,sc,io)                plastic work
c      epseq    (dp,sc,io)                equivalent plastic strain
c      tsstif   (dp,ar(2),io)             transverse shear stiffness
c                                         tsstif(1) - Gxz
c                                         tsstif(2) - Gyz
c                                         tsstif(1) is also used to calculate hourglass
c                                         stiffness, this value must be defined when low
c                                         order element, such as 181, 182, 185 with uniform 
c                                         integration is used.
c      var?     (dp,sc,io)                not used, they are reserved arguments 
c                                         for further development
c
c     output arguments
c     ================
c      keycut   (int,sc,io)               loading bisect/cut control
c                                         0 - no bisect/cut
c                                         1 - bisect/cut 
c                                         (factor will be determined by ANSYS solution control)
c      dsdePl   (dp,ar(ncomp,ncomp),io)   material jacobian matrix
c      epsZZ    (dp,sc,o)                 strain epsZZ for plane stress,
c                                         define it when accounting for thickness change 
c                                         in shell and plane stress states
c
c*************************************************************************
c
c      ncomp   6   for 3D  (nshear=3)
c      ncomp   4   for plane strain or axisymmetric (nShear = 1)
c      ncomp   3   for plane stress (nShear = 1)
c      ncomp   3   for 3d beam      (nShear = 2)
c      ncomp   1   for 1D (nShear = 0)
c
c      stresss and strains, plastic strain vectors
c          11, 22, 33, 12, 23, 13    for 3D
c          11, 22, 33, 12            for plane strain or axisymmetry
c          11, 22, 12                for plane stress
c          11, 13, 12                for 3d beam
c          11                        for 1D
c
c      material jacobian matrix
c        3D
c           dsdePl    |  1111   1122   1133   1112   1123   1113 |
c           dsdePl    |  2211   2222   2233   2212   2223   2213 |
c           dsdePl    |  3311   3322   3333   3312   3323   3313 |
c           dsdePl    |  1211   1222   1233   1212   1223   1213 |
c           dsdePl    |  2311   2322   2333   2312   2323   2313 |
c           dsdePl    |  1311   1322   1333   1312   1323   1313 |
c        plane strain or axisymmetric (11, 22, 33, 12)
c           dsdePl    |  1111   1122   1133   1112 |
c           dsdePl    |  2211   2222   2233   2212 |
c           dsdePl    |  3311   3322   3333   3312 |
c           dsdePl    |  1211   1222   1233   1212 |
c        plane stress (11, 22, 12)
c           dsdePl    |  1111   1122   1112 |
c           dsdePl    |  2211   2222   2212 |
c           dsdePl    |  1211   1222   1212 |
c        3d beam (11, 13, 12)
c           dsdePl    |  1111   1113   1112 |
c           dsdePl    |  1311   1313   1312 |
c           dsdePl    |  1211   1213   1212 |
c        1d
c           dsdePl    |  1111 |
c
c*************************************************************************
#include "impcom.inc"
c
      INTEGER          
     &                 matId, elemId,
     &                 kDomIntPt, kLayer, kSectPt,
     &                 ldstep,isubst,keycut,
     &                 nDirect,nShear,ncomp,nStatev,nProp
      DOUBLE PRECISION 
     &                 Time,    dTime,   Temp,    dTemp,
     &                 sedEl,   sedPl,   epseq,   epsZZ
      DOUBLE PRECISION 
     &                 stress  (ncomp  ), ustatev (nStatev),
     &                 dsdePl  (ncomp,ncomp),
     &                 Strain  (ncomp  ), dStrain (ncomp  ), 
     &                 epsPl   (ncomp  ), prop    (nProp  ), 
     &                 coords  (3),
     &                 defGrad (3,3),     defGrad_t(3,3),
     &                 tsstif  (2)
c
c***************** User defined part *************************************
c
c --- parameters
c
      INTEGER          mcomp
      DOUBLE PRECISION ZERO, HALF, ONE, TWO, SMALL
      PARAMETER       (ZERO       = 0.d0,
     &                 HALF       = 0.5d0,
     &                 ONE        = 1.d0,
     &                 TWO        = 2.d0,
     &                 SMALL      = 1.d-08,
     &                 mcomp      = 1
     &                 )
c
c --- local variables
c
c      sigElp   (dp,ar(6  ),l)            trial stress
c      dsdeEl   (dp,ar(6,6),l)            elastic moduli
c      sigDev   (dp,ar(6  ),l)            deviatoric stress tensor
c      dfds     (dp,ar(6  ),l)            derivative of the yield function 
c      JM       (dp,ar(6,6),l)            2D matrix for a 4 order tensor
c      pEl      (dp,sc     ,l)            hydrostatic pressure stress
c      qEl      (dp,sc     ,l)            von-mises stress
c      pleq_t   (dp,sc     ,l)            equivalent plastic strain at beginnig of time increment
c      pleq     (dp,sc     ,l)            equivalent plastic strain at end of time increment
c      dpleq    (dp,sc     ,l)            incremental equivalent plastic strain
c      sigy_t   (dp,sc     ,l)            yield stress at beginnig of time increments
c      sigy     (dp,sc     ,l)            yield stress at end of time increment
c      young    (dp,sc     ,l)            Young's modulus
c      posn     (dp,sc     ,l)            Poiss's ratio
c      sigy0    (dp,sc     ,l)            initial yield stress
c      dsigdep  (dp,sc     ,l)            plastic slop
c      twoG     (dp,sc     ,l)            two time of shear moduli
c
c
      DOUBLE PRECISION sigElp(mcomp), dsdeEl(mcomp,mcomp)

      DOUBLE PRECISION var0, var1, var2, var3, var4, var5,
     &                 var6, var7, var8

      DOUBLE PRECISION qEl,   pleq_t,  sigy_t , sigy,
     &                 dpleq, pleq,    signTens,
     &                 young, posn,    sigy0,   dsigdep, 
     &                 twoG,  fratio
c*************************************************************************
c
      keycut   = 0
      dsigdep  = ZERO 
      pleq_t   = ustatev(1)
      pleq     = pleq_t
c *** get Young's modulus and Poisson's ratio, initial yield stress and others
      young    = prop(1)
      posn     = prop(2)
      sigy0    = prop(3)
c *** calculate plastic slope
      dsigdep  = young*prop(4)/(young-prop(4))
      twoG     = young / (ONE+posn)
c *** define tsstif(1) since it is used for calculation of hourglass stiffness
      tsstif(1) = HALF * twoG
c
c *** calculate elastic stiffness matrix 
c
      dsdeEl(1,1)= young
c
c *** calculate the trial stress and 
c     copy elastic moduli dsdeEl to material Jacobian matrix
      sigElp(1)   = stress(1)
      dsdePl(1,1) = dsdeEl(1,1)
      sigElp(1)   = sigElp(1) + dsdeEl(1,1) * dStrain(1)
c *** sign of predicted stress
      signTens = sign (ONE, sigElp(1))
c *** compute von-mises equivalent stress
      qEl = abs(sigElp(1))
c *** compute current yield stress
      sigy    = sigy0 + dsigdep * pleq
c
      fratio = qEl / sigy - ONE
c *** check for yielding
      IF (sigy .LE. ZERO.or.fratio .LE. -SMALL) GO TO 500
c
      sigy_t   = sigy
c *** initial guess of incremental equivalent plastic strain   
      dpleq    = (qEl - sigy) / young
      pleq     = pleq_t + dpleq
      sigy     = sigy0 + dsigdep * pleq
c
c ***  update plastic strains, stresses
      epsPl(1) = epsPl(1) + dpleq * signTens
      stress(1) =  signTens * sigy
c
c ***  update plastic strains
      epseq  = pleq
c *** Update state variables
      ustatev(1) = pleq
      ustatev(2) = epsPl(1)
c *** Update plastic work
      sedPl = sedPl + HALF * (sigy_t + sigy) * dpleq
c
c *** Material Jcobian matrix
c
      dsdePl(1,1) = dsdeEl(1,1) * dsigdep /(dsdeEl(1,1) + dsigdep)
c *** Allow a small number for Jcobian matrix if it is an ideal plasticity
      if(dsdePl(1,1).LE.ZERO) dsdePl(1,1) = SMALL*dsdeEl(1,1)
c
      goto 600
  500 continue

c *** Update stress in case of elastic/unloading
      stress(1) = sigElp(1)

  600 continue
c *** elastic strain energy
      sedEl = HALF * stress(1) * (Strain(1)+dStrain(1)-epsPl(1))
c *** update state variables
      ustatev(nStatev) = sigy
c
      return
      end
*deck,usermat3d    USERDISTRIB  parallel                                gal
      subroutine usermat3d(
     &                   matId, elemId,kDomIntPt, kLayer, kSectPt,
     &                   ldstep,isubst,keycut,
     &                   nDirect,nShear,ncomp,nStatev,nProp,
     &                   Time,dTime,Temp,dTemp,
     &                   stress,ustatev,dsdePl,sedEl,sedPl,epseq,
     &                   Strain,dStrain, epsPl, prop, coords, 
     &                   var0, defGrad_t, defGrad,
     &                   tsstif, epsZZ,
     &                   var1, var2, var3, var4, var5,
     &                   var6, var7, var8)
c*************************************************************************
c     *** primary function ***
c
c           user defined material constitutive model
c
c      Attention:
c           User must define material constitutive law properly
c           according to the stress state such as 3D, plane strain
c           and axisymmetry, plane stress and beam.
c
c           a 3D material constitutive model can use for
c           plane strain and axisymmetry cases.
c
c           When using shell elements, a plane stress algorithm
c           must be use.
c
c                                             gal July, 1999
c
c       The following demonstrates a USERMAT subroutine for
c       a plasticity model of 3D solid elements or plane elements
c       in plane strain or axisymmetric stress state. The plasticity
c       model is the same as TB, BISO.
c       See "ANSYS user material subroutine USERMAT" for detailed
c       description of how to write a USERMAT routine.
c
c*************************************************************************
c
c     input arguments
c     ===============
c      matId     (int,sc,i)               material #
c      elemId    (int,sc,i)               element #
c      kDomIntPt (int,sc,i)               "k"th domain integration point
c      kLayer    (int,sc,i)               "k"th layer
c      kSectPt   (int,sc,i)               "k"th Section point
c      ldstep    (int,sc,i)               load step number
c      isubst    (int,sc,i)               substep number
c      nDirect   (int,sc,in)              # of direct components
c      nShear    (int,sc,in)              # of shear components
c      ncomp     (int,sc,in)              nDirect + nShear
c      nstatev   (int,sc,l)               Number of state variables
c      nProp     (int,sc,l)               Number of material ocnstants
c
c      Temp      (dp,sc,in)               temperature at beginning of
c                                         time increment
c      dTemp     (dp,sc,in)               temperature increment 
c      Time      (dp,sc,in)               time at beginning of increment (t)
c      dTime     (dp,sc,in)               current time increment (dt)
c
c      Strain   (dp,ar(ncomp),i)          Strain at beginning of time increment
c      dStrain  (dp,ar(ncomp),i)          Strain increment
c      prop     (dp,ar(nprop),i)          Material constants defined by TB,USER
c      coords   (dp,ar(3),i)              current coordinates
c      defGrad_t(dp,ar(3,3),i)            Deformation gradient at time t
c      defGrad  (dp,ar(3,3),i)            Deformation gradient at time t+dt
c
c     input output arguments              
c     ======================             
c      stress   (dp,ar(nTesn),io)         stress
c      ustatev   (dp,ar(nstatev),io)      user state variable
c            ustatev(1)                     - equivalent plastic strain
c            ustatev(2) - statev(1+ncomp)   - plastic strain vector
c            ustatev(nStatev)               - von-Mises stress
c      sedEl    (dp,sc,io)                elastic work
c      sedPl    (dp,sc,io)                plastic work
c      epseq    (dp,sc,io)                equivalent plastic strain
c      tsstif   (dp,ar(2),io)             transverse shear stiffness
c                                         tsstif(1) - Gxz
c                                         tsstif(2) - Gyz
c                                         tsstif(1) is also used to calculate hourglass
c                                         stiffness, this value must be defined when low
c                                         order element, such as 181, 182, 185 with uniform 
c                                         integration is used.
c      var?     (dp,sc,io)                not used, they are reserved arguments 
c                                         for further development
c
c     output arguments
c     ================
c      keycut   (int,sc,io)               loading bisect/cut control
c                                         0 - no bisect/cut
c                                         1 - bisect/cut 
c                                         (factor will be determined by ANSYS solution control)
c      dsdePl   (dp,ar(ncomp,ncomp),io)   material jacobian matrix
c      epsZZ    (dp,sc,o)                 strain epsZZ for plane stress,
c                                         define it when accounting for thickness change 
c                                         in shell and plane stress states
c
c*************************************************************************
c
c      ncomp   6   for 3D  (nshear=3)
c      ncomp   4   for plane strain or axisymmetric (nShear = 1)
c      ncomp   3   for plane stress (nShear = 1)
c      ncomp   3   for 3d beam      (nShear = 2)
c      ncomp   1   for 1D (nShear = 0)
c
c      stresss and strains, plastic strain vectors
c          11, 22, 33, 12, 23, 13    for 3D
c          11, 22, 33, 12            for plane strain or axisymmetry
c          11, 22, 12                for plane stress
c          11, 13, 12                for 3d beam
c          11                        for 1D
c
c      material jacobian matrix
c        3D
c           dsdePl    |  1111   1122   1133   1112   1123   1113 |
c           dsdePl    |  2211   2222   2233   2212   2223   2213 |
c           dsdePl    |  3311   3322   3333   3312   3323   3313 |
c           dsdePl    |  1211   1222   1233   1212   1223   1213 |
c           dsdePl    |  2311   2322   2333   2312   2323   2313 |
c           dsdePl    |  1311   1322   1333   1312   1323   1313 |
c        plane strain or axisymmetric (11, 22, 33, 12)
c           dsdePl    |  1111   1122   1133   1112 |
c           dsdePl    |  2211   2222   2233   2212 |
c           dsdePl    |  3311   3322   3333   3312 |
c           dsdePl    |  1211   1222   1233   1212 |
c        plane stress (11, 22, 12)
c           dsdePl    |  1111   1122   1112 |
c           dsdePl    |  2211   2222   2212 |
c           dsdePl    |  1211   1222   1212 |
c        3d beam (11, 13, 12)
c           dsdePl    |  1111   1113   1112 |
c           dsdePl    |  1311   1313   1312 |
c           dsdePl    |  1211   1213   1212 |
c        1d
c           dsdePl    |  1111 |
c
c*************************************************************************
#include "impcom.inc"
c
      INTEGER          
     &                 matId, elemId,
     &                 kDomIntPt, kLayer, kSectPt,
     &                 ldstep,isubst,keycut,
     &                 nDirect,nShear,ncomp,nStatev,nProp
      DOUBLE PRECISION 
     &                 Time,    dTime,   Temp,    dTemp,
     &                 sedEl,   sedPl,   epseq,   epsZZ
      DOUBLE PRECISION 
     &                 stress  (ncomp  ), ustatev (nStatev),
     &                 dsdePl  (ncomp,ncomp),
     &                 Strain  (ncomp  ), dStrain (ncomp  ), 
     &                 epsPl   (ncomp  ), prop    (nProp  ), 
     &                 coords  (3),
     &                 defGrad (3,3),     defGrad_t(3,3),
     &                 tsstif  (2)
c
c***************** User defined part *************************************
c
c --- parameters
c
      INTEGER          mcomp
      DOUBLE PRECISION HALF, THIRD, ONE, TWO, SMALL, ONEHALF,
     &                 ZERO, TWOTHIRD, ONEDM02, ONEDM05, sqTiny
      PARAMETER       (ZERO       = 0.d0,
     &                 HALF       = 0.5d0,
     &                 THIRD      = 1.d0/3.d0,
     &                 ONE        = 1.d0,
     &                 TWO        = 2.d0,
     &                 SMALL      = 1.d-08,
     &                 sqTiny     = 1.d-20,
     &                 ONEDM02    = 1.d-02,
     &                 ONEDM05    = 1.d-05,
     &                 ONEHALF    = 1.5d0,
     &                 TWOTHIRD   = 2.0d0/3.0d0,
     &                 mcomp      = 6
     &                 )
c
c --- local variables
c
c      sigElp   (dp,ar(6  ),l)            trial stress
c      dsdeEl   (dp,ar(6,6),l)            elastic moduli
c      sigDev   (dp,ar(6  ),l)            deviatoric stress tensor
c      dfds     (dp,ar(6  ),l)            derivative of the yield function 
c      JM       (dp,ar(6,6),l)            2D matrix for a 4 order tensor
c      pEl      (dp,sc     ,l)            hydrostatic pressure stress
c      qEl      (dp,sc     ,l)            von-mises stress
c      pleq_t   (dp,sc     ,l)            equivalent plastic strain at beginnig of time increment
c      pleq     (dp,sc     ,l)            equivalent plastic strain at end of time increment
c      dpleq    (dp,sc     ,l)            incremental equivalent plastic strain
c      sigy_t   (dp,sc     ,l)            yield stress at beginnig of time increments
c      sigy     (dp,sc     ,l)            yield stress at end of time increment
c      young    (dp,sc     ,l)            Young's modulus
c      posn     (dp,sc     ,l)            Poiss's ratio
c      sigy0    (dp,sc     ,l)            initial yield stress
c      dsigdep  (dp,sc     ,l)            plastic slop
c      twoG     (dp,sc     ,l)            two time of shear moduli
c      threeG   (dp,sc     ,l)            three time of shear moduli
c      elast1   (dp,sc     ,l)            lamé moduli
c      elast2   (dp,sc     ,l)            shear moduli
c      G                                  Kronecker delta G =[11 22 33 12 23 13]   
c --- temperary variables for solution purpose
c      i, j
c      threeOv2qEl, oneOv3G, qElOv3G, con1, con2, fratio
c
      EXTERNAL         vzero, vmove, get_ElmData
      DOUBLE PRECISION sigElp(mcomp), dsdeEl(mcomp,mcomp), G(mcomp),
     &                 sigDev(mcomp), JM    (mcomp,mcomp), dfds(mcomp),
     &                 sigi  (mcomp), strainEl(mcomp)

      DOUBLE PRECISION var0, var1, var2, var3, var4, var5,
     &                 var6, var7, var8,El(3,3),Ela(3,3),S(6)

      DATA G/1.0D0,1.0D0,1.0D0,0.0D0,0.0D0,0.0D0/
c
      INTEGER          i, j
      DOUBLE PRECISION pEl,   qEl,     pleq_t,  sigy_t , sigy,
     &                 dpleq, pleq, 
     &                 young, posn,    sigy0,   dsigdep, 
     &                 elast1,elast2,sM,bM,etrial(6),N(6),D(6,6),
     &                 twoG,  threeG,  oneOv3G, qElOv3G, threeOv2qEl, 
     &                 fratio,  con1,    con2, dperr(3),H,strial(6),
     &                 qtrial,phitrial
c*************************************************************************
c
      keycut   = 0
      dsigdep  = ZERO 
      pleq_t   = ustatev(1)
      pleq     = pleq_t
c *** get Young's modulus and Poisson's ratio, initial yield stress and others
      sM    = prop(1)
      bM     = prop(2)
      sigy0    = prop(3)
	  H        = prop(4)
	  
c *** plastic strain tensor
      call vmove(ustatev(2), epsPl(1), ncomp)
	  
c calculate stress
      call Fn12Etrial(defGrad,defGrad_t,etrial)
      call vec_print(etrial,6,"etrial")
	  etrial = Strain + dStrain - epsPl
      call vec_print(etrial,6,"etrial")
	  
      call 	hookeLaw(etrial,sM,bM,D,strial)  

c *** check for yielding
      call  j2calc(strial,N,qtrial)
      sigy    = sigy0 + H * pleq
      phitrial=qtrial-sigy
	  print *,"phitrial=",phitrial
      sigElp=strial
      if (phitrial .gt. 0.0d0) then
c *** plastic load
		dpleq=phitrial/(3.0d0*sM)
		pleq  = pleq_t + dpleq
        call plasticLaw(etrial,sM,bM,H,dpleq,qtrial,N,D,strial)
c ***  update stresses
        sigy    = sigy0 + H * pleq
        DO i = 1 , ncomp
            strial(i) =  sigElp(i) - TWOTHIRD * (qtrial-sigy) * N(i)
        END DO
      end if 	  
                 	  
	  
	  
! update Ansys Variables
      
c *** compute derivative of the yield function

c *** Update state variables
      ustatev(1) = pleq
      DO i = 1 , nDirect
         ustatev(i+1) = epsPl(i) + N(i) * dpleq
      END DO
      DO i = nDirect + 1 , ncomp
         ustatev(i+1) = epsPl(i) + TWO * N(i) * dpleq
      END DO
	  stress=strial

	  dsdePl=D
	  ustatev(nStatev) = sigy
	  call tensor4ord_print(dsdePl,"Jacobian Matrix")
      tsstif(:)=sM
	  
	  sedEl = ZERO
      DO i = 1 , ncomp
         sedEl = sedEl + stress(i)*(Strain(i)+dStrain(i)-epsPl(i))
      END DO
      sedEl    = sedEl * HALF
      ustatev(nStatev) = sigy
      return
      end
*deck,usermatbm    USERDISTRIB  parallel                                gal
      subroutine usermatbm(
     &                   matId, elemId,kDomIntPt, kLayer, kSectPt,
     &                   ldstep,isubst,keycut,
     &                   nDirect,nShear,ncomp,nStatev,nProp,
     &                   Time,dTime,Temp,dTemp,
     &                   stress,ustatev,dsdePl,sedEl,sedPl,epseq,
     &                   Strain,dStrain, epsPl, prop, coords, 
     &                   var0, defGrad_t, defGrad,
     &                   tsstif, epsZZ,
     &                   var1, var2, var3, var4, var5,
     &                   var6, var7, var8)
c*************************************************************************
c     *** primary function ***
c
c           user defined material constitutive model
c
c      Attention:
c           User must define material constitutive law properly
c           according to the stress state such as 3D, plane strain
c           and axisymmetry, plane stress and beam.
c
c           a 3D material constitutive model can use for
c           plane strain and axisymmetry cases.
c
c           When using shell elements, a plane stress algorithm
c           must be use.
c
c                                             gal July, 1999
c
c       The following demonstrates a USERMAT subroutine for 
c       a plasticity model in 3D beam(188, 189). The plasticity
c       model is the same as TB, BISO.
c       See "ANSYS user material subroutine USERMAT" for detailed
c       description of how to write a USERMAT routine.
c
c*************************************************************************
c
c     input arguments
c     ===============
c      matId     (int,sc,i)               material #
c      elemId    (int,sc,i)               element #
c      kDomIntPt (int,sc,i)               "k"th domain integration point
c      kLayer    (int,sc,i)               "k"th layer
c      kSectPt   (int,sc,i)               "k"th Section point
c      ldstep    (int,sc,i)               load step number
c      isubst    (int,sc,i)               substep number
c      nDirect   (int,sc,in)              # of direct components
c      nShear    (int,sc,in)              # of shear components
c      ncomp     (int,sc,in)              nDirect + nShear
c      nStatev   (int,sc,l)               Number of state variables
c      nProp     (int,sc,l)               Number of material ocnstants
c
c      Temp      (dp,sc,in)               temperature at beginning of
c                                         time increment
c      dTemp     (dp,sc,in)               temperature increment 
c      Time      (dp,sc,in)               time at beginning of increment (t)
c      dTime     (dp,sc,in)               current time increment (dt)
c
c      Strain   (dp,ar(ncomp),i)          Strain at beginning of time increment
c      dStrain  (dp,ar(ncomp),i)          Strain increment
c      prop     (dp,ar(nprop),i)          Material constants defined by TB,USER
c      coords   (dp,ar(3),i)              current coordinates
c      defGrad_t(dp,ar(3,3),i)            Deformation gradient at time t
c      defGrad  (dp,ar(3,3),i)            Deformation gradient at time t+dt
c
c     input output arguments              
c     ======================             
c      stress   (dp,ar(nTesn),io)         stress
c      ustatev   (dp,ar(nStatev),io)       statev
c           ustatev(1)                     - equivalent plastic strain
c           ustatev(2) - ustatev(1+ncomp)  - plastic strain vector
c           ustatev(nStatev)               - von-Mises stress
c      sedEl    (dp,sc,io)                elastic work
c      sedPl    (dp,sc,io)                plastic work
c      epseq    (dp,sc,io)                equivalent plastic strain
c      tsstif   (dp,ar(2),io)             transverse shear stiffness
c                                         tsstif(1) - Gxz
c                                         tsstif(2) - Gyz
c                                         tsstif(1) is also used to calculate hourglass
c                                         stiffness, this value must be defined when low
c                                         order element, such as 181, 182, 185 with uniform 
c                                         integration is used.
c      var?     (dp,sc,io)                not used, they are reserved arguments 
c                                         for further development
c
c     output arguments
c     ================
c      keycut   (int,sc,io)               loading bisect/cut control
c                                         0 - no bisect/cut
c                                         1 - bisect/cut 
c                                         (factor will be determined by ANSYS solution control)
c      dsdePl   (dp,ar(ncomp,ncomp),io)   material jacobian matrix
c      epsZZ    (dp,sc,o)                 strain epsZZ for plane stress,
c                                         define it when accounting for thickness change 
c                                         in shell and plane stress states
c
c*************************************************************************
c
c      ncomp   6   for 3D
c      ncomp   4   for plane strain, axisymmetric (nShear = 1)
c      ncomp   3   for plane stress (nShear = 1)
c      ncomp   3   for 3D beam (nShear = 2), beam188/189
c      ncomp   1   for 1D beam, link180
c
c      stresss and strains, plastic strain vectors
c          11, 22, 33, 12, 23, 13    for 3D
c          11, 22, 33, 12            for Plane strain and axisymmetry
c          11, 22, 12                for Plane stress
c          11, 13, 12                for 3d beam
c          11                        for 1D
c
c      material jacobian matrix
c        3D
c           dsdePl    |  1111   1122   1133   1112   1123   1113 |
c           dsdePl    |  2211   2222   2233   2212   2223   2213 |
c           dsdePl    |  3311   3322   3333   3312   3323   3313 |
c           dsdePl    |  1211   1222   1233   1212   1223   1213 |
c           dsdePl    |  2311   2322   2333   2312   2323   2313 |
c           dsdePl    |  1311   1322   1333   1312   1323   1313 |
c        plane strain, axisymmetric
c           dsdePl    |  1111   1122   1133   1112 |
c           dsdePl    |  2211   2222   2233   2212 |
c           dsdePl    |  3311   3322   3333   3312 |
c           dsdePl    |  1211   1222   1233   1212 |
c        plane stress
c           dsdePl    |  1111   1122   1112 |
c           dsdePl    |  2211   2222   2212 |
c           dsdePl    |  1211   1222   1212 |
c        3d beam plasticity
c           dsdePl    |  1111   1113   1112 |
c           dsdePl    |  1311   1313   1312 |
c           dsdePl    |  1211   1213   1212 |
c        1d
c           dsdePl    |  1111 |
c
c*************************************************************************
#include "impcom.inc"
c
      INTEGER          
     &                 matId, elemId,
     &                 kDomIntPt, kLayer, kSectPt,
     &                 ldstep,isubst,keycut,
     &                 nDirect,nShear,ncomp,nStatev,nProp
      DOUBLE PRECISION 
     &                 Time,    dTime,   Temp,    dTemp,
     &                 sedEl,   sedPl,   epseq,   epsZZ
      DOUBLE PRECISION 
     &                 stress  (ncomp  ), ustatev (nStatev),
     &                 dsdePl  (ncomp,ncomp), sigi(ncomp),
     &                 Strain  (ncomp  ), dStrain (ncomp  ), 
     &                 epsPl   (ncomp  ), prop    (nProp  ), 
     &                 coords  (3),       
     &                 defGrad (3,3),     defGrad_t(3,3),
     &                 tsstif  (2)
c
c***************** User defined part *************************************
c
c --- parameters
c
      INTEGER          NEWTON, mcomp
      DOUBLE PRECISION HALF, ONE, TWO, SMALL, SQTWOTHIRD,
     &                 ZERO, TWOTHIRD, ONEDM02, ONEDM05, sqTiny
      PARAMETER       (ZERO       = 0.d0,
     &                 HALF       = 0.5d0,
     &                 ONE        = 1.d0,
     &                 TWO        = 2.d0,
     &                 SMALL      = 1.d-08,
     &                 sqTiny     = 1.d-20,
     &                 ONEDM02    = 1.d-02,
     &                 ONEDM05    = 1.d-05,
     &                 TWOTHIRD   = 2.0d0/3.0d0,
     &                 SQTWOTHIRD = 0.816496580927726030d0,
     &                 NEWTON     = 20,
     &                 mcomp      = 3
     &                 )
c
c --- local variables
c
c      sigElp   (dp,ar(3  ),l)            trial stress
c      dsdeEl   (dp,ar(3,3),l)            elastic moduli
c      pleq_t   (dp,sc     ,l)            equivalent plastic strain at beginnig of time increment
c      pleq     (dp,sc     ,l)            equivalent plastic strain at end of time increment
c      dpleq    (dp,sc     ,l)            incremental equivalent plastic strain
c      gamma    (dp,sc     ,l)            variable for solving incremental equivalent plastic strain
c      dgamma   (dp,sc     ,l)            correction of gamma
c      sigy_t   (dp,sc     ,l)            yield stress at beginnig of time increments
c      sigy     (dp,sc     ,l)            yield stress at end of time increment
c      young    (dp,sc     ,l)            Young's modulus
c      posn     (dp,sc     ,l)            Poiss's ratio
c      sigy0    (dp,sc     ,l)            initial yield stress
c      dsigdep  (dp,sc     ,l)            plastic slop
c      twoG     (dp,sc     ,l)            two time of shear moduli
c      funcf    (dp,sc     ,l)            nonlinear function to be solved for gamma
c      dFdep    (dp,sc     ,l)            derivative of nonlinear function over gamma
c
c --- temperary variables for solution purpose
c      i, j
c      c1, c2, c3, fratio
c      wk1(3), wk2(3), wk3(3), wk4(3) vector working arrays
c
      EXTERNAL         vmove, vzero, vapb1, vamb1,get_ElmData
      DOUBLE PRECISION sigElp(mcomp), dsdeEl(mcomp,mcomp), 
     &                 wk1(3), wk2(3), wk3(3), wk4(3)

      DOUBLE PRECISION var0, var1, var2, var3, var4, var5,
     &                 var6, var7, var8

      INTEGER          i, j, k
      DOUBLE PRECISION pleq_t,  sigy_t , sigy,
     &                 cpleq, dpleq,   pleq,    twoG,    et,
     &                 young, posn,    sigy0,   dsigdep, 
     &                 gamma, dgamma,  dfdga,   dplga,   fratio,
     &                 funcFb,funcFb2, funcf,   dFdep,
     &                 c1, c2, c3, c4, c5
      DOUBLE PRECISION pv(3)
      data pv/TWOTHIRD, TWO, TWO/
c*************************************************************************
c
      keycut   = 0
c *** equivalent plastic strain at beginning of time step
      pleq_t   = ustatev(ncomp+1)
      pleq     = pleq_t
c *** get Young's modulus and Poisson's ratio, initial yield stress and slope of stress-strain
      young    = prop(1)
      posn     = prop(2)
      sigy0    = prop(3)
      et       = prop(4)
c *** calculate plastic slope
      dsigdep  = young * et/(young - et)
      twoG     = young / (ONE+posn)
c *** define tsstif(1) since it is used for calculation of hourglass stiffness
      tsstif(1) = HALF * twoG
c
c *** calculate elastic stiffness matrix
c
      call vzero(dsdeEl(1,1), ncomp * ncomp)
      c1 = twoG * HALF
      dsdeEl (1,1) = young
      dsdeEl (2,2) = c1
      dsdeEl (3,3) = c1
      DO i = 1, ncomp
         wk3(i) = dsdeEl(i,i)
      END DO
c *** calculate predicted strain 
      call vmove(Strain(1), wk1(1), ncomp)
      call vapb1(wk1(1), dStrain(1), ncomp)
      call vamb1(wk1(1), ustatev(1), ncomp)

c
c *** get initial stress
      call vzero(sigi(1),ncomp)
      i = ncomp
      call get_ElmData ('ISIG', elemId,kDomIntPt, i, sigi)

c
c *** calculate the trial stress and 
c     copy elastic moduli dsdeEl to material Jacobian matrix
      call vmove(dsdeEl(1,1), dsdePl(1,1), ncomp * ncomp)
      do i=1,ncomp
         sigElp(i) = wk3(i) * wk1(i) + sigi(i)
      end do
c
      funcFb2 = ZERO
      DO i = 1, ncomp
        funcFb2 = funcFb2 + pv(i) * sigElp(i) * sigElp(i)
      END DO
      funcFb = sqrt(funcFb2)

c *** compute current yield stress
      sigy    = sigy0 + dsigdep * pleq
c
      fratio = funcFb/sigy - SQTWOTHIRD
c *** check for yielding
      IF (fratio .LE. -SMALL) GO TO 500
      sigy_t  = sigy

      DO i = 1, ncomp
         wk3(i) = wk3(i) * pv(i)
      END DO

      gamma    = ZERO
      dplga    = ZERO
      dfdga    = ZERO
      dpleq    = ZERO
      pleq     = pleq_t 

c *** Local New-Raphson procedure for solving the gamma and 
c     thus the incremental equivalent plastic strain
      DO k=1,NEWTON
         funcFb2 = ZERO
         dfdga   = ZERO
         DO j = 1 , ncomp
            c1 = ONE + gamma * wk3(j)
            c1 = ONE / c1
            c2 = sigElp(j) * c1
            wk4(j) = c2
            funcFb2 = funcFb2 + pv(j) * c2 * c2
            c2 = c2 * c2 * c1 * wk3(j) * pv(j)
            dfdga   = dfdga - c2
         END DO
         funcFb   = sqrt(funcFb2)
c ***    derivative of funcFb w.r.t. gamma
         dfdga   = dfdga / funcFb

c ***    calculate the incremental equivalent plastic strain
         dpleq    = gamma * SQTWOTHIRD * funcFb
c ***    update the total equivalent plastic strain
         pleq     = pleq_t + dpleq
c ***    current yield stress
         sigy     = sigy0 + dsigdep * pleq
c ***    calculate the residual
         funcf    = funcFb - SQTWOTHIRD * sigy
c ***    derivative of incremental equivalent plastic strain w.r.t. gamma
         dplga    = SQTWOTHIRD * (gamma * dfdga + funcFb)
c ***    derivative of residual function w.r.t. gamma
         dFdep    = dfdga - SQTWOTHIRD * dsigdep * dplga
c ***    correction of gamma
         dgamma   = -funcf / dFdep 
         gamma    = gamma   + dgamma
c ***    check for negative gamma
         gamma    = max (gamma, sqTiny)
         fratio   = funcf/ sigy
c
c ***    Check for convergence of local New-Raphson iteration
         IF (((abs(fratio) .LT. ONEDM05 ) .AND.
     &        (abs(dgamma) .LT. ONEDM02*gamma)) .OR.
     &       ((abs(fratio) .LT. ONEDM05 ) .AND.
     &        (dgamma      .LE. sqTiny  ) .AND.
     &        ( gamma      .LE. sqTiny  )))  GO TO 100

      END DO
c
c *** Uncovergence, set keycut to 1 for bisection/cutback the load increment
      keycut   = 1
      GO TO 990
 100  CONTINUE
c
c *** update stresses
      call vmove(wk4(1), stress(1), ncomp)

c *** calculate incremental plastic strain 
      DO j = 1, ncomp
         wk2(j) = gamma * pv(j) * wk4(j)
      END DO
c *** update plastic strains
      call vapb1(epsPl(1),wk2(1),ncomp)

c *** Update state variables
      ustatev(ncomp+1) = pleq
      do i=1,ncomp
         ustatev(i) = epsPl(i)
      end do

c *** update plastic work
      sedPl     = sedPl + HALF * (sigy_t+sigy) * dpleq

      c1     = TWOTHIRD * dsigdep
      c3     = c1 * funcFb2 / (ONE - c1 * gamma)
      DO j = 1 , ncomp
         c1 = ONE / (ONE + gamma * wk3(j))
         wk3(j) = wk3(j) * c1 / pv(j)
      END DO
      DO j = 1 , ncomp
         wk4(j) = wk4(j) * pv(j)
      END DO
      DO j = 1 , ncomp
         c3 = c3 + wk4(j) * wk4(j) * wk3(j)
      END DO
      DO j = 1 , ncomp
         wk4(j) = wk4(j) * wk3(j)
      END DO

      c3     = ONE / c3
      DO i=1,ncomp
         dsdePl(i,i) = wk3(i)
      END DO

c *** Calculate the plastic Jacobian

      DO i=1,ncomp
         DO j=1,ncomp
            dsdePl(i,j) =    dsdePl(i,j) - c3 * wk4(i) * wk4(j)
         END DO
      END DO

      goto 600

  500 continue

c *** Update stress in case of elastic/unloading
      call vmove(sigElp(1),stress(1),ncomp)

  600 continue
c *** elastic strain energy
      sedEl = ZERO
      DO i = 1 , ncomp
         sedEl = sedEl + stress(i)*(Strain(i)+dStrain(i)-epsPl(i))
      END DO
      sedEl = sedEl * HALF
      ustatev(nStatev) = funcFb / SQTWOTHIRD

 990  CONTINUE
c
      return
      end
*deck,usermatps    USERDISTRIB  parallel                                gal
      subroutine usermatps(
     &                   matId, elemId,kDomIntPt, kLayer, kSectPt,
     &                   ldstep,isubst,keycut,
     &                   nDirect,nShear,ncomp,nStatev,nProp,
     &                   Time,dTime,Temp,dTemp,
     &                   stress,ustatev,dsdePl,sedEl,sedPl,epseq,
     &                   Strain,dStrain, epsPl, prop, coords, 
     &                   var0, defGrad_t, defGrad,
     &                   tsstif, epsZZ, 
     &                   var1, var2, var3, var4, var5, 
     &                   var6, var7, var8)
c*************************************************************************
c     *** primary function ***
c
c           user defined material constitutive model
c
c      Attention:
c           User must define material constitutive law properly
c           according to the stress state such as 3D, plane strain
c           and axisymmetry, plane stress and beam.
c
c           a 3D material constitutive model can use for
c           plane strain and axisymmetry cases.
c
c           When using shell elements, a plane stress algorithm
c           must be use.
c
c                                             gal July, 1999
c
c       The following demonstrates a USERMAT subroutine for
c       a plasticity model of plane stress state (such as PLANE182,
c       PLANE183 or SHELL181). The plasticity model is the same 
c       as ANSYS TB, BISO.
c
c       See "ANSYS user material subroutine USERMAT" for detailed
c       description of how to write a USERMAT routine.
c
c*************************************************************************
c
c     input arguments
c     ===============
c      matId     (int,sc,i)               material #
c      elemId    (int,sc,i)               element #
c      kDomIntPt (int,sc,i)               "k"th domain integration point
c      kLayer    (int,sc,i)               "k"th layer
c      kSectPt   (int,sc,i)               "k"th Section point
c      ldstep    (int,sc,i)               load step number
c      isubst    (int,sc,i)               substep number
c      nDirect   (int,sc,in)              # of direct components
c      nShear    (int,sc,in)              # of shear components
c      ncomp     (int,sc,in)              nDirect + nShear
c      nStatev   (int,sc,l)               Number of state variables
c      nProp     (int,sc,l)               Number of material ocnstants
c
c      Temp      (dp,sc,in)               temperature at beginning of
c                                         time increment
c      dTemp     (dp,sc,in)               temperature increment 
c      Time      (dp,sc,in)               time at beginning of increment (t)
c      dTime     (dp,sc,in)               current time increment (dt)
c
c      Strain   (dp,ar(ncomp),i)          Strain at beginning of time increment
c      dStrain  (dp,ar(ncomp),i)          Strain increment
c      prop     (dp,ar(nprop),i)          Material constants defined by TB,USER
c      coords   (dp,ar(3),i)              current coordinates
c      defGrad_t(dp,ar(3,3),i)            Deformation gradient at time t
c      defGrad  (dp,ar(3,3),i)            Deformation gradient at time t+dt
c
c     input output arguments              
c     ======================             
c      stress   (dp,ar(nTesn),io)         stress
c      ustatev   (dp,ar(nStatev),io)      user state variables
c            ustatev(1)                     - equivalent plastic strain
c            ustatev(2) - ustatev(1+ncomp)   - plastic strain vector
c            ustatev(nStatev)               - von-Mises stress
c      sedEl    (dp,sc,io)                elastic work
c      sedPl    (dp,sc,io)                plastic work
c      epseq    (dp,sc,io)                equivalent plastic strain
c      tsstif   (dp,ar(2),io)             transverse shear stiffness
c                                         tsstif(1) - Gxz
c                                         tsstif(2) - Gyz
c                                         tsstif(1) is also used to calculate hourglass
c                                         stiffness, this value must be defined when low
c                                         order element, such as 181, 182, 185 with uniform 
c                                         integration is used.
c      var?     (dp,sc,io)                not used, they are reserved arguments 
c                                         for further development
c
c     output arguments
c     ================
c      keycut   (int,sc,io)               loading bisect/cut control
c                                         0 - no bisect/cut
c                                         1 - bisect/cut 
c                                         (factor will be determined by ANSYS solution control)
c      dsdePl   (dp,ar(ncomp,ncomp),io)   material jacobian matrix
c      epsZZ    (dp,sc,o)                 strain epsZZ for plane stress, 
c                                         define it when accounting for thickness change 
c                                         in shell and plane stress states
c
c*************************************************************************
c
c      ncomp   6   for 3D  (nshear=3)
c      ncomp   4   for plane strain or axisymmetric (nShear = 1)
c      ncomp   3   for plane stress (nShear = 1)
c      ncomp   3   for 3d beam      (nShear = 2)
c      ncomp   1   for 1D (nShear = 0)
c
c      stresss and strains, plastic strain vectors
c          11, 22, 33, 12, 23, 13    for 3D
c          11, 22, 33, 12            for plane strain or axisymmetry
c          11, 22, 12                for plane stress
c          11, 13, 12                for 3d beam
c          11                        for 1D
c
c      material jacobian matrix
c        3D
c           dsdePl    |  1111   1122   1133   1112   1123   1113 |
c           dsdePl    |  2211   2222   2233   2212   2223   2213 |
c           dsdePl    |  3311   3322   3333   3312   3323   3313 |
c           dsdePl    |  1211   1222   1233   1212   1223   1213 |
c           dsdePl    |  2311   2322   2333   2312   2323   2313 |
c           dsdePl    |  1311   1322   1333   1312   1323   1313 |
c        plane strain or axisymmetric (11, 22, 33, 12)
c           dsdePl    |  1111   1122   1133   1112 |
c           dsdePl    |  2211   2222   2233   2212 |
c           dsdePl    |  3311   3322   3333   3312 |
c           dsdePl    |  1211   1222   1233   1212 |
c        plane stress (11, 22, 12)
c           dsdePl    |  1111   1122   1112 |
c           dsdePl    |  2211   2222   2212 |
c           dsdePl    |  1211   1222   1212 |
c        3d beam (11, 13, 12)
c           dsdePl    |  1111   1113   1112 |
c           dsdePl    |  1311   1313   1312 |
c           dsdePl    |  1211   1213   1212 |
c        1d
c           dsdePl    |  1111 |
c
c*************************************************************************
#include "impcom.inc"
c
      INTEGER          
     &                 matId, elemId,
     &                 kDomIntPt, kLayer, kSectPt,
     &                 ldstep,isubst,keycut,
     &                 nDirect,nShear,ncomp,nStatev,nProp
      DOUBLE PRECISION 
     &                 Time,    dTime,   Temp,    dTemp,
     &                 sedEl,   sedPl,   epseq,   epsZZ
      DOUBLE PRECISION 
     &                 stress  (ncomp  ), ustatev (nStatev),
     &                 dsdePl  (ncomp,ncomp), sigi(ncomp),
     &                 Strain  (ncomp  ), dStrain (ncomp  ), 
     &                 epsPl   (ncomp  ), prop    (nProp  ), 
     &                 coords  (3),       
     &                 defGrad (3,3),     defGrad_t(3,3),
     &                 tsstif  (2)
c
c***************** User defined part *************************************
c
c --- parameters
c
      INTEGER          NEWTON, mcomp
      DOUBLE PRECISION HALF, THIRD, ONE, TWO, SMALL, 
     &                 SQTWOTHIRD, SQTWO1,
     &                 ZERO, TWOTHIRD, ONEDM02, ONEDM05, sqTiny
      PARAMETER       (ZERO       = 0.d0,
     &                 HALF       = 0.5d0,
     &                 THIRD      = 1.d0/3.d0,
     &                 ONE        = 1.d0,
     &                 TWO        = 2.d0,
     &                 SMALL      = 1.d-08,
     &                 sqTiny     = 1.d-20,
     &                 ONEDM02    = 1.d-02,
     &                 ONEDM05    = 1.d-05,
     &                 TWOTHIRD   = 2.0d0/3.0d0,
     &                 SQTWOTHIRD = 0.816496580927726030d0,
     &                 SQTWO1     = 0.707106769084930420d0,
     &                 NEWTON     = 20,
     &                 mcomp      = 6
     &                 )
c
c --- local variables
c
c      sigElp   (dp,ar(6  ),l)            trial stress
c      dsdeEl   (dp,ar(6,6),l)            elastic moduli
c      pleq_t   (dp,sc     ,l)            equivalent plastic strain at beginnig of time increment
c      pleq     (dp,sc     ,l)            equivalent plastic strain at end of time increment
c      dpleq    (dp,sc     ,l)            incremental equivalent plastic strain
c      gamma    (dp,sc     ,l)            variable for solving incremental equivalent plastic strain
c      dgamma   (dp,sc     ,l)            correction of gamma
c      sigy_t   (dp,sc     ,l)            yield stress at beginnig of time increments
c      sigy     (dp,sc     ,l)            yield stress at end of time increment
c      young    (dp,sc     ,l)            Young's modulus
c      posn     (dp,sc     ,l)            Poiss's ratio
c      sigy0    (dp,sc     ,l)            initial yield stress
c      dsigdep  (dp,sc     ,l)            plastic slop
c      twoG     (dp,sc     ,l)            two time of shear moduli
c      funcf    (dp,sc     ,l)            nonlinear function to be solved for dpleq
c      dFdep    (dp,sc     ,l)            derivative of nonlinear function over dpleq
c
c --- temperary variables for solution purpose
c      i, j
c      con1eOv2qEl, oneOv3G, qElOv3G, con1, con2, fratio
c
      EXTERNAL         vmove, vzero, vapb1, get_ElmData
      DOUBLE PRECISION sigElp(mcomp), dsdeEl(mcomp,mcomp), 
     &                 wk1(3), wk2(3), wk3(3), wk4(3)

      DOUBLE PRECISION var0, var1, var2, var3, var4, var5,
     &                 var6, var7, var8

      INTEGER          i, j, k
      DOUBLE PRECISION pleq_t,  sigy_t , sigy,
     &                 dpleq,   pleq,    twoG,    et,
     &                 young, posn,    sigy0,   dsigdep, tEo1pm,
     &                 gamma, dgamma,  dfdga,   dplga, 
     &                 funcFb,funcFb2, funcf,   dFdep,   fratio,
     &                 con1,  con2,    con3,  con4,
     &                 con2p1, ocon2p1,
     &                 ocon2p2, con4p1, ocon4p1, ocon4p2,
     &                 c1, c2, c3,c4, c5,dperr(3)
c*************************************************************************
c
      keycut   = 0
c *** equivalent plastic strain at beginning of time step
      pleq_t   = ustatev(1)
      pleq     = pleq_t
c *** get Young's modulus and Poisson's ratio, initial yield stress and slope of stress-strain
      young    = prop(1)
      posn     = prop(2)
      sigy0    = prop(3)
      et       = prop(4)
c *** calculate plastic slope
      dsigdep  = young * et/(young - et)
      twoG     = young / (ONE+posn)
      tEo1pm   = THIRD * young /(ONE - posn)
c *** define tsstif(1) since it is used for calculation of hourglass stiffness
      tsstif(1) = HALF * twoG
c
c *** calculate elastic stiffness matrix (3d)
c
      c1 = ONE - posn * posn
      c2 = young / c1
      c3 = posn * c2
      dsdeEl (1,1) = c2
      dsdeEl (1,2) = c3
      dsdeEl (1,3) = ZERO
      dsdeEl (2,2) = c2
      dsdeEl (2,3) = ZERO
      dsdeEl (3,3) = HALF * twoG
      do i=1,ncomp-1
        do j=i+1,ncomp
          dsdeEl(j,i)=dsdeEl(i,j)
        end do
      end do
      call vmove(ustatev(2), epsPl(1), ncomp)
c *** calculate elastic strain
      do i=1,ncomp
         wk1(i) = Strain(i) - epsPl(i) + dStrain(i)
      end do

c
c *** get initial stress
      call vzero(sigi(1),ncomp)
      i = ncomp
      call get_ElmData ('ISIG', elemId,kDomIntPt, i, sigi)

c
c *** calculate the trial stress and 
c     copy elastic moduli dsdeEl to material Jacobian matrix
      do i=1,ncomp
         sigElp(i) = ZERO
         do j=1,ncomp
            dsdePl(j,i) = dsdeEl(j,i)
            sigElp(i) = sigElp(i) + dsdeEl(j,i) * wk1(j)
         end do
         sigElp(i) = sigElp(i) + sigi(i)
      end do
c
      wk1(1)   = SQTWO1 * ( sigElp(1) +  sigElp(2) )
      wk1(2)   = SQTWO1 * (-sigElp(1) +  sigElp(2) )
      wk1(3)   = sigElp(3)
c
      funcFb2 =  THIRD * wk1(1)  * wk1(1) +
     &   wk1(2)  * wk1(2) + TWO *  wk1(3)  * wk1(3)

      funcFb  = sqrt(funcFb2)

c *** compute current yield stress
      sigy    = sigy0 + dsigdep * pleq
c
      fratio = funcFb/sigy - SQTWOTHIRD
c *** check for yielding
      IF (fratio .LE. -SMALL) GO TO 500
      sigy_t  = sigy

      gamma    = ZERO
      dplga    = ZERO
      dfdga    = ZERO
      con1     = THIRD * wk1(1) * wk1(1)
      con2     = wk1(2) * wk1(2) + TWO * wk1(3) * wk1(3)
      con3     = - con1 * tEo1pm
      con4     = - con2 * twoG
      funcf    = funcFb - SQTWOTHIRD * sigy
      dfdga    =  (con3 + con4) / funcFb
      dFdep    = dfdga - TWOTHIRD * dsigdep * funcFb
      dgamma   = -funcf / dFdep
      gamma    = gamma   + dgamma
      gamma    = max (gamma, sqTiny)

      DO k=1,NEWTON

         con2p1   = ONE + tEo1pm * gamma
         con4p1   = ONE + twoG * gamma
         ocon2p1  = ONE / con2p1
         ocon2p2  = ocon2p1 * ocon2p1
         ocon4p1  = ONE / con4p1
         ocon4p2  = ocon4p1 * ocon4p1
         funcFb2  =  con1 * ocon2p2 + con2 * ocon4p2
         funcFb   = sqrt(funcFb2)
c
c ***    calculate the incremental equivalent plastic strain
         dpleq    = gamma * SQTWOTHIRD * funcFb
c ***    update the total equivalent plastic strain
         pleq     = pleq_t + dpleq
c ***    current yield stress
         sigy     = sigy0 + dsigdep * pleq

c ***    calculate the residual
         funcf    = funcFb - SQTWOTHIRD * sigy
c ***    derivative of funcFb w.r.t. gamma
         dfdga    =  ( con3 * ocon2p2 * ocon2p1
     &               + con4 * ocon4p2 * ocon4p1) / funcFb
c ***    derivative of incremental equivalent plastic strain w.r.t. gamma
         dplga    = SQTWOTHIRD * (gamma * dfdga + funcFb)
c ***    derivative of residual function w.r.t. gamma
         dFdep  = dfdga - SQTWOTHIRD * dsigdep * dplga
         fratio   = funcf/ sigy

         dgamma   = -funcf / dFdep
         gamma    = gamma   + dgamma
         gamma    = max (gamma, sqTiny)
c
         IF (((abs(fratio) .LT. ONEDM05 ) .AND.
     &        (abs(dgamma) .LT. ONEDM02*gamma)) .OR.
     &       ((abs(fratio) .LT. ONEDM05 ) .AND.
     &        (dgamma      .LE. sqTiny  ) .AND.
     &        ( gamma      .LE. sqTiny  )))  GO TO 100

      END DO
c
c *** Uncovergence, set keycut to 1 for bisect/cut
      keycut   = 1
      GO TO 990
 100  CONTINUE
c
c *** update stresses
      c1 = ONE / ( ONE + gamma * tEo1pm)
      c2 = ONE / ( ONE + gamma * twoG  )
      stress(1) = SQTWO1 * (c1 * wk1(1) - c2 * wk1(2))
      stress(2) = SQTWO1 * (c1 * wk1(1) + c2 * wk1(2))
      stress(3) = c2 * wk1(3)
c
      wk1(1)    = SQTWO1 * ( stress(1) +  stress(2))
      wk1(2)    = SQTWO1 * (-stress(1) +  stress(2))
      wk1(3)    = stress(3)
      con1  = SQTWO1 * gamma
c *** calculate incremental plastic strain, wk2(i)
      wk2(1)    = con1 * ( THIRD * wk1(1) -  wk1(2))
      wk2(2)    = con1 * ( THIRD * wk1(1) +  wk1(2))
      wk2(3)    = TWO * gamma * wk1 (3)

c ***  update plastic strains
      call vapb1(epsPl(1),wk2(1),ncomp)
      epseq     = pleq
      ustatev(1) = pleq
c *** calculate plastic work
      sedPl     = sedPl + HALF * (sigy_t+sigy) * dpleq

c *** consistent tangent stiffness matrix
      con3      = TWOTHIRD*dsigdep*funcFb2/(ONE-TWOTHIRD*dsigdep*gamma)
      call vmove(wk1(1),wk4(1),3)
      c1 = wk4(1)
      c2 = wk4(2)
      c3 = wk4(3)
      c4 = tEo1pm / (ONE + tEo1pm * gamma)
      c5 = twoG / (ONE + twoG * gamma)
      wk4(1) = SQTWO1 * ( c4 * c1 - c5 * c2)
      wk4(2) = SQTWO1 * ( c4 * c1 + c5 * c2)
      wk4(3) = c5 * c3

      call vmove(wk4(1),wk3(1),3)
      c1 = wk3(1)
      c2 = wk3(2)
      wk3(1) = SQTWO1 * ( c1+c2)
      wk3(2) = SQTWO1 * (-c1+c2)

      con3     = con3 + THIRD*wk1(1)*wk3(1)+
     &           wk1(2)*wk3(2)+TWO*wk1(3)*wk3(3)
      c1  = 3.0d0 * tEo1pm / (ONE + tEo1pm * gamma)
      c2  = twoG / (ONE + twoG * gamma)
      call vzero(dsdePl(1,1),9)
      dsdePl(1,1) = HALF * (c1 + c2)
      dsdePl(2,1) = HALF * (c1 - c2)
      dsdePl(1,2) = dsdePl(2,1)
      dsdePl(2,2) = dsdePl(1,1)
      dsdePl(3,3) = HALF * c2
      con3     = ONE / con3
      DO i=1,3
         DO j=1,3
            dsdePl(i,j) = dsdePl(i,j) - con3 * wk4(i) * wk4(j)
         END DO
      END DO

      goto 600

  500 continue

c *** Update stress in case of elastic/unloading
      call vmove(sigElp(1),stress(1),ncomp)

  600 continue
c *** elastic strain and put into epsZZ
      epsZZ    = -posn/young * (stress(1) + stress(2))
c *** add plastic strain to total strain epsZZ
      epsZZ    = epsZZ - (epsPl(1) + epsPl(2))

c *** elastic strain energy
      sedEl = ZERO
      DO i = 1 , ncomp
         sedEl = sedEl + stress(i)*(Strain(i)+dStrain(i)-epsPl(i))
      END DO
      sedEl    = sedEl * HALF
      ustatev(nStatev) = sigy
      do i=1,ncomp
         ustatev(i+1) = epsPl(i)
      end do

 990  CONTINUE
c
      return
      end
