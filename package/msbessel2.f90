!******************************************************************
!*       Purpose: This program computes the modified spherical    *
!*                Bessel functions kn(x) and kn'(x) using         *
!*                subroutine SPHK                                 *
!*       Input :  x --- Argument of kn(x)  ( x > 0 )              *
!*                n --- Order of kn(x) ( n = 0 to 250 )           *
!*       Output:  SK(n) --- kn(x)                                 *
!*                DK(n) --- kn'(x)                                *
!*       Example: x= 10.0                                         *
!*                  n          kn(x)               kn'(x)         *
!*                --------------------------------------------    *
!*                  0     .7131404291D-05    -.7844544720D-05     * 
!*                  1     .7844544720D-05    -.8700313235D-05     *
!*                  2     .9484767707D-05    -.1068997503D-04     * 
!*                  3     .1258692857D-04    -.1451953914D-04     *
!*                  4     .1829561771D-04    -.2173473743D-04     *
!*                  5     .2905298451D-04    -.3572740841D-04     *
!* -------------------------------------------------------------- *
!* REFERENCE: "Fortran Routines for Computation of Special        *
!*             Functions,                                         *
!*             jin.ece.uiuc.edu/routines/routines.html".          *
!*                                                                *
!*                            F90 Release By J-P Moreau, Paris.   *
!*                                   (www.jpmoreau.fr)            *
!******************************************************************
!        PROGRAM MSPHK
!        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!        DIMENSION SK(0:250),DK(0:250)
!        WRITE(*,*)'Please enter n and x '
!        READ(*,*)N,X
!        WRITE(*,30)N,X
!        IF (N.LE.10) THEN
!           NS=1
!        ELSE
!           WRITE(*,*) 'Please enter order step Ns'
!           READ(*,*) NS
!        ENDIF
!        CALL SPHK(N,X,NM,SK,DK)
!        WRITE(*,*)
!        WRITE(*,*)'  n          kn(x)               kn''(x)'
!        WRITE(*,*)'--------------------------------------------'
!        DO 10 K=0,NM,NS
!10         WRITE(*,20)K,SK(K),DK(K)
!20      FORMAT(1X,I3,2D20.10)
!30      FORMAT(3X,'Nmax =',I3,',     ','x =',F6.1)
!        END


        SUBROUTINE MSBESSEL2(N,X,SK)

!       =====================================================
!       Purpose: Compute modified spherical Bessel functions
!                of the second kind, kn(x) and kn'(x)
!       Input :  x --- Argument of kn(x)  ( x ò 0 )
!                n --- Order of kn(x) ( n = 0,1,2,... )
!       Output:  SK(n) --- kn(x)
!                DK(n) --- kn'(x)
!                NM --- Highest order computed
!       =====================================================
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION SK(0:N),DK(0:N)
        PI=3.141592653589793D0
        NM=N
        IF (X.LT.1.0D-60) THEN
           DO 10 K=0,N
              SK(K)=1.0D+300
10            DK(K)=-1.0D+300
           RETURN
        ENDIF
        !SK(0)=0.5D0*PI/X*DEXP(-X)
        SK(0)=DEXP(-X)/X
        SK(1)=SK(0)*(1.0D0+1.0D0/X)
        F0=SK(0)
        F1=SK(1)
        DO 15 K=2,N
           F=(2.0D0*K-1.0D0)*F1/X+F0
           SK(K)=F
           IF (DABS(F).GT.1.0D+300) GO TO 20
           F0=F1
15         F1=F
20      NM=K-1
        DK(0)=-SK(1)
        DO 25 K=1,NM
25         DK(K)=-SK(K-1)-(K+1.0D0)/X*SK(K)
        RETURN
        END

!end of file msphk.f90
