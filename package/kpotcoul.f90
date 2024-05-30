
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: potcoul
! !INTERFACE:
subroutine kpotcoul(m,res,kre)
! !USES:
use modmain
use modplmix
use modomp
! !DESCRIPTION:
!   Calculates the Coulomb potential of the real charge density stored in the
!   global variables {\tt rhomt} and {\tt rhoir} by solving Poisson's equation.
!   These variables are coverted to complex representations and passed to the
!   routine {\tt zpotcoul}.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!EOP
!BOC
implicit none
integer, intent(in) :: m
real(8), intent(in) :: res(m)
real(8), intent(out) :: kre(m)
! local variables
integer is,ias,nthd
integer nr,nri,ir,i,j,k
real(8) :: ts0,ts1,pot
! allocatable arrays
complex(8), allocatable :: zrhomt(:,:),zrhoir(:)
complex(8), allocatable :: zvclmt(:,:),zvclir(:)
real(8), allocatable :: kdmt(:,:),kdir(:)
call timesec(ts0)
allocate(zrhomt(npmtmax,natmtot))
allocate(kvclmt(npmtmax,natmtot))
allocate(kvclir(ngtot))

allocate(kdmt(npmtmax,natmtot))
allocate(kdir(ngtot))

j=1
do i = 1,natmtot
  do k = 1,npmtmax
    kdmt(k,i)=-res(j)/fourpi
    j=j+1
  end do
end do
do i = 1,ngtot
  kdir(i)=-res(j)/fourpi
  j=j+1
end do


! Generating the first and second kind of modified spherical Bessel function
call kgenrik
! convert real muffin-tin charge density to complex spherical harmonic expansion
call holdthd(natmtot,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(is) &
!$OMP NUM_THREADS(nthd)
do ias=1,natmtot
  is=idxis(ias)
  call rtozfmt(nrmt(is),nrmti(is),kdmt(:,ias),zrhomt(:,ias))
end do
!$OMP END PARALLEL DO
call freethd(nthd)
! solve the complex Poisson's equation in the muffin-tins
allocate(zvclmt(npmtmax,natmtot))
call kgenzvclmt(nrmt,nrmti,nrmtmax,ilrmt,klrmt,rr2,wprmt,npmtmax,zrhomt,zvclmt)
deallocate(zrhomt)
allocate(zrhoir(ngtot),zvclir(ngtot))
! store real interstitial charge density in complex array
zrhoir(:)=kdir(:)
! solve Poisson's equation in the entire unit cell
call kzpotcoul(nrmt,nrmti,npmt,npmti,nrmtmax,rlmt,ngridg,igfft,ngvec,gc,gclg, &
 ngvec,jlgrmt,ylmg,sfacg,zrhoir,npmtmax,zvclmt,zvclir)
! convert complex muffin-tin potential to real spherical harmonic expansion
call holdthd(natmtot,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(is) &
!$OMP NUM_THREADS(nthd)
do ias=1,natmtot
  is=idxis(ias)
  call ztorfmt(nrmt(is),nrmti(is),zvclmt(:,ias),kvclmt(:,ias))
end do
!$OMP END PARALLEL DO
call freethd(nthd)
! store complex interstitial potential in real array
kvclir(:)=dble(zvclir(:))

j=1
do i = 1,natmtot
  do k = 1,npmtmax
    kre(j)=kvclmt(k,i)
    j=j+1
  end do
end do
do i = 1,ngtot
  kre(j)=kvclir(i)
  j=j+1
end do
deallocate(kdmt,kdir)
deallocate(zrhoir,zvclmt,zvclir)
deallocate(kvclmt,kvclir)
call timesec(ts1)
pot=ts1-ts0
end subroutine
!EOC

