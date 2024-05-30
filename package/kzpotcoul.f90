
subroutine kzpotcoul(nr,nri,np,npi,ld1,rl,ngdg,igf,ngp,gpc,gclgp,ld2,jlgprmt, &
 ylmgp,sfacgp,zrhoir,ld3,zvclmt,zvclir)
use modplmix
use modmain
implicit none
! arguments
integer, intent(in) :: nr(nspecies),nri(nspecies)
integer, intent(in) :: np(nspecies),npi(nspecies)
integer, intent(in) :: ld1
real(8), intent(in) :: rl(ld1,-lmaxo-1:lmaxo+2,nspecies)
integer, intent(in) :: ngdg(3),igf(*),ngp
real(8), intent(in) :: gpc(ngp),gclgp(ngp)
integer, intent(in) :: ld2
real(8), intent(in) :: jlgprmt(0:lnpsd,ld2,nspecies)
complex(8), intent(in) :: ylmgp(lmmaxo,ngp),sfacgp(ld2,natmtot)
complex(8), intent(in) :: zrhoir(*)
integer, intent(in) :: ld3
complex(8), intent(inout) :: zvclmt(ld3,natmtot)
complex(8), intent(out) :: zvclir(*)
! local variables
integer is,ia,ias
integer ir,l,m,lm,it
integer ig,jg,i,j
real(8) t1,t2,t3
real(8) tt1,tt2,tt3
complex(8) z1,z2
! automatic arrays
real(8), allocatable :: rmti(:,:)
real(8), allocatable :: rmtk(:,:)
complex(8) qlm(lmmaxo,natmtot)
complex(8) zl(0:lmaxo),zlm(lmmaxo)
complex(8) zhmt(ld3)
! external functions
real(8), external :: factnm
!kerker
complex(8) :: zzvclir(ngtot)



it=max(lmaxo+3,lnpsd)
allocate(rmti(0:it,nspecies))
allocate(rmtk(0:it,nspecies))
! compute (R_mt)^l
do is=1,nspecies
!  rmtl(0,is)=1.d0
!  do l=0,it
!    rmtl(l,is)=rmtl(l-1,is)*rmt(is)
  call msbessel1(it,rmt(is)*klambda,rmti(:,is))
  call msbessel2(it,rmt(is)*klambda,rmtk(:,is))
!  end do
end do
! compute the multipole moments from the muffin-tin potentials
t1=1.d0/(fourpi*klambda)
do ias=1,natmtot
  is=idxis(ias)
  i=np(is)-lmmaxo
  lm=0
  do l=0,lmaxo
!    t2=t1*dble(2*l+1)*rmtl(l+1,is)
    t2=t1/rmtk(l,is)*dble(factnm(2*l+1,2))/(klambda**l)
    do m=-l,l
      lm=lm+1
      i=i+1
      qlm(lm,ias)=t2*zvclmt(i,ias)
    end do
  end do
end do
! Fourier transform density to G-space and store in zvclir
call zcopy(ngdg(1)*ngdg(2)*ngdg(3),zrhoir,1,zvclir,1)
call zfftifc(3,ngdg,-1,zvclir)
! subtract the multipole moments of the interstitial charge density
do is=1,nspecies
  do l=0,lmaxo
    zl(l)=fourpi*zil(l)*dble(factnm(2*l+1,2))/(klambda**l)
  end do
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    zlm(:)=0.d0
    do ig=1,ngp
      jg=igf(ig)
      if (gpc(ig).gt.epslat) then
        z1=zvclir(jg)*sfacgp(ig,ias)*rmt(is)**2/(gpc(ig)**2+klambda**2)
       ! zlm(1)=zlm(1)+jlgprmt(1,ig,is)*z1*zl(0)*y00
        zlm(1)=zlm(1)+z1*zl(0)*y00*(gpc(ig)*rmti(0,is)*jlgprmt(1,ig,is)+klambda*rmti(1,is)*jlgprmt(0,ig,is))
        lm=1
        do l=1,lmaxo
       !   z2=jlgprmt(l+1,ig,is)*z1*zl(l)
          z2=z1*zl(l)*(gpc(ig)*rmti(l,is)*jlgprmt(l+1,ig,is)+klambda*rmti(l+1,is)*jlgprmt(l,ig,is))
          do m=-l,l
            lm=lm+1
            zlm(lm)=zlm(lm)+z2*conjg(ylmgp(lm,ig))
          end do
        end do
      else
       ! t1=(fourpi/3.d0)*rmtl(3,is)*y00
        t1=(fourpi)*rmt(is)**2*rmti(1,is)*y00/klambda
        zlm(1)=zlm(1)+t1*zvclir(jg)
      end if
    end do
    qlm(:,ias)=qlm(:,ias)-zlm(:)
  end do
end do
! find the smooth pseudocharge within the muffin-tin whose multipoles are the
! difference between the real muffin-tin and interstitial multipoles
t1=(fourpi/omega)*(klambda**lnpsd)!*factnm(2*lnpsd+1,2)
do ias=1,natmtot
  is=idxis(ias)
  lm=0
  tt1=t1/rmti(lnpsd,is)
  do l=0,lmaxo
    t2=tt1/(factnm(2*l+1,2))
    z1=t2*zilc(l)
    do m=-l,l
      lm=lm+1
      zlm(lm)=z1*qlm(lm,ias)
    end do
  end do
! add the pseudocharge and real interstitial densities in G-space
  do ig=1,ngp
    jg=igf(ig)
    if (gpc(ig).gt.epslat) then
      !t2=gpc(ig)*rmt(is)
      !t3=1.d0/t2**lnpsd
      t3=1.d0/gpc(ig)**lnpsd
      z1=t3*zlm(1)*y00
      lm=1
      do l=1,lmaxo
        lm=lm+1
        z2=zlm(lm)*ylmgp(lm,ig)
        do m=1-l,l
          lm=lm+1
          z2=z2+zlm(lm)*ylmgp(lm,ig)
        end do
        !t3=t3*t2
        t3=t3*gpc(ig)
        z1=z1+t3*z2
      end do
      z2=jlgprmt(lnpsd,ig,is)*conjg(sfacgp(ig,ias))
      zvclir(jg)=zvclir(jg)+z1*z2
    else
      !t2=y00/factnm(2*lnpsd+1,2)/rmti(lnpsd,is)
      t2=y00*(fourpi)/omega*(klambda*rmt(is))**lnpsd/factnm(2*lnpsd+1,2)/rmti(lnpsd,is)
      !zvclir(jg)=zvclir(jg)+t2*zlm(1)
      zvclir(jg)=zvclir(jg)+t2*qlm(1,ias)
    end if
  end do
end do
zzvclir(1:ngtot)=zvclir(1:ngtot)
! solve Poisson's equation in G+p-space for the pseudocharge
do ig=1,ngp
  jg=igf(ig)
!  write(2111,*)zvclir(jg)
  if (Resta > 0) then
    if (gclgp(ig) < 0.000001) then
      zvclir(jg)=zvclir(jg)*Resta*gclgp(ig+1)
!      write(2116,*)zvclir(jg),ig
!      write(2116,*)
    else
      t1=1/(1/gclgp(ig)+(klambda**2/fourpi))/gclgp(ig)
      t1=max(t1,Resta)
      zvclir(jg)=t1*zvclir(jg)*gclgp(ig)
!      write(2116,*)ig,zvclir(jg)
    end if
  else
    t1=1/(1/gclgp(ig)+(klambda**2/fourpi))  
    zvclir(jg)=t1*zvclir(jg)
  end if
end do
!kerker
do ig=1,ngp
  jg=igf(ig)
  t1=1/(1/gclgp(ig)+(klambda**2/fourpi))  
  zzvclir(jg)=t1*zzvclir(jg)
end do
! match potentials at muffin-tin boundary by adding homogeneous solution
do ias=1,natmtot
  is=idxis(ias)
! find the spherical harmonic expansion of the interstitial potential at the
! muffin-tin radius
  zlm(:)=0.d0
  do ig=1,ngp
    z1=fourpi*zzvclir(igf(ig))*sfacgp(ig,ias)
    zlm(1)=zlm(1)+jlgprmt(0,ig,is)*z1*y00
    lm=1
    do l=1,lmaxo
      z2=jlgprmt(l,ig,is)*z1*zil(l)
      do m=-l,l
        lm=lm+1
        zlm(lm)=zlm(lm)+z2*conjg(ylmgp(lm,ig))
      end do
    end do
  end do
! calculate the homogenous solution
  i=np(is)-lmmaxo
  lm=0
  do l=0,lmaxi
    t1=1.d0/rmti(l,is)
    do m=-l,l
      lm=lm+1
      i=i+1
      z1=t1*(zlm(lm)-zvclmt(i,ias))
      j=lm
      do ir=1,nri(is)
        zhmt(j)=z1*ilrmt(ir,l,is)
        j=j+lmmaxi
      end do
      do ir=nri(is)+1,nr(is)
        zhmt(j)=z1*ilrmt(ir,l,is)
        j=j+lmmaxo
      end do
    end do
  end do
  do l=lmaxi+1,lmaxo
    t1=1.d0/rmti(l,is)
    do m=-l,l
      lm=lm+1
      i=i+1
      z1=t1*(zlm(lm)-zvclmt(i,ias))
      j=npi(is)+lm
      do ir=nri(is)+1,nr(is)
        zhmt(j)=z1*ilrmt(ir,l,is)
        j=j+lmmaxo
      end do
    end do
  end do
  zvclmt(1:np(is),ias)=zvclmt(1:np(is),ias)+zhmt(1:np(is))
end do
!open(2116,file='kflm')
!write(2116,*)zvclmt(1:ngp,1)
!close(2116)
! Fourier transform interstitial potential to real-space
call zfftifc(3,ngdg,1,zvclir)
deallocate(rmti,rmtk)
end subroutine
!EOC

