
subroutine kgenrik
use modmain
use modvars
use modplmix
implicit none
! local variables
integer is,nr,nrc,i
integer ir,irc,l
real(8) t1
real(8), allocatable :: ii(:,:,:),kk(:,:,:)
if (allocated(ii)) deallocate(ii)
allocate(ii(nrmtmax,0:lmaxo+10,nspecies))
if (allocated(kk)) deallocate(kk)
allocate(kk(nrmtmax,0:lmaxo+10,nspecies))

if (allocated(ilrmt)) deallocate(ilrmt)
allocate(ilrmt(nrmtmax,0:lmaxo+2,nspecies))
if (allocated(klrmt)) deallocate(klrmt)
allocate(klrmt(nrmtmax,0:lmaxo+2,nspecies))
if (allocated(rr2)) deallocate(rr2)
allocate(rr2(nrmtmax,nspecies))
! generate the radial meshes
do is=1,nspecies
! calculate i_l(\lambda r) on the fine radial mesh
  nr=nrmt(is)
  do i = 1,nr
    rr2(i,is)=rsp(i,is)**2
    t1=rsp(i,is)*klambda
  !  do l=0,lmaxo+2
  !    t1=rsp(i,is)*klambda
  !    call msbessel1(t1,l,ilrmt(i,l,is))
  !    call msbessel2(t1,l,klrmt(i,l,is))
  !  end do
    call msbessel1(lmaxo+10,t1,ii(i,:,is))
    call msbessel2(lmaxo+10,t1,kk(i,:,is))
  end do
end do
ilrmt(:,0:lmaxo+2,:)=ii(:,0:lmaxo+2,:)
klrmt(:,0:lmaxo+2,:)=kk(:,0:lmaxo+2,:)
deallocate(ii,kk)

end subroutine
!EOC

