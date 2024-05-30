      subroutine mix_init(maxsubsp,m,x,mtypes,alpha,rpara)
      use modplmix
      use modmain
      implicit none

      integer, intent(in) :: maxsubsp,m,mtypes
      real*8,intent(in) :: alpha,rpara
      real*8,intent(in) :: x(m)

      integer :: i 
      
      plalpha = alpha
      plmixtype = mtypes 
      plmaxsub = maxsubsp
      repara=rpara

      if (plmixtype == 1) then
         if (allocated(x_m)) deallocate(x_m)
         allocate(x_m(m,2))
         if (allocated(res_m)) deallocate(res_m)
         allocate(res_m(m,2))

      else if (plmixtype == 2) then
         if (allocated(x_m)) deallocate(x_m)
         allocate(x_m(m,2))
         if (allocated(res_m)) deallocate(res_m)
         allocate(res_m(m,2))
         if (allocated(y_m1)) deallocate(y_m1)
         allocate(y_m1(m,1))
         if (allocated(s_m1)) deallocate(s_m1)
         allocate(s_m1(m,1))
         if (allocated(sh1y)) deallocate(sh1y)
         allocate(sh1y(m,1))
         if (allocated(yty)) deallocate(yty)
         allocate(yty(1,1))

      else if (plmixtype == 3) then
         if (allocated(x_m)) deallocate(x_m)
         allocate(x_m(m,2))
         if (allocated(res_m)) deallocate(res_m)
         allocate(res_m(m,2))
         if (allocated(y_m1)) deallocate(y_m1)
         allocate(y_m1(m,plmaxsub))
         if (allocated(s_m1)) deallocate(s_m1)
         allocate(s_m1(m,plmaxsub))
         if (allocated(sh1y)) deallocate(sh1y)
         allocate(sh1y(m,plmaxsub))
         if (allocated(yty)) deallocate(yty)
         allocate(yty(plmaxsub,plmaxsub))

      end if
      call dcopy(m,x(:),1,x_m,1)
      call dcopy(m,x(:),1,x_m(:,2),1)
      res_m(:,1)=0.d0
      y_m1(:,1)=0.d0
      s_m1(:,1)=0.d0
      sh1y=0.d0
      yty(:,:)=0.d0

      if (allocated(krhomt)) deallocate(krhomt)
      allocate(krhomt(npmtmax,natmtot))
      if (allocated(krhoir)) deallocate(krhoir)
      allocate(krhoir(ngtot))

      krhomt(1:npmtmax,1:natmtot)=rhomt(1:npmtmax,1:natmtot)
      krhoir(1:ngtot)=rhoir(1:ngtot)

      end subroutine



