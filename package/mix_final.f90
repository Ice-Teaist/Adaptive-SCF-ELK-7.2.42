      subroutine mix_final
      use modplmix
      implicit none

      if (allocated(res_m)) deallocate(res_m)
      if (allocated(y_m1)) deallocate(y_m1)
      if (allocated(s_m1)) deallocate(s_m1)
      if (allocated(x_m)) deallocate(x_m)
      if (allocated(sh1y)) deallocate(sh1y)
      if (allocated(yty)) deallocate(yty)
      if (allocated(ilrmt)) deallocate(ilrmt)
      if (allocated(klrmt)) deallocate(klrmt)
      if (allocated(rr2)) deallocate(rr2)
      if (allocated(kvclmt)) deallocate(kvclmt)
      if (allocated(kvclir)) deallocate(kvclir)
      if (allocated(krhomt)) deallocate(krhomt)
      if (allocated(krhoir)) deallocate(krhoir)

      end subroutine
