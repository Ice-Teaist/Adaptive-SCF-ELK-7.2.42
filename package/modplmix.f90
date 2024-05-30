      module modplmix
! Pulay mixing variable      
      !Maximum dimension of subspace; mixed type 
      integer :: plmaxsub, plmixtype, premix
      !alpha
      real*8 :: plalpha
      !Nonsingular parameter
      real*8 :: repara
      !x vector  array
      real*8, allocatable :: x_m(:,:)
      !r vector  array
      real*8, allocatable :: res_m(:,:)
      !\delta x vector  array
      real*8, allocatable :: s_m1(:,:)
      !\delta r vector  array
      real*8, allocatable :: y_m1(:,:)
      !Y^TY, S+H_1Y  array
      real*8, allocatable :: sh1y(:,:)
      real*8, allocatable :: yty(:,:)
          
! Kerker precondition variable 
      real(8) :: klambda
      real(8) :: Resta
      ! muffin-tin and interstitial Coulomb potential
      real(8), allocatable :: kvclmt(:,:),kvclir(:)
      real(8), allocatable :: krhomt(:,:),krhoir(:)
      ! i_l(r) on fine radial mesh
      real(8), allocatable :: ilrmt(:,:,:)
      ! k_l(r) on fine radial mesh
      real(8), allocatable :: klrmt(:,:,:)
      ! k_l(r) on fine radial mesh
      real(8), allocatable :: rr2(:,:)





      end module
