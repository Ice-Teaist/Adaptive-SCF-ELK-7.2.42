      subroutine pmix_scf(nums,m,x,d,ttra)
      use modplmix
      implicit none
!     Describe:
!     \textbf{x}_{m+1}&=\textbf{x}_m +H_m\textbf{r}(\textbf{x}_m) \\
!     & = \textbf{x}_m +H_1\textbf{r}(\textbf{x}_m) âˆ?(s_{mâˆ?}+H_1y_{mâˆ?})\\
!     &\quad\times( y_{m-1}^Ty_{m-1})^{-1}y^T_{m-1}\textbf{r}(\textbf{x}_m)

      integer, intent(in) :: nums,m
      real*8,intent(inout) :: x(m)
      real*8,intent(out) :: d
      integer, intent(out) :: ttra

      integer :: numscf,i,j,k,info,rank,lwork,liwork,jc,kp,kc,mm
      integer:: ipiv(plmaxsub)
      real*8,allocatable:: work(:)
      real*8, allocatable :: singular(:)
      real*8, allocatable :: svd_u(:,:),svd_vt(:,:),svd_s(:,:)
      real*8, allocatable :: svd_vs(:,:),svd_y(:,:)
      real*8, allocatable :: yyy(:,:)
      real*8 :: dummy(1,1)
      real*8 :: lw(1)
      integer :: liw(1)
      real*8 :: t1,t2
      real*8, external :: dnrm2
      real*8 :: c(plmaxsub),gamma(plmaxsub),beta(plmaxsub,plmaxsub)

      real*8 :: kre(m)        
      real*8, allocatable ::  indic(:,:)        
      Integer, Parameter               :: nb = 64, nin = 5, nout = 6
!     .. Local Scalars ..
      Real (Kind=8)               :: abnrm
      Integer                          :: ihi, ilo, ldvl,ldvr 
!     .. Local Arrays ..
      Real (Kind=8), Allocatable  :: rconde(:), rcondv(:), scale(:)
      Real (Kind=8), Allocatable  :: vl(:,:), vr(:,:),wi(:), wr(:)
      Integer, Allocatable             :: iwork(:)






      d=1.d0

      !simple mixing
      !if (plmixtype == 1.or.nums==1 ) then
      if (plmixtype == 1 ) then
         ! previous index modulo x_m and res_m
         kp=mod(nums-1,2)+1
         ! current index modulo x_m and res_m
         kc=mod(nums,2)+1
 
         !residual
         do i = 1,m    
            res_m(i,kc)=x(i)-x_m(i,kp)
         end do

         ! residual module
         d=sum(res_m(:,kc)**2)
         d=sqrt(d/dble(m))

         if (nums == 1) then
            do i = 1,m
               x(i) = x_m(i,kp)+plalpha * res_m(i,kc)
            end do
         else
            call kpotcoul(m,res_m(:,kc),kre)
            do i = 1,m
               x(i) = x_m(i,kp)+plalpha*(res_m(i,kc)+klambda**2*kre(i))
            end do
         end if

         !x_{m+1}
         do i = 1,m
            x(i) = x_m(i,kp)+plalpha * res_m(i,kc)
         end do        
         !x_{m+1} to x_m array
         call dcopy(m,x,1,x_m(:,kc),1)        


      !Anderson mixing    
      else if (plmixtype == 2 .and. nums >1 ) then
         ! previous index modulo x_m and res_m
         kp=mod(nums-1,2)+1
         ! current index modulo x_m and res_m
         kc=mod(nums,2)+1
 
         !residual
         do i = 1,m    
            res_m(i,kc)=x(i)-x_m(i,kp)
         end do

         ! residual module
         d=sum(res_m(:,kc)**2)
         d=sqrt(d/dble(m))
         !y_m1
         do i = 1,m    
            y_m1(i,1)=res_m(i,kc)-res_m(i,kp)
         end do
         t1=dnrm2(m,y_m1(:,1),1)
         if (t1.gt.1.d-8) t1=1.d0/t1
         !Normalized y_m1
         call dscal(m,t1,y_m1(:,1),1)

         !s_m1
         do i = 1,m
            s_m1(i,1)=t1*(x_m(i,kp)-x_m(i,kc))
         end do

         !s_m1 + H_1 * y_m1
         do i=1,m
            sh1y(i,1)=plalpha*y_m1(i,1)+s_m1(i,1)
         end do

         !inverse
         yty(1,1)=dot_product(y_m1(:,1),y_m1(:,1))
         if (yty(1,1).gt.1.d-8) yty(1,1)=1.d0/yty(1,1)
         !{y_m1}^{T}*res(x_m)
         c(1)=dot_product(y_m1(:,1),res_m(:,kc))
         !(y_m1^Ty_m1)^{-1}{y_m1}^T*res(x_m)     
         gamma(1)=c(1)*yty(1,1)

         !x_{m+1}
         do i = 1,m
            x(i) = x_m(i,kp)+plalpha * res_m(i,kc)
         end do        
         call daxpy(m,-gamma(1),sh1y(:,1),1,x,1)
         !x_{m+1} to x_m array
         call dcopy(m,x,1,x_m(:,kc),1)        

      !Pulay mixing
      !else if (plmixtype == 3 .and. nums>1 ) then
      else if (plmixtype == 3  ) then

         ! current subspace dimension
         mm=min(nums+1,plmaxsub)
         ! current index modulo subspace
         jc=mod(nums,mm)+1
         ! previous index modulo x_m and res_m
         kp=mod(nums-1,2)+1
         ! current index modulo x_m and res_m
         kc=mod(nums,2)+1

         !residual
         do i = 1,m    
            res_m(i,kc)=x(i)-x_m(i,kp)
         end do

         ! residual module
         d=sum(res_m(:,kc)**2)
         d=sqrt(d/dble(m))
         !y_m1
         do i = 1,m    
            y_m1(i,jc)=res_m(i,kc)-res_m(i,kp)
         end do
         t1=dnrm2(m,y_m1(:,jc),1)
         if (t1.gt.1.d-8) t1=1.d0/t1
         !Normalized y_m1
         call dscal(m,t1,y_m1(:,jc),1)

         !s_m1
         do i = 1,m
            s_m1(i,jc)=t1*(x_m(i,kp)-x_m(i,kc))
         end do

         !s_m1 + H_1 * y_m1
         !call kpotcoul(m,y_m1(:,jc),kre)
         do i=1,m
            sh1y(i,jc)=plalpha*y_m1(i,jc)+s_m1(i,jc)
         !   sh1y(i,jc)=plalpha*(y_m1(i,jc)+klambda**2*kre(i))+s_m1(i,jc)
         end do

         ldvl = mm
         ldvr = mm
         lwork = (2+nb)*mm
         Allocate(rconde(mm),rcondv(mm),scale(mm),vl(mm,mm),vr(mm,mm), &
         wi(mm),wr(mm),iwork(2*mm-2))
         !The construction of the a posteriori indicator
         if (allocated(indic)) deallocate(indic)
         allocate(indic(mm,mm))
         !call dgemm('T','N',mm,mm,m,1.0,y_m1(:,:),m,-sh1y(:,:),m, &
         !     0.0,indic,mm)
         indic=0
         do i = 1, mm
            do j = 1,mm
               indic(i,j)=dot_product(y_m1(:,i),-sh1y(:,j))
            end do
         end do
         !Use routine workspace query to get optimal workspace.
         lwork = -1
         !The NAG name equivalent of dgeevx is f08nbf
         Call dgeevx('N','N','N','N',mm,indic,mm,wr,wi,vl,ldvl,vr,ldvr,ilo, &
         ihi,scale,abnrm,rconde,rcondv,dummy(1,1),lwork,iwork,info)

         !Make sure that there is enough workspace for block size nb.
         lwork = max((nb+2)*mm,nint(dummy(1,1)))
         if (allocated(work)) deallocate(work)
         Allocate (work(lwork))
         !Solve the eigenvalue problem
         !The NAG name equivalent of dgeevx is f08nbf
         Call dgeevx('N','N','N','N',mm,indic,mm,wr,wi,vl,ldvl,vr,ldvr,ilo, &
         ihi,scale,abnrm,rconde,rcondv,dummy(1,1),lwork,iwork,info)
         t2=10
         do i = 1,mm
            t2=min(ABS(wr(i)+1),t2)
         end do
         if (t2<0.2 ) then
           ttra=nums 
         end if


         if (allocated(indic)) deallocate(indic)
         if (allocated(work)) deallocate(work)
         deallocate(rconde,rcondv,scale,vl,vr, &
         wi,wr,iwork)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         if (repara == 0) then
            
            if (allocated(yyy)) deallocate(yyy)
            allocate(yyy(m,plmaxsub))
            if (allocated(singular)) deallocate(singular)
            allocate(singular(mm))
            if (allocated(svd_u)) deallocate(svd_u)
            allocate(svd_u(m,mm))
            if (allocated(svd_vt)) deallocate(svd_vt)
            allocate(svd_vt(mm,m))
            if (allocated(svd_s)) deallocate(svd_s)
            allocate(svd_s(mm,mm))
            if (allocated(svd_vs)) deallocate(svd_vs)
            allocate(svd_vs(mm,mm))
            if (allocated(svd_y)) deallocate(svd_y)
            allocate(svd_y(mm,m))
           
            !yyy(:,:)=y_m1(:,:)
            do i = 1,plmaxsub
               call dcopy(m,y_m1(:,i),1,yyy(:,i),1)
            end do
            !Use routine workspace query to get optimal workspace.
            lwork = -1
            call dgesvd('S','S',m,mm,yyy,m,singular,svd_u,m,&
               svd_vt,mm,dummy,lwork,info)
            !Make sure that there is enough workspace
            lwork = max(m+4*mm+64*(m+mm),nint(dummy(1,1)))
            if (allocated(work)) deallocate(work)
            allocate (work(lwork))
            !Compute the singular values and left and right singular vectors
            !of y_m1.
            call dgesvd('S','S',m,mm,yyy,m,singular,svd_u,m,&
               svd_vt,mm,work,lwork,info)
            svd_s(:,:)=0
            do i = 1,mm
               if (singular(i).gt.1.d-8) svd_s(i,i)=1.d0/singular(i)
            end do  
           !{y_m1}^{-1}=V\sigma^+U^T 
            call dgemm('T','N',mm,mm,mm,1.d0,svd_vt(:,:), &
               mm,svd_s,mm,0.d0,svd_vs(:,:),mm)
            call dgemm('N','T',mm,m,mm,1.d0,svd_vs,mm, &
               svd_u(:,:),m,0.d0,svd_y(:,:),mm)

            !{y_m1}^{-1}*res(x_m)
            call dgemv('N',mm,m,1.d0,svd_y(:,:),mm,res_m(:,kc),1,0.d0,&
               c(:),1)
            !(y_m1^Ty_m1)^{-1}{y_m1}^T*res(x_m)     
            do i=1,mm
               gamma(i)=c(i)
            end do

            if (allocated(yyy)) deallocate(yyy)
            if (allocated(singular)) deallocate(singular)
            if (allocated(svd_u)) deallocate(svd_u)
            if (allocated(svd_vt)) deallocate(svd_vt)
            if (allocated(svd_s)) deallocate(svd_s)
            if (allocated(svd_vs)) deallocate(svd_vs)
            if (allocated(svd_y)) deallocate(svd_y)
         
         else
            
            if (allocated(work)) deallocate(work)
            allocate(work(plmaxsub))

            !{y_m1}^{T}y_m1
            do k=1,mm
               yty(k,jc)=dot_product(y_m1(:,jc),y_m1(:,k))
               yty(jc,k)=yty(k,jc)
            end do      
            !beta(:,:)=yty(:,:)
            do i = 1,plmaxsub
               call dcopy(plmaxsub,yty(:,i),1,beta(:,i),1)
            end do
            do k=1,mm
               beta(k,k)=beta(k,k)+repara**2
            end do
            ! invert 
            call dgetrf(mm,mm,beta,plmaxsub,ipiv,info)
            if (info==0) call dgetri(mm,beta,plmaxsub,ipiv,work,mm,info)
            if (info.ne.0) then
               write(*,*)
               write(*,'("Error(mix): could not invert matrix")')
               write(*,*)
               stop
            end if        

            !{y_m1}^{T}*res(x_m)
            call dgemv('C',m,mm,1.d0,y_m1(:,:),m,res_m(:,kc),1,0.d0,&
               c(:),1)

            !(y_m1^Ty_m1)^{-1}{y_m1}^T*res(x_m)     
            do i=1,mm
               if(abs(c(i))<1E-6) c(i)=0
               gamma(i)=0.d0
               gamma(i)=dot_product(c(:),beta(:,i))
            end do
         end if
         
         kre=0
         !x_{m+1}
         !call kpotcoul(m,res_m(:,kc),kre)
         do i = 1,m
            x(i) = x_m(i,kp)+plalpha * res_m(i,kc)
         !   x(i) = x_m(i,kp)+plalpha * (res_m(i,kc)+klambda**2*kre(i))
         end do        
         do i=1,mm
            call daxpy(m,-gamma(i),sh1y(:,i),1,x,1)
         end do
         !x_{m+1} to x_m array
         call dcopy(m,x,1,x_m(:,kc),1)        

         if (allocated(work)) deallocate(work)


      end if 


      end subroutine
