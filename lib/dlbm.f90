!
! Deformable Lattice Boltzmann Module
!
module dlbm

   use, intrinsic :: iso_c_binding, only: &
      int => c_int, 
      dp => c_double

   implicit none ! ~ Carpe Diem ~
   private

   public :: dlbm_init
   public :: dlbm_step
   public :: dlbm_integrate

   ! Density of water and solids
   real(dp), parameter :: dw = 998.0_dp
   real(dp), parameter :: ds = 1510.0_dp

   ! Math constants
   real(dp), parameter :: pi = 4.0_dp*atan(1.0_dp)
   real(dp), parameter :: twopi = 2*pi

contains

! Initialize LBM method
subroutine dlbm_init(n,csqr,dt,X0,L0,X,us,r,Npdf) bind(c)
    integer(int), intent(in), value :: n
    real(dp), intent(in), value :: csqr, dt, X0, L0
    real(dp), intent(out) :: X(n), us(n), r(n), Npdf(n,3) 

end subroutine

! Time-stepping
subroutine dlbm_step(nsteps,D,csqr,n,X,us,r,N) bind(c)

   integer(int), intent(in), value :: nsteps, n
   real(dp), intent(in), value :: D, csqr
   real(dp), intent(out) :: X(n), us(n), r(n)
   real(dp), intent(inout) :: Npdf(n,3)

   ! Automatic buffers
   real(dp) :: c(0:2,n), Ntmp(n,3)

   ! Initialize particle velocites

   c(0,:) = us

   do step = 0, nsteps - 1

      if (mod(step,2) == 0) then
         call stream_collide(Npdf,Ntmp)
      else
         call stream_collide(Ntmp,Npdf)
      end if
   
   end do

   ! Results for output

   us = c(0,:)

contains

   subroutine stream_collide(n,c,Nold,Nnew)
      integer(int), intent(in) :: n
      real(dp), intent(inout) :: c(0:2,n)
      real(dp), intent(in) :: Nold(n,0:2)
      real(dp), intent(out) :: Nnew(n,0:2)

      use dlbm_kernels, only: collide, &
         stream => streaming_and_particle_exchange
         
         call stream( &
            n = n, &
            c = c, &
            N0  = Nold(:,0), N1  = Nold(:,1), N2  = Nold(:,2), &
            N0p = Nnew(:,0), N1p = Nnew(:,1), N2p = Nnew(:,2), &
            us = us)

         call collide( &
            n = n,
            N0 = Nnew(:,0), N1 = Nnew(:,1), N2 = Nnew(:,2), &
            c=c,nw=nw)

   end subroutine

end subroutine


! Midpoint integration rule for cell-centered quantities
function dlbm_integrate(n,rw,X,R) result(avg) bind(c)
   integer, intent(in) :: n
   real(dp), intent(in) :: rw(n), X(n)
   real(dp), intent(out) :: R

   real(dp) :: avg
   integer :: i

   avg = 0
   R = 0

   rcenter = 0.5_dp*rw(1)
   do i = 1, n-1
      avg = avg + rw(i)*(rcenter*X(i))
      R = R + rw(i)
      rcenter = rcenter + 0.5_dp*(rw(i) + rw(i+1))
   end do

   avg = avg + rw(n)*(rcenter*X(n))
   R = R + rw(n)

end function

end module