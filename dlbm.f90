module dlbm
implicit none
private

public :: wp, CSQR, dw, ds, pi, twopi

public :: collide, stream 

! Numeric work precision
integer, parameter :: wp = kind(1.0d0)

    ! Density of water
    real(wp), parameter :: dw = 998.0_wp
    ! Density of solids
    real(wp), parameter :: ds = 1510.0_wp
    ! Diffusivity of water
    real(wp), parameter :: diff = 1.0e-10_wp


    ! Canonical speed of sound in LBM
    real(wp) :: csqr = 1._wp/3._wp

    ! Math constants
    real(wp), parameter :: pi = 4.0_wp*atan(1.0_wp)
    real(wp), parameter :: twopi = 2*pi


! 
! --- work arrays ---
!
! PDF distributions
    real(wp), allocatable :: Nf1(:,:), Nf2(:,:)
! Velocities
    real(wp), allocatable :: c(:,:)
! Cell widths
    real(wp), allocatable :: L0(:)
! Solid velocities (defined at edges!)
    real(wp), allocatable :: us(:)
! Moisture concentration
    real(wp), allocatable :: nw(:), X(:), 

! Surface quantities
    real(wp) :: nws, Xs, dnwdx

contains

    pure function solid_density(nw)
        real(wp), intent(in) :: nw
        real(wp) :: solid_density
        solid_density = ds*(1.0_wp - nw(i)/dw)
    end function


subroutine initialize

subroutine collide(N0,N1,N2,c0,c1,c2,nw,X)

    real(wp), intent(inout) :: N0(n),N1(n),N2(n)
    real(wp), intent(inout) :: c0(n),c1(n),c2(n)
    real(wp), intent(out) :: nw(n), X(n)

    real(wp) :: t(0:2)

    ! Relaxation rate (constant diffusivity only)
    omega = dt/(diff/csqr + 0.5_wp*dt)

    do i = 1, nlatt

        ! Calculate pre-collision moments
        m0 = N0(i) + (N1(i) + N2(i))
        m1 = ((c(1,i)*N1(i) - c(2,i)*N2(i)) + c(0,i)*N0(i)) - m0*(c(0,i))

        ! Calculate cell edge positions
        r 
        rl =
        ru = 


        ! Update macroscopic fields
        cell_vol = shell_vol(rl=rloc - 0.5_wp*self%L(i),ru=rloc + 0.5_wp*self%L(i))
        nw(i) = m0/cell_vol
        X(i) = nw(i)/solid_density(nw(i))

        ! Update cell lengths and velocities
        L(i) = L(i) + (usl - usr)*dt
        c(0,i) = 0.5_wp * (usl + usr)
        c(1,i) = abs(c(0,i) + length/dt)
        c(2,i) = abs(c(0,i) - length/dt)

        ! Calculate equilibrium weights
        call weights(c(:,i),t)

        ! Calculate source term
        src = 0.5_wp*dt*cell_vol*(diff/r)*(gradnw)
        m0 = m0 + src
    
        ! Relaxation
        m1 = m1*(1.0_wp - omega)

        ! Post-collision PDFs
        associate(csum => (c(1,i) + c(2,i)))
            N0(i) = t(1)*m0
            N1(i) = t(2)*m0 + m1/csum
            N2(i) = t(3)*m0 - m1/csum
        end associate
    
    end do

contains

    subroutine weights(c,t)
        real(wp), intent(in) :: c(0:2)
        real(wp), intent(out) :: t(0:2)

        associate(num => csqr/(c(1) + c(2)))
            t(0) = 1.0_wp - csqr/(c(2)*c(3))
            t(1) = num/c(2)
            t(2) = num/c(3)
        end associate
    end subroutine

end subroutine collide

! Streaming and particle exchange
!
! Since we have two copies, we can write into the result array
! without fear of over-writing needed values.
!
subroutine streaming_and_particle_exchange(n,c,N0,N1,N2,N0p,N1p,N2p,us)
    integer, intent(in) :: n

! Particle velocites
    real(wp), intent(in) :: c(0:2,n)
! Post-collision PDFs
    real(wp), intent(in) :: N0(n),N1(n),N2(n)
! After streaming (p for prime -> N')
    real(wp), intent(out) :: N0p(n),N1p(n),N2p(n)

! Solid velocity
    real(wp), intent(out) :: us(0:n)

    ! Copy rest particles
    N0p = N0

    ! Bounceback at the center
    N1p(1) = N2(1)

    do i = 2, n

        csum = c(1,i) + c(2,i)

        ! Weights in particle
        a = (c(1,i) - c(2,i))/csum
        b1 = 2.0_wp*c(1,i)/csum
        b2 = 2.0_wp*c(2,i)/csum

        ! Population of the right cell
        N1p(i)   = b1*N2(i) + a*N1(i-1)
        N2p(i-1) = b2*N1(i-1) - a*N2(i)

        ! Solid velocity at left cell edge
        rw = rloc - 0.5_wp*self%L(i)
        us(i-1) = -(N2(i-1) - N1p(i))/(dw*(twopi*rw)*dt)

    end do

    ! Boundary condition at surface
    call surface_boundary(J)

    N2p(n) = N1(n) - J*dt

    nws = nweq
    Xs  = nws / solid_density(nws)

    ! Solid velocity at last/right cell edge
    edgepos
    us(n) = (N1(n) - N2p(n))/(dw*twopi*edgepos*self%dt)

end subroutine stream




end module