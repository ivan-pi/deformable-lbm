program dlbm_main

use dlbm, only: wp, dw, ds, collide, stream

use dlbm_cmd
use dlbm_io

implicit none 

type(cmd_settings) :: cmd

! default settings
cmd = 

cmd = read_from_namelist(cmd)
cmd = parse_command_line(cmd)

! priority
! 1) defaults
! 2) namelist
! 3) command line variables

call dlbm_initialize(L0,X0,n,csqr,dt)

call dlbm_init_output_units(mean=.true.,full=.true.,
    regular_regime=-1,penetration_period=-1)

! Main timeloop
do step = 1, cmd%nsteps

    ! Periodic output
    if (is_output_step(step,cmd%nout)) then

    end if

    call dlbm_collide()
    call dlbm_stream(nwb=cmd%nwb)

end do

call dlbm_close_output_units()

call print_statistics()

contains

    logical function is_output_step(step,output_step)
        integer, intent(in) :: step, output_step
        is_output_step = mod(step,cmd%nout) == 0
    end function

    subroutine print_statistics(unit)
        integer, intent(in), optional :: unit
        use, intrinsic :: iso_fortran_env, only: output_unit
        
        integer :: unit_

        unit_ = output_unit
        if (present(unit)) unit_ = unit

        write(unit_,'(A)') "Final size [m] = ", 

    end subroutine

    subroutine compare_against_theoretical_length(p,length,X0,Xeq)
        integer, intent(in) :: p 
        real(wp), intent(in) :: length

        real(wp) :: volume_ratio, Leq

        volume_ratio = (1 + (ds/dw)*Xeq)/(1 + (ds/dw)*X0)

            print *, "Length/radius    = ", length

        select case(p)
        case(0)
            Leq = L0*volume_ratio
            print *, "L_eq (slab)      = ", Leq
        case(1)
            Leq = L0*sqrt(volume_ratio)
            print *, "R_eq (cylinder)  = ", Leq
        case(2)
            Leq = L0*(volume_ratio)**(1._wp/3._wp)
            print *, "R_eq (sphere)    = ", Leq
        case default
            error stop "invalid value for p, p = ", p
        end select
            
            print *, "Rel. error       = ", abs(length - Leq)/Leq

    end subroutine

end program