!-----------------------------------------------------------------------
!Module: mechanics
!-----------------------------------------------------------------------
!! By: David Wilkins
!!
!! The subroutines and functions in this module are where the equations
!! for calculating the physical properties of the system are defined.
!!----------------------------------------------------------------------
!! Included subroutines:
!!
!! calculate_energy
!!----------------------------------------------------------------------
!! Included functions:
!!
!! planets_ode
!-----------------------------------------------------------------------
module mechanics
use types
implicit none
private
public :: planets_ode, calculate_energy

contains


!-----------------------------------------------------------------------
!! function: planets_ode
!-----------------------------------------------------------------------
!! By: David Wilkins
!!
!! planets_ode is where the system of eight differential equations is defined.
!! The position and velocity values are sent to planets_ode, where the new
!! values of the differential equations are determined.
!!----------------------------------------------------------------------
!! Input:
!!
!! q        real(dp)    array containing current values of position and velocity
!! t        real(dp)    time value
!! work     real(dp)    array containing masses
!!----------------------------------------------------------------------
!! Output:
!!
!! f        real(dp)    array containing values of eight differential equations
!-----------------------------------------------------------------------
function planets_ode(q, t, work) result(f)
    implicit none
    real(dp), intent(in) :: q(:), t, work(:)
    real(dp), allocatable :: f(:)
    real(dp) :: r1, r2, r12
    ! This is the function that will be sent to 
    ! solve_runge_kutta_4 as an argument.

    allocate(f(size(q)))

    ! defining distance terms:
    r1 = sqrt(q(1)**2.0_dp + q(2)**2.0_dp)
    r2 = sqrt(q(3)**2.0_dp + q(4)**2.0_dp)
    r12 = sqrt((q(1)-q(3))**2.0_dp + (q(2)-q(4))**2.0_dp)

    ! system of differential equations defined here:
    f(1) = q(5)
    f(2) = q(6)
    f(3) = q(7)
    f(4) = q(8)
    f(5) = (-work(1) / (r1**3.0_dp)) * q(1) - (work(3) / (r12**3.0_dp)) * (q(1) - q(3))
    f(6) = (-work(1) / (r1**3.0_dp)) * q(2) - (work(3) / (r12**3.0_dp)) * (q(2) - q(4))
    f(7) = (-work(1) / (r2**3.0_dp)) * q(3) - (work(2) / (r12**3.0_dp)) * (q(3) - q(1))
    f(8) = (-work(1) / (r2**3.0_dp)) * q(4) - (work(2) / (r12**3.0_dp)) * (q(4) - q(2))

    ! Make sure that the definition here matches
    ! your interface in the ode_solver module
end function planets_ode

!-----------------------------------------------------------------------
!! Subroutine: calculate_energy
!-----------------------------------------------------------------------
!! By: David Wilkins
!!
!! calculate_energy receives a work array containing masses, and the
!! solutions to the differential equations, and it returns an array
!! containing the energy as a function of time.
!!----------------------------------------------------------------------
!! Input:
!!
!! work     real(dp)    array containing masses
!! solution real(dp)    2D array containing solutions to the differential equations
!!----------------------------------------------------------------------
!! Output:
!!
!! energy   real(dp)    array containing energy values as a function of time
!-----------------------------------------------------------------------
subroutine calculate_energy(work, solution, energy)
    implicit none
    real(dp), intent(in) :: work(:), solution(:,:)
    real(dp), allocatable, intent(out) :: energy(:)
    real(dp) :: r1, r2, r12
    integer :: i

    allocate(energy(size(solution,2)))
    ! the following do loop calculates energy at each interval
    do i = 1, size(solution,2)
        ! defining distance terms:
        r1 = sqrt(solution(1,i)**2.0_dp + solution(2,i)**2.0_dp)
        r2 = sqrt(solution(3,i)**2.0_dp + solution(4,i)**2.0_dp)
        r12 = sqrt((solution(1,i)-solution(3,i))**2.0_dp + (solution(2,i)-solution(4,i))**2.0_dp)
        !calculating energy
        energy(i) = 0.5_dp*work(2)*((solution(5,i)**2.0_dp)+(solution(6,i)**2.0_dp)) &
                & + 0.5_dp*work(3)*((solution(7,i)**2.0_dp)+(solution(8,i)**2.0_dp)) &
                & - (work(2)*work(1))/r1 - (work(3)*work(1))/r2 - (work(2)*work(2))/r12

    enddo
end subroutine calculate_energy

   
end module mechanics