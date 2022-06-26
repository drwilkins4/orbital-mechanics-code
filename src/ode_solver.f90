!-----------------------------------------------------------------------
!Module: ode_solver
!-----------------------------------------------------------------------
!! By: David Wilkins
!!
!! The subroutine in this module is responsible for solving an arbitrary
!! set of differential equations.
!!----------------------------------------------------------------------
!! Included subroutines:
!!
!! solve_runge_kutta_4
!-----------------------------------------------------------------------
module ode_solver
use types
implicit none
private

public :: solve_runge_kutta_4

interface
    function func(q, t, work) result(f)
        use types, only : dp
        implicit none
        real(dp), intent(in) :: q(:), t, work(:)
        real(dp), allocatable :: f(:)
    end function func
end interface

contains

!-----------------------------------------------------------------------
!! Subroutine: solve_runge_kutta_4
!-----------------------------------------------------------------------
!! By: David Wilkins
!!
!! solve_runge_kutta_4 receives a work array containing the masses and an
!! array containing the initial conditions and it uses the fourth order
!! Runge-Kutta method to solve the system of differential equations by
!! calling planets_ode, which is passed as an argument. The solutions to
!! the system of differential equations are returned in a 2D array.
!!----------------------------------------------------------------------
!! Input:
!!
!! f            procedure   planets_ode passed as an argument
!! work_array   real(dp)    array containing masses
!! r_i          real(dp)    array containing initial conditions
!! t_i          real(dp)    initial value of independent variable
!! t_f          real(dp)    final value of independent variable
!! n_steps      integer     number of lattice points for calculation
!!----------------------------------------------------------------------
!! Output:
!!
!! time         real(dp)    array containing values of independent variable
!! solution     real(dp)    2D array containing solutions
!-----------------------------------------------------------------------
subroutine solve_runge_kutta_4(f, work_array, r_i, t_i, t_f, n_steps, time, solution)
    ! You can use the class example as a starting
    ! point for your fourth order Runge Kutta
    implicit none
    procedure(func) :: f
    real(dp), intent(in) :: t_i, r_i(:), t_f, work_array(:)
    integer, intent(in) :: n_steps
    real(dp), intent(out), allocatable :: time(:), solution(:, :)
    integer :: n_variables, i
    real(dp) :: h, t
    real(dp), allocatable :: k1(:), k2(:), k3(:), k4(:), q(:)

    n_variables = size(r_i)

    allocate(time(n_steps), solution(n_variables, n_steps))
    allocate(k1(n_variables), k2(n_variables), k3(n_variables), k4(n_variables), q(n_variables))

    h = (t_f - t_i)/(real(n_steps, kind=dp))
    q = r_i
    t = t_i
    do i=1,n_steps
        k1 = h*f(q, t, work_array)
        k2 = h*f(q + 0.5_dp*k1, t + 0.5_dp*h, work_array)
        k3 = h*f(q + 0.5_dp*k2, t + 0.5_dp*h, work_array)
        k4 = h*f(q + k3, t + h, work_array)
        q = q + (1.0_dp/6.0_dp) * (k1 + 2.0_dp*k2 + 2.0_dp*k3 + k4)
        t = t + h
        time(i) = t
        solution(:, i) = q
    enddo
end subroutine solve_runge_kutta_4


end module ode_solver