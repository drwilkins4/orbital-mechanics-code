! Program: planets
!
! By: David Wilkins
!
! This program solves a system of eight differential equations using the fourth
! order Runge-Kutta method in order to compute the orbits of planets around a
! primary mass and the energy of the system as a function of time. The initial
! positions and velocities of two planets are read into the program from
! namelist files in order to split the calculations into two cases: one in which
! two distant planets have circular orbits around the primary mass, and one in
! which the smaller planet orbits the larger planet, which orbits the primary
! mass. The solutions to the differential equations and the energies of the
! systems are written to output files as a function of time in order to be plotted
! and analyzed.
!-----------------------------------------------------------------------------
program planets

use types
use read_write, only : read_input, write_results
use ode_solver, only : solve_runge_kutta_4
use mechanics, only : calculate_energy, planets_ode

implicit none

real(dp) :: work_array(1:3), initial_condition(1:8)
real(dp) :: initial_time, final_time
integer :: n_steps
real(dp), allocatable :: time(:), solution(:,:), energy(:)
character(len=1024) :: output_file

call read_input(work_array, initial_condition, initial_time, final_time, n_steps, output_file)

call solve_runge_kutta_4(planets_ode, work_array, initial_condition, initial_time, final_time, n_steps, time, solution)
call calculate_energy(work_array, solution, energy)

call write_results(time, solution, energy, output_file)


end program planets