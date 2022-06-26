!-----------------------------------------------------------------------
!Module: read_write
!-----------------------------------------------------------------------
!! By: David Wilkins
!!
!! The subroutines in this module are responsible for reading input from
!! namelist files and writing output to output files.
!!----------------------------------------------------------------------
!! Included subroutines:
!!
!! read_input
!! write_results
!-----------------------------------------------------------------------
module read_write
use types
implicit none

private
public :: read_input, write_results

contains

!-----------------------------------------------------------------------
!! Subroutine: read_input
!-----------------------------------------------------------------------
!! By: David Wilkins
!!
!! read_input uses namelist files in order to set up the initial postions
!! and velocities of the planets as well as the masses, the initial and
!! final times, and the number of steps.
!!----------------------------------------------------------------------
!! Output:
!!
!! work_array           real(dp)    array containing masses
!! initial_condition    real(dp)    array containing initial positions and velocities
!! initial_time         real(dp)    start time
!! final_time           real(dp)    end time
!! n_steps              integer     number of lattice points for calculation
!! output_file          char        name of output file
!-----------------------------------------------------------------------
subroutine read_input(work_array, initial_condition, initial_time, final_time, n_steps, output_file)
    implicit none
    real(dp), intent(out) :: work_array(1:3)
    real(dp), intent(out) :: initial_condition(1:8)
    real(dp), intent(out) :: initial_time, final_time
    integer, intent(out) :: n_steps
    character(len=*) output_file
    real(dp) :: primary_mass, planet_mass_1, planet_mass_2
    real(dp) :: initial_pos_1(1:2), initial_pos_2(1:2)
    real(dp) :: initial_vel_1(1:2), initial_vel_2(1:2)
    integer :: n_arguments, file_unit, ierror
    character(len=200) :: namelist_file
    logical :: file_exists


    namelist /masses/ primary_mass, planet_mass_1, planet_mass_2
    namelist /initial_conditions/ initial_pos_1, initial_pos_2, initial_vel_1, initial_vel_2
    namelist /solution_parameters/ final_time, n_steps
    namelist /output/ output_file


    ! Set default values
    primary_mass  = 0.0_dp
    planet_mass_1 = 1.0_dp
    planet_mass_2 = 0.0_dp
    initial_pos_1 = (/1.0_dp, 0.0_dp/)
    initial_pos_2 = (/0.0_dp, 0.0_dp/)
    initial_vel_1 = (/0.0_dp, 1.0_dp/)
    initial_vel_2 = (/0.0_dp, 0.0_dp/)
    initial_time = 0.0_dp
    final_time = 10.0_dp
    n_steps = 100
    output_file = 'default.dat'


    ! get namelist file name from command line and read namelists

    n_arguments = command_argument_count()
 
    if (n_arguments == 1) then
        call get_command_argument(1, namelist_file)
        inquire(file = trim(namelist_file), exist = file_exists)
        if (file_exists) then
            open(newunit=file_unit, file=namelist_file)
            read(file_unit, nml=masses, iostat=ierror)
            if(ierror /= 0) then
                print*, "Error reading masses namelist"
                stop
            endif
            read(file_unit, nml=initial_conditions, iostat=ierror)
            if(ierror /= 0) then
                print*, "Error reading initial_conditions namelist"
                stop
            endif
            read(file_unit, nml=solution_parameters, iostat=ierror)
            if(ierror /= 0) then
                print*, "Error reading solution_parameters namelist"
                stop
            endif
            read(file_unit, nml=output, iostat=ierror)
            if(ierror /= 0) then
                print*, "Error reading output namelist"
                stop
            endif
            close(file_unit)
        else
            print*, namelist_file, 'not found'
            stop
        endif
    elseif (n_arguments /= 0) then
        print*, 'Incorrect number of arguments. Program takes either 0 or 1 argument only'
        print*, 'See details in README.md'
        stop
    endif

    ! You can follow what we did in class during the namelist example
    ! The code is in the class repository

    work_array = [primary_mass, planet_mass_1, planet_mass_2]
    initial_condition = [initial_pos_1, initial_pos_2, initial_vel_1, initial_vel_2]

end subroutine read_input

!-----------------------------------------------------------------------
!! Subroutine: write_results
!-----------------------------------------------------------------------
!! By: David Wilkins
!!
!! write_results receives the solution and energy arrays and writes them
!! to an output file specified by the namelist file as a function of time
!!----------------------------------------------------------------------
!! Input:
!!
!! time         real(dp)    array containing time values
!! solution     real(dp)    2D array containing solutions to the differential equations
!! energy       real(dp)    array containing energy values as a function of time
!! output_file  real(dp)    name of output file to be written to
!-----------------------------------------------------------------------
subroutine write_results(time, solution, energy, output_file)
    implicit none
    real(dp), intent(in) :: time(:), solution(:,:), energy(:)
    character(len=*), intent(in) :: output_file
    integer :: unit, i

    open(newunit=unit,file=output_file)
    write(unit,'(6a28)') 't', 'x_1', 'y_1', 'x_2', 'y_2', 'E'

    do i = 1, size(time)
        ! write the t value and corresponding coordinates and energy
        write(unit,'(6g28.16)') time(i), solution(1,i), solution(2,i), &
                              & solution(3,i), solution(4,i), energy(i)
    enddo
    
    close(unit)


end subroutine write_results


end module read_write
