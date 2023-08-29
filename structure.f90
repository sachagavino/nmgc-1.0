!******************************************************************************
! MODULE: structure
!******************************************************************************
!
! DESCRIPTION: 
!> @brief Module that contains all the routines linked to the physical structure
!! (cloud and such) and its properties
!
!******************************************************************************

module structure

use iso_fortran_env
use numerical_types
use global_variables

implicit none

integer :: structure_sample !< the number of points from structure_evolution.dat, about the time evolution of the physical structure properties

real(double_precision), allocatable, dimension(:) :: structure_time !< dim(structure_sample) time [s] read from structure_evolution.dat
real(double_precision) :: structure_sample_step !< time [s] between each sample point for the structure evolution

real(double_precision), allocatable, dimension(:) :: structure_log_Av !< dim(structure_sample) log10(Av) [log10(mag)] read from structure_evolution.dat
real(double_precision), allocatable, dimension(:) :: structure_log_density !< dim(structure_sample) log10(density) [log10(part/cm^3)] read from structure_evolution.dat
real(double_precision), allocatable, dimension(:) :: structure_log_gas_temperature !< dim(structure_sample) log10(gas temperature) [log10(K)] read from structure_evolution.dat

logical :: read_dust !< If dust temperature must be read directly from the data file (or not)
real(double_precision), allocatable, dimension(:) :: structure_log_dust_temperature !< dim(structure_sample) log10(dust temperature) [log10(K)] read from structure_evolution.dat


contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Christophe Cossou
!
!> @date 2014
!
! DESCRIPTION: 
!> @brief Read time evolution for the physical structure in structure_evolution.dat
!! to initialize arrays. Thus, interpolation of physical structure properties
!! throughout the simulation be possible.
!! \n
!! structure_evolution.dat will have a two lines header. Then columns will be as follow :
!! \n Time (Myr) ; Number density (part/cm3) ; Temperature (K) ; Av (mag)
!
!> @warning Time sample MUST be linearly and equally spaced
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine init_structure_evolution()

use utilities

implicit none


character(len=80) :: filename = 'structure_evolution.dat' !< name of the file in which time evolution is stored
character(len=80) :: line
character(len=1), parameter :: comment_character = '!' !< character that will indicate that the rest of the line is a comment
integer :: comment_position !< the index of the comment character on the line. if zero, there is none on the current string
integer :: error !< to store the state of a read instruction
integer :: nb_columns

integer :: i !< loop index
logical :: isDefined
!------------------------------------------------------------------------------

inquire(file=filename, exist=isDefined)
if (isDefined) then
  
  call get_linenumber(filename=filename, nb_lines=structure_sample)
  nb_columns = get_nb_columns(filename)
  
  if (nb_columns.ge.5) then
    read_dust = .true.
  else
    read_dust = .false.
  endif
  
  ! If by any chance the routine is run several times, this test is here to avoid bugs
  if (allocated(structure_time)) then
    deallocate(structure_time)
    deallocate(structure_log_density)
    deallocate(structure_log_gas_temperature)
    deallocate(structure_log_Av)
  end if
  allocate(structure_time(structure_sample))
  allocate(structure_log_density(structure_sample))
  allocate(structure_log_gas_temperature(structure_sample))
  allocate(structure_log_Av(structure_sample))
  
  
  structure_time(1:structure_sample) = 0.d0
  structure_log_density(1:structure_sample) = 0.d0
  structure_log_gas_temperature(1:structure_sample) = 0.d0
  structure_log_Av(1:structure_sample) = 0.d0
  
  if (read_dust) then
    if ((GRAIN_TEMPERATURE_TYPE.ne.'table_evolv').and.(GRAIN_TEMPERATURE_TYPE.ne.'table_1D')) then
      write(Error_unit, *) 'Error: The grain temperature column exist'
      write(Error_unit, *) 'in structure_evolution.dat or 1D_static.dat but '
      write(Error_unit, *) 'GRAIN_TEMPERATURE_TYPE = ',trim(GRAIN_TEMPERATURE_TYPE)
      call exit(10)
    endif
    
    if (allocated(structure_log_dust_temperature)) then
      deallocate(structure_log_dust_temperature)
    endif
    
    allocate(structure_log_dust_temperature(structure_sample))
    structure_log_dust_temperature(1:structure_sample) = 0.d0
  endif
  
  open(10, file=filename, status='old')
  i = 1
  do 
    read(10, '(a)', iostat=error) line
    if (error /= 0) exit
      
    ! We get only what is on the left of an eventual comment parameter
      comment_position = index(line, comment_character)
    
    ! if there are comments on the current line, we get rid of them
    if (comment_position.ne.0) then
      line = line(1:comment_position - 1)
    end if
    
    if (line.ne.'') then
      if (read_dust) then
        read(line, *) structure_time(i), structure_log_Av(i), structure_log_density(i), structure_log_gas_temperature(i), &
                      structure_log_dust_temperature(i)
      else
        read(line, *) structure_time(i), structure_log_Av(i), structure_log_density(i), structure_log_gas_temperature(i)
      endif
      i = i + 1
    endif
  enddo
  close(10)
  
  ! We convert time in years to time in seconds
  structure_time(1:structure_sample) = YEAR * structure_time(1:structure_sample)
  
  ! We get the space between all times for the structure evolution. This has some sense only if times are linearly and equally spaced.
  structure_sample_step = (structure_time(structure_sample) - structure_time(1)) / dfloat(structure_sample - 1)
  
else
  write (Error_Unit,*) 'Error: The file "structure_evolution.dat" does not exist.'
  call exit(7)
end if

return
end subroutine init_structure_evolution

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Christophe Cossou
! V. Wakelam
!
!> @date 2014
! !> @date 2015
!
! DESCRIPTION: 
!> @brief Routine to get properties of the structure (Av, temperature, density)
!! in function of the time. We must launch 'init_structure_evolution' 
!! beforehand to initialize arrays
! The interpolation from Christophe has been removed and now the code
! uses only the points provided in the input file. This is up to the 
! user to check the time dependent structure that he/she is using.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine get_structure_properties_table(time, Av, density, gas_temperature)

  implicit none

  ! Inputs
  real(double_precision), intent(in) :: time !<[in] Current time of the simulation [s]
  
  ! Outputs
  real(double_precision), intent(out) :: Av !<[out] Visual extinction [mag]
  real(double_precision), intent(out) :: gas_temperature !<[out] gas temperature [K]
  real(double_precision), intent(out) :: density !<[out] gas density [part/cm^3]
  
  ! Local
  integer :: closest_low_id ! the index of the first closest lower value of radius regarding the radius value given in parameter of the subroutine. 
  !------------------------------------------------------------------------------

    closest_low_id = minloc(abs(structure_time - time),1)
    density = 10.0d0**(structure_log_density(closest_low_id))
    Av = 10.0d0**(structure_log_Av(closest_low_id))
    gas_temperature = 10.0d0**(structure_log_gas_temperature(closest_low_id))

  return
end subroutine get_structure_properties_table

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Christophe Cossou
!
!> @date 2014
!
! DESCRIPTION: 
!> @brief Routine to get properties of the structure (Av, temperature, density)
!! In this routine, everything is constant.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine get_structure_properties_fixed(time, Av, density, gas_temperature)

  implicit none

  ! Inputs
  real(double_precision), intent(in) :: time !<[in] Current time of the simulation [s]
  
  ! Outputs
  real(double_precision), intent(out) :: Av !<[out] Visual extinction [mag]
  real(double_precision), intent(out) :: gas_temperature !<[out] gas temperature [K]
  real(double_precision), intent(out) :: density !<[out] gas density [part/cm^3]

  !------------------------------------------------------------------------------

  density = initial_gas_density
  av = initial_visual_extinction 
  gas_temperature = initial_gas_temperature
  
  return
end subroutine get_structure_properties_fixed

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Christophe Cossou
!
!> @date 2014
!
! DESCRIPTION: 
!> @brief Routine to get the grain temperature, here assumed to be equal to 
!! the gas temperature
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine get_grain_temperature_gas(space,time, gas_temperature, av, grain_temperature)

  implicit none

  ! Inputs
  integer, intent(in) :: space !<[in] current spatial point in 1D.
  real(double_precision), intent(in) :: time !<[in] current time of the simulation [s]
  real(double_precision), intent(in) :: gas_temperature !<[in] gas temperature [K]
  real(double_precision), intent(in) :: av !<[in] visual extinction [mag]
  
  ! Outputs
  real(double_precision),dimension(:), intent(out) :: grain_temperature !<[out] grain temperature [K]
  !------------------------------------------------------------------------------
  
  grain_temperature = gas_temperature
  
  return
end subroutine get_grain_temperature_gas

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Christophe Cossou
!
!> @date 2014
!
! DESCRIPTION: 
!> @brief Routine to get the grain temperature, here fixed to the initial 
!! dust temperature
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine get_grain_temperature_fixed(space,time, gas_temperature, av, grain_temperature)

  implicit none

  ! Inputs
  integer, intent(in) :: space !<[in] current spatial point in 1D.
  real(double_precision), intent(in) :: time !<[in] current time of the simulation [s]
  real(double_precision), intent(in) :: gas_temperature !<[in] gas temperature [K]
  real(double_precision), intent(in) :: av !<[in] visual extinction [mag]
  
  ! Outputs
  real(double_precision),dimension(:), intent(out) :: grain_temperature !<[out] grain temperature [K]
  !------------------------------------------------------------------------------
  if(multi_grain==1) then
    write(*,*) 'Please check parameter.in:'
    write(*,*) 'if multi_grain = 1, grain_temperature_type = fixed is not valid.'
    call exit(1)
  elseif((structure_type.eq.'1D_no_diff').or.(structure_type.eq.'1D_diff')) then
    write(*,*) 'Please check parameter.in:'
    write(*,*) 'if structure_type is 1D, grain_temperature_type = fixed is not valid.'
    call exit(1)
  else
    grain_temperature = initial_dust_temperature
  endif
  
  return
end subroutine get_grain_temperature_fixed



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Wasim Iqbal
!
!> @date 2017
!
! DESCRIPTION: 
!> @brief Routine to get the grain temperature, here fixed to the initial 
!! dust temperature for each grain sizes
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine get_grain_temperature_fixed_to_dust_size(space,time, gas_temperature, av, grain_temperature)

  implicit none

  ! Inputs
  integer, intent(in) :: space !<[in] current spatial point in 1D.
  real(double_precision), intent(in) :: time !<[in] current time of the simulation [s]
  real(double_precision), intent(in) :: gas_temperature !<[in] gas temperature [K]
  real(double_precision), intent(in) :: av !<[in] visual extinction [mag]
  
  ! Outputs
  real(double_precision),dimension(:), intent(out) :: grain_temperature !<[out] grain temperature [K]
  !------------------------------------------------------------------------------
  integer::i
  
  ! if(is_dust_1D ==0 ) then
  !   do i=1,nb_grains
  !     grain_temperature(i) = grain_temp(i)
  !   enddo
  ! elseif(nb_grains/=1) then
  !   do i=1,nb_grains
  !     grain_temperature(i) = grain_temp_1D(i,space)
  !   enddo
  ! else
  !   write(*,*) 'Please check: if nb_grains = 1, then'
  !   write(*,*) 'in parameter.in file grain_temperature_type = fixed_to_dust_size is not valid'
  !   call exit(1)
  ! endif

  if(multi_grain==0) then
    write(*,*) 'Please check parameter.in:'
    write(*,*) 'if multi_grain = 0, grain_temperature_type = fixed_to_dust_size is not valid.'
    call exit(1)
  elseif(STRUCTURE_TYPE.eq.'0D') then
    do i=1,nb_grains
      grain_temperature(i) = grain_temp(i)
    enddo
  elseif((STRUCTURE_TYPE.eq.'1D_no_diff').or.(STRUCTURE_TYPE.eq.'1D_diff')) then   
    do i=1,nb_grains
      grain_temperature(i) = grain_temp_1D(i,space)
    enddo
  endif
  
  return
end subroutine get_grain_temperature_fixed_to_dust_size

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author
!> Valentine Wakelam
!
!> @date 2015
!
! DESCRIPTION:
!> @brief Routine to get the grain temperature from the table in 1D static case
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine get_grain_temperature_table_1D(space, time, gas_temperature, av, grain_temperature)

implicit none

! Inputs
integer, intent(in) :: space !<[in] current spatial point in 1D.
real(double_precision), intent(in) :: time !<[in] current time of the simulation [s]
real(double_precision), intent(in) :: gas_temperature !<[in] gas temperature [K]
real(double_precision), intent(in) :: av !<[in] visual extinction [mag]

! Outputs
real(double_precision),dimension(:), intent(out) :: grain_temperature !<[out] grain temperature [K]
!------------------------------------------------------------------------------


if(STRUCTURE_TYPE.eq.'0D') then
  write(*,*) 'Please check parameter.in:'
  write(*,*) 'in 0D mode grain_temperature_type = table_1D is not valid.'
  call exit(1)
endif
!if it is single grain model then we read values from 1D_static file else we read it from 1D_grain_sizes.in file
if (multi_grain == 0 ) then 
  grain_temperature = tmp_grain_temperature(space)
else
  grain_temperature(:) = grain_temp_1D(:,space)
endif
!grain_temperature = 10.D0

return
end subroutine get_grain_temperature_table_1D

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Christophe Cossou
!
!> @date 2014
!
! DESCRIPTION: 
!> @brief Routine to get the grain temperature, here interpolated from 
!! the data read in structure_evolution.dat
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine get_grain_temperature_table_evolv(space, time, gas_temperature, av, grain_temperature)

  implicit none

  ! Inputs
  integer, intent(in) :: space !<[in] current spatial point in 1D.
  real(double_precision), intent(in) :: time !<[in] current time of the simulation [s]
  real(double_precision), intent(in) :: gas_temperature !<[in] gas temperature [K]
  real(double_precision), intent(in) :: av !<[in] visual extinction [mag]
  
  ! Outputs
  real(double_precision),dimension(:), intent(out) :: grain_temperature !<[out] grain temperature [K]
  !------------------------------------------------------------------------------
  
  ! Local
  integer :: closest_low_id ! the index of the first closest lower value of radius regarding the radius value given in parameter of the subroutine.
  !------------------------------------------------------------------------------

   closest_low_id = minloc(abs(structure_time - time),1)
   grain_temperature = 10.0d0**structure_log_dust_temperature(closest_low_id)
    
  return
end subroutine get_grain_temperature_table_evolv

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Christophe Cossou
!
!> @date 2014
!
! DESCRIPTION: 
!> @brief Version of the diffusion routine for 0D structure. Thus, 
!! nothing is done here, just to point toward a "do nothing" procedure
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine structure_no_diffusion(timestep, temp_abundances)
  
  implicit none
  ! Inputs
  real(double_precision), intent(in) :: timestep !<[in] timestep for the diffusion process [s]
  
  ! Inputs/Outputs
  real(double_precision), dimension(:,:), intent(inout) :: temp_abundances !<[in,out] dim(nb_species, spatial_resolution) 
  !! The abundances for all species, and 
  !! all 1D mesh points (relative to H) [number ratio]
  
  ! Abundances stay equal here
  
end subroutine structure_no_diffusion

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Franck Hersant & Christophe Cossou
!
!> @date 2014
!
! DESCRIPTION: 
!> @brief 1D diffusion using a Crank Nicholson Scheme. Routine inspired by
!! Numerical Recipes. In particular the tridiag routine to solve tridiagonal
!! matrices
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine crank_nicholson_1D(f,ny,dt,dy,nu,rho, ibc)
! A Crank-Nicholson scheme
! ibc is a flag for boundary conditions
! ibc = 0 -> no flux boundaries (bc is not used then)
! ibc = 1 -> user supplied boundary conditions (bc)
! ibc > 1 -> stops the code, insulting the user
! Computes the evolution for a single timestep of the equation
! df/dt=a d/dy b d/dy c f + S
! On an homogeneous mesh, for user specified a,b,c
! The discretized equation takes the form
! xf(i+1)+yf(i)+zf(i-1)=W
!
! Tu ne dois pas avoir de a, b et c dans ta version. Si je te dis pas de
! betise, avec les notations de ce header :
! a = 1/rho
! b = nu * rho
! c = 1.
! S = 0.
implicit none

! Inputs
integer, intent(in) :: ny
real(double_precision), intent(in), dimension(ny+1) :: rho
real(double_precision), intent(in), dimension(ny+1) :: nu
real(double_precision), intent(in) :: dt
real(double_precision), intent(in) :: dy
integer, intent(in) :: ibc

! Outputs
real(double_precision), intent(out), dimension(ny+1) :: f

! Locals
integer :: ind
real(double_precision), dimension(ny+1) :: s,Q,W,x,y,z,u,v, dd1d
real(double_precision) :: d !, nu


dd1d(1:ny+1) = nu(1:ny+1) * rho(1:ny+1)

d = dt / dy**2
Q(1:ny+1) = rho(1:ny+1)

s(1:ny+1) = 0.d0

do ind = 2, ny
  W(ind) = s(ind) * dt + d / 4 * (dd1d(ind+1) + dd1d(ind)) * f(ind+1) / rho(ind) + (Q(ind) - d / 4 * (dd1d(ind+1) + 2 * dd1d(ind) &
  + dd1d(ind-1))) * f(ind) / rho(ind) + d / 4 * (dd1d(ind) + dd1d(ind-1)) * f(ind-1) / rho(ind)
enddo

do ind = 2,ny
  x(ind) = - d / 4 * (dd1d(ind+1) + dd1d(ind)) / rho(ind)
  y(ind) = Q(ind) / rho(ind) + d / 4 * (dd1d(ind+1) + 2 * dd1d(ind) + dd1d(ind-1)) / rho(ind)
  z(ind) = - d / 4 * (dd1d(ind) + dd1d(ind-1)) / rho(ind)
enddo

! Test
u(ny+1) = 1.d0
v(ny+1) = 0.d0
x(ny+1) = -d / 2 * dd1d(ny+1) / rho(ny+1)
y(ny+1) = Q(ny+1) / rho(ny+1) + d / 4 * (3 * dd1d(ny+1) + dd1d(ny)) / rho(ny+1)
z(ny+1) = -d / 4 * (dd1d(ny+1) + dd1d(ny)) / rho(ny+1)
W(ny+1) = d / 2 * dd1d(ny+1) * f(ny+1) / rho(ny+1) + (Q(ny+1) - d / 4 * (3 * dd1d(ny+1) &
+ dd1d(ny))) * f(ny+1) / rho(ny+1) + d / 4 * (dd1d(ny+1) + dd1d(ny)) * f(ny) / rho(ny+1)

do ind = ny+1, 2, -1
  u(ind-1) = -z(ind) / (x(ind) * u(ind) + y(ind))
  v(ind-1) = (W(ind) - x(ind) * v(ind)) / (x(ind) * u(ind) + y(ind))
enddo

if (ibc.eq.0) then
  f(1) =  v(1) / (1. - u(1))
endif

do ind = 1, ny
  f(ind+1) = u(ind) * f(ind) + v(ind)
enddo

return
end subroutine crank_nicholson_1D

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Christophe Cossou
!
!> @date 2014
!
! DESCRIPTION: 
!> @brief Version of the diffusion routine for 1D structure of a disk object. 
!! The diffusion is done vertically, on the 'z' dimension. The first point is on the outside 
!! TODO (check that once the procedure is written)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine structure_diffusion_1D (timestep, temp_abundances)
  
  implicit none
  ! Inputs
  real(double_precision), intent(in) :: timestep !<[in] timestep for the diffusion process [s]
  
  ! Inputs/Outputs
  real(double_precision), dimension(:,:), intent(inout) :: temp_abundances !<[in,out] dim(nb_species, spatial_resolution) 
  !! The abundances for all species, and 
  !! all 1D mesh points (relative to H) [number ratio]
  
  ! Locals
  integer :: reaction !< For loops
  
  do reaction=1,nb_species
    call crank_nicholson_1D(f=temp_abundances(reaction, 1:spatial_resolution), &
    ny=spatial_resolution-1, dt=timestep, dy=grid_cell_size, nu=diffusion_coefficient(1:spatial_resolution),&
    rho=H_number_density(1:spatial_resolution), ibc=0)
  enddo  
  
end subroutine structure_diffusion_1D

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Christophe Cossou
!
!> @date 2014
!
! DESCRIPTION: 
!> @brief Version of the timestep routine for 0D structure.
!! Nothing is done here, there will be only one timestep, equal to the output one. 
!! Indeed, subtimestep are defined by the diffusion process, and there's none in 0D
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine get_timestep_0D(current_time, final_time, next_timestep)

  implicit none

  ! Inputs
  real(double_precision), intent(in) :: current_time !<[in] current time [s]
  real(double_precision), intent(in) :: final_time !<[in] Final output time of the current 
  !! loop occurence. The last sub-step must lead exactly to this time [s]

  ! Outputs
  real(double_precision), intent(out) :: next_timestep !<[out] The next integration sub timestep withing an output integration step [s]

  next_timestep = final_time - current_time

end subroutine get_timestep_0D

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Christophe Cossou
!
!> @date 2014
!
! DESCRIPTION: 
!> @brief Version of the timestep routine for 1D structure of a disk object 
!! whose 1D diffusion is on the vertical dimension. 
!! This routine provide a sub-timestep that allow accurate computation of the diffusion.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine get_timestep_1D_diff(current_time, final_time, next_timestep)

  implicit none

  ! Inputs
  real(double_precision), intent(in) :: current_time !<[in] current time [s]
  real(double_precision), intent(in) :: final_time !<[in] Final output time of the current 
  !! loop occurence. The last sub-step must lead exactly to this time [s]

  ! Outputs
  real(double_precision), intent(out) :: next_timestep !<[out] The next integration sub timestep withing an output integration step [s]

  next_timestep = grid_cell_size**2 / maxval(diffusion_coefficient)

  ! TODO write the calculation of the diffusion timestep before that test
  if (current_time+next_timestep.gt.final_time) then
    next_timestep = final_time - current_time
  endif
end subroutine get_timestep_1D_diff

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Christophe Cossou
!
!> @date 2014
!
! DESCRIPTION: 
!> @brief Read 1D structure in 1D_static.dat
!! to initialize arrays. Thus, interpolation of physical structure properties
!! throughout the simulation be possible.
!! \n
!! 1D_static.dat will have a two lines header. Then columns will be as follow :
!! \n Z position (AU) ; Number density (part/cm3) ; Temperature (K) ; diffusion coeff (cm^2/s) ; Av (mag)
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine init_1D_static()

use utilities

implicit none


character(len=80) :: filename = '1D_static.dat' !< name of the file in which time evolution is stored
character(len=200) :: line
character(len=1), parameter :: comment_character = '!' !< character that will indicate that the rest of the line is a comment
integer :: error !< to store the state of a read instruction

real(double_precision), dimension(:), allocatable :: tmp_grid !< Distance [cm] /!\ But read column is in AU
real(double_precision), dimension(:), allocatable :: tmp_density !< Gas density [g/cm^3]
real(double_precision), dimension(:), allocatable :: tmp_gas_temperature !< Gas temperature [K]
!real(double_precision), dimension(:), allocatable :: tmp_grain_temperature !< Grain temperature [K]
real(double_precision), dimension(:), allocatable :: tmp_av !< Visual extinction [mag]
real(double_precision), dimension(:), allocatable :: tmp_kappa !< Diffusion coefficient [cm^2/s]
real(double_precision), dimension(:), allocatable :: tmp_GTODN_1D !< Abundance of dust [no unit]
real(double_precision), dimension(:), allocatable :: tmp_grain_radius_1D !< Abundance of dust [no unit]
real(double_precision), dimension(:), allocatable :: tmp_AV_NH_1D  !< AV to NH conversion factor [??]
real(double_precision), dimension(:), allocatable :: tmp_UV_flux_1D  !< UV flux in unit of the standard flux.

integer :: closest_low_id, nb_values
real(double_precision) :: x1, x2, y1, y2

integer :: i !< loop index
logical :: isDefined
!------------------------------------------------------------------------------

inquire(file=filename, exist=isDefined)
if (isDefined) then
  
  ! We get the total lines of the file
  open(10, file=filename, status='old')
  i = 0
  do
    read(10, '(a)', iostat=error) line
    if (error /= 0) exit
    
    if (line(1:1).ne.comment_character) then
      i = i + 1
    end if
  end do
  close(10)
  
  ! We define the sizes of the arrays
  nb_values = i
  if (allocated(tmp_grid)) then
    deallocate(tmp_grid)
    deallocate(tmp_density)
    deallocate(tmp_gas_temperature)
    deallocate(tmp_grain_temperature)
    deallocate(tmp_av)
    deallocate(tmp_kappa)
    deallocate(tmp_GTODN_1D)
    deallocate(tmp_grain_radius_1D)
    deallocate(tmp_AV_NH_1D)
    deallocate(tmp_UV_flux_1D)
  end if
  allocate(tmp_grid(nb_values))
  allocate(tmp_density(nb_values))
  allocate(tmp_gas_temperature(nb_values))
  allocate(tmp_grain_temperature(nb_values))
  allocate(tmp_av(nb_values))
  allocate(tmp_kappa(nb_values))
  allocate(tmp_GTODN_1D(nb_values))
  allocate(tmp_grain_radius_1D(nb_values))
  allocate(tmp_AV_NH_1D(nb_values))
  allocate(tmp_UV_flux_1D(nb_values))

  if (allocated(grid_sample)) then
    deallocate(grid_sample)
    deallocate(H_number_density)
    deallocate(gas_temperature)
    deallocate(visual_extinction)
    deallocate(diffusion_coefficient)
    deallocate(GTODN_1D)
    deallocate(grain_radius_1D)
    deallocate(AV_NH_1D)
    deallocate(UV_flux_1D)
    deallocate(NH_z)
    deallocate(NH2_z)
    deallocate(NN2_z)
    deallocate(NCO_z)
    deallocate(NH2O_z)
    deallocate(NCO2_z)
    deallocate(NN2O_z)
    deallocate(NCH_z)
    deallocate(NCH3_z)
    deallocate(NCH4_z)
    deallocate(NOH_z)
    deallocate(NHCO_z)
    deallocate(NH2CO_z)
    deallocate(NCN_z)
    deallocate(NHCN_z)
    deallocate(NHNC_z)
    deallocate(NNH_z)
    deallocate(NNH2_z)
    deallocate(NNH3_z)
  end if
  allocate(grid_sample(nb_values))
  allocate(H_number_density(nb_values))
  allocate(gas_temperature(nb_values))
  allocate(visual_extinction(nb_values))
  allocate(diffusion_coefficient(nb_values))
  allocate(GTODN_1D(nb_grains,nb_values))
  allocate(grain_radius_1D(1:nb_grains,1:nb_values))
  allocate(AV_NH_1D(nb_values))
  allocate(UV_flux_1D(nb_values))
  allocate(NH_z(nb_values))
  allocate(NH2_z(nb_values))
  allocate(NN2_z(nb_values))
  allocate(NCO_z(nb_values))
  allocate(NH2O_z(nb_values))
  allocate(NCO2_z(nb_values))
  allocate(NN2O_z(nb_values))
  allocate(NCH_z(nb_values))
  allocate(NCH3_z(nb_values))
  allocate(NCH4_z(nb_values))
  allocate(NOH_z(nb_values))
  allocate(NHCO_z(nb_values))
  allocate(NH2CO_z(nb_values))
  allocate(NCN_z(nb_values))
  allocate(NHCN_z(nb_values))
  allocate(NHNC_z(nb_values))
  allocate(NNH_z(nb_values))
  allocate(NNH2_z(nb_values))
  allocate(NNH3_z(nb_values))

  spatial_resolution = nb_values
  NH_z(1:nb_values) = 0.d0
  NH2_z(1:nb_values) = 0.d0
  NN2_z(1:nb_values) = 0.d0
  NCO_z(1:nb_values) = 0.d0
  NH2O_z(1:nb_values) = 0.d0
  NCO2_z(1:nb_values) = 0.d0
  NN2O_z(1:nb_values) = 0.d0
  NCH_z(1:nb_values) = 0.d0
  NCH3_z(1:nb_values) = 0.d0
  NCH4_z(1:nb_values) = 0.d0
  NOH_z(1:nb_values) = 0.d0
  NHCO_z(1:nb_values) = 0.d0
  NH2CO_z(1:nb_values) = 0.d0
  NCN_z(1:nb_values) = 0.d0
  NHCN_z(1:nb_values) = 0.d0
  NHNC_z(1:nb_values) = 0.d0
  NNH_z(1:nb_values) = 0.d0
  NNH2_z(1:nb_values) = 0.d0
  NNH3_z(1:nb_values) = 0.d0


  ! We get the values of the torque profile in the file
  open(10, file=filename, status='old')
  i = 0
  do
    read(10, '(a)', iostat=error) line
    if (error /= 0) exit
    
    if(line(1:1).ne.comment_character) then
      i = i + 1
      read(line, *, iostat=error) tmp_grid(i), tmp_density(i), tmp_gas_temperature(i), tmp_av(i),&
!           tmp_kappa(i), tmp_grain_temperature(i),tmp_GTODN_1D(i),tmp_grain_radius_1D(i)
           tmp_kappa(i), tmp_grain_temperature(i),tmp_GTODN_1D(i),tmp_AV_NH_1D(i),tmp_grain_radius_1D(i),&
           tmp_UV_flux_1D(i)
    end if
  end do
  
  ! Convert distances from AU to cm
  tmp_grid(1:nb_values) = tmp_grid(1:nb_values) * AU

  grid_sample(1:nb_values) = tmp_grid(1:nb_values)
  H_number_density(1:nb_values) = tmp_density(1:nb_values)
  gas_temperature(1:nb_values) = tmp_gas_temperature(1:nb_values)
  visual_extinction(1:nb_values) = tmp_av(1:nb_values)
  diffusion_coefficient(1:nb_values) = tmp_kappa(1:nb_values)
  GTODN_1D(1,1:nb_values)=tmp_GTODN_1D(1:nb_values)
  grain_radius_1D(1,1:nb_values)=tmp_grain_radius_1D(1:nb_values)
  AV_NH_1D(1:nb_values)=tmp_AV_NH_1D(1:nb_values)
  UV_flux_1D(1:nb_values)=tmp_UV_flux_1D(1:nb_values)

  
  !  we reassign values read from 1D_grain_sizes.in if nb_grains more than 1, that is in case of multi grain model  
  if (nb_grains /=1 ) then  
    GTODN_1D(1:nb_grains,1:nb_values)=GTODN_1D_temp(1:nb_grains,1:nb_values)
    grain_radius_1D(1:nb_grains,1:nb_values)=grain_radii_1D(1:nb_grains,1:nb_values)
  endif  
  

else
  write (Error_Unit,*) 'Error: The file "',trim(filename),'" does not exist.'
  call exit(25)
end if

return
end subroutine init_1D_static

end module structure
