program nautilus_rates

use numerical_types
use iso_fortran_env
use utilities
use nautilus_main

implicit none

! Locals
character(len=80) :: filename_output
integer :: species, output, reaction, i, k ! index for loops
integer :: nb_points !< number of spatial points the user wants rates at.

logical :: isDefined
real(double_precision) :: tmp !< temporary variable

character(len=80) :: rate_format, time_format !< string to store specific format used to output datas

! /!\ Variable names with _out are variable that already exist in the nautilus code, but here they are arrays, one value per output.

real(double_precision), dimension(:,:,:), allocatable :: abundances_out !< abundances_out over time for each species. (spatial_resolution,nb_outputs, nb_species)

real(double_precision), dimension(:,:), allocatable :: densities_spec_out !< species densities over time for each species. (nb_outputs, nb_species)

real(double_precision), dimension(:), allocatable :: time
real(double_precision), dimension(:,:), allocatable :: gas_temperature_out
real(double_precision), dimension(:,:), allocatable :: dust_temperature_out
real(double_precision), dimension(:,:), allocatable :: density
real(double_precision), dimension(:,:), allocatable :: visual_extinction_out !< visual extinction
real(double_precision), dimension(:), allocatable :: zeta

! For rates
real(double_precision), allocatable, dimension(:,:,:) :: reaction_rates_out ! (spatial_resolution,nb_outputs, nb_reactions)

! Output of the code
real(double_precision), allocatable, dimension(:,:) :: reaction_fluxes

! User asked values
logical :: change_space = .true. !< If true, ask the user for a value
logical :: wrong_1D !< Flags for while loops when asking the user something
integer :: user_1D_id !< designed spatial id by the user
integer, dimension(:), allocatable :: user_1D_id_list !< designed spatial id_s by the user
character(len=5) :: str_user_id

! Initialise all variables from the global_variables module. Only some of them are used here.
call initialisation()

! We calculate the total number of outputs by checking for each file if it exist or not.
nb_outputs = 0
isDefined = .true.
do while(isDefined)
  nb_outputs = nb_outputs + 1
  write(filename_output, '(a,i0.6,a)') 'abundances.',nb_outputs,'.out'
  inquire(file=filename_output, exist=isDefined)

enddo
nb_outputs = nb_outputs - 1

write(*,'(a, i0)') 'Spatial resolution: ', spatial_resolution
write(*,'(a, i0)') 'Number of time outputs: ', nb_outputs
write(*,'(a, i0)') 'Number of species: ', nb_species
write(*,'(a, i0)') 'Number of reactions: ', nb_reactions

! We allocate the output arrays

allocate(time(nb_outputs))
allocate(gas_temperature_out(spatial_resolution, nb_outputs))
allocate(dust_temperature_out(spatial_resolution, nb_outputs))
allocate(density(spatial_resolution, nb_outputs))
allocate(visual_extinction_out(spatial_resolution, nb_outputs))
allocate(zeta(nb_outputs))

allocate(abundances_out(nb_outputs, nb_species+1, spatial_resolution)) ! We create an extra species that will always have an abundance of 1

allocate(densities_spec_out(nb_outputs, nb_species+1)) ! We create an extra species that will always have an abundance of 1

allocate(reaction_rates_out(spatial_resolution,nb_outputs, nb_reactions))

! Outputs
allocate(reaction_fluxes(nb_outputs, nb_reactions))

! What spatial point ?
40 if (change_space) then
  if (spatial_resolution.ne.1) then
      wrong_1D = .true.
      write(*,'(a,i0,a)') 'select a number of points needed (from 1 to ', spatial_resolution, '):'
      read(*,*) nb_points
      allocate(user_1D_id_list(nb_points))
      do while (wrong_1D)
          write(*,'(a,i0,a)') 'Select ',nb_points,' spatial points (separated by space):'
          read(*,*) user_1D_id_list

          
          if (ANY(user_1D_id_list.gt.spatial_resolution).or.ANY(user_1D_id_list.lt.0)) then
              write(*,'(a,i0)') "Choose values between 1 and ", spatial_resolution
          else
              wrong_1D = .false.
          endif
      enddo
  else
  ! Default value if only one point
      write(*,*) "We are in 0D, skipping spatial point choosing"
      allocate(user_1D_id_list(1))
      user_1D_id_list = 1
  endif
  change_space = .false.
endif


  ! The next write will be written in the same line
  write(*,'(a)', advance='no') 'Reading unformatted outputs...'
  ! We read output files
  do output=1,nb_outputs
    write(filename_output, '(a,i0.6,a)') 'abundances.',output,'.out'

    open(10, file=filename_output, status='old', form='unformatted')
    read(10) time(output)
    read(10) gas_temperature_out(1:spatial_resolution, output), dust_temperature_out(1:spatial_resolution, output), &
             density(1:spatial_resolution, output), visual_extinction_out(1:spatial_resolution, output), zeta(output)
    read(10) abundances_out(output,1:nb_species, 1:spatial_resolution)
    close(10)
  
    write(filename_output, '(a,i0.6,a)') 'rates.',output,'.out'

    open(10, file=filename_output, status='old', form='unformatted')
    read(10) reaction_rates_out(1:spatial_resolution,output,1:nb_reactions)
    close(10)

  enddo
  ! achar(13) is carriage return '\r'. Allow to go back to the beginning of the line
  write(*,'(a,a)') achar(13), 'Reading unformatted outputs... Done'


do k=1,nb_points
  ! transforming abundances to densities
  do output=1,nb_outputs
  do i=1,nb_species
  densities_spec_out(output, i) = abundances_out(output, i,user_1D_id_list(k))*density(user_1D_id_list(k),output)
  enddo
  enddo

  ! For non existing reactants (whose index is 'nb_species+1') in reactions, we create a new species whose abundance is always 1, so that we can calculate the fluxes
  !! more easily.
  densities_spec_out(1:nb_outputs, nb_species+1) = 1.d0

  ! We replace blanck species by 'XXX' for the outputs constrains
  do reaction=1,nb_reactions
    do species=1,MAX_COMPOUNDS
      if (REACTION_COMPOUNDS_NAMES(species, reaction).eq.'   ') then
        REACTION_COMPOUNDS_NAMES(species, reaction) = 'XXX'
      endif
    enddo
  enddo

  ! We write fluxes for all reactions and all output times
  do output=1, nb_outputs
    do reaction=1, nb_reactions
      tmp = 1.d0
      do i=1, MAX_REACTANTS
        tmp = tmp * densities_spec_out(output, REACTION_COMPOUNDS_ID(i, reaction))
      enddo
      reaction_fluxes(output, reaction) = reaction_rates_out(user_1D_id_list(k),output, reaction) * tmp
    enddo
  enddo

  write(str_user_id, '(I5)' ) user_1D_id_list(k) !convert int into string to use it for variable names.
  ! The next write will be written in the same line
  write(*,'(a)', advance='no') 'Writing rates ASCII files ',str_user_id,'...'

  write(rate_format, '(a,i0,a,i0,a)') '(',MAX_COMPOUNDS,'(a11," "),', nb_outputs, '(es12.4e2)," ",i7," ",i3," ")'
  write(time_format, '(a,i0,a,i0,a)') '(',12*(MAX_COMPOUNDS-1),'(" "),a12,', nb_outputs, '(es12.4e2))'

  ! We write ASCII output file
  
  open(10, file='rates' // trim(adjustl(str_user_id)) // '.out')
  ! all MAX_COMPOUNDS species involved ('XXX' if no species) ; Each column is the flux for several output times'
  ! The first line list time for each column of flux
  write(10, time_format) 'Time (year)', time(1:nb_outputs)*3.17e-8
  do reaction=1, nb_reactions
    write(10,rate_format) REACTION_COMPOUNDS_NAMES(1:MAX_COMPOUNDS, reaction), reaction_fluxes(1:nb_outputs, reaction), &
                          REACTION_ID(reaction), REACTION_TYPE(reaction)
  enddo
  close(10)

  ! We write ASCII output file
  open(10, file='rate_coefficients' // trim(adjustl(str_user_id)) // '.out')
  ! all MAX_COMPOUNDS species involved ('XXX' if no species) ; Each column is the flux for several output times'
  ! The first line list time for each column of flux
  write(10, time_format) 'Time (year)', time(1:nb_outputs)*3.17e-8
  do reaction=1, nb_reactions
  write(10,rate_format) REACTION_COMPOUNDS_NAMES(1:MAX_COMPOUNDS, reaction), &
  reaction_rates_out(user_1D_id_list(k),1:nb_outputs, reaction), REACTION_ID(reaction), REACTION_TYPE(reaction)
  enddo
  close(10)

end do

! achar(13) is carriage return '\r'. Allow to go back to the beginning of the line
write(*,'(a,a)') achar(13), '... Done'

end program nautilus_rates
