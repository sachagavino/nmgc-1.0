program nautilus_major_reactions

use numerical_types
use iso_fortran_env
use utilities
use nautilus_main

implicit none

! Parameters
real(double_precision) :: PERCENTAGE_THRESHOLD !< Percentage below which the reaction will not be displayed

! Locals
character(len=80) :: filename_output
integer :: output, reaction, i ! index for loops
logical :: isDefined
real(double_precision) :: tmp !< temporary variable
character(len=80) :: output_format

! /!\ Variable names with _out are variable that already exist in the nautilus code, but here they are arrays, one value per output.

real(double_precision), dimension(:,:,:), allocatable :: abundances_out !< abundances_out over time for each species. (nb_outputs, nb_species)

real(double_precision), dimension(:), allocatable :: time
real(double_precision), dimension(:,:), allocatable :: gas_temperature_out
real(double_precision), dimension(:,:), allocatable :: dust_temperature_out
real(double_precision), dimension(:,:), allocatable :: density
real(double_precision), dimension(:,:), allocatable :: visual_extinction_out !< visual extinction
real(double_precision), dimension(:), allocatable :: zeta

! For rates
real(double_precision), allocatable, dimension(:,:,:) :: reaction_rates_out ! (spatial_resolution,nb_outputs, nb_reactions)

! User asked values
logical :: change_time = .true. !< If true, ask the user for a value
logical :: change_species = .true. !< If true, ask the user for a value
logical :: change_space = .true. !< If true, ask the user for a value
logical :: print_types = .true. !< print different reaction types the first time
logical :: wrong_species, wrong_output, wrong_1D !< Flags for while loops when asking the user something
character(len=1) :: user_action !< Ask the user what he wants to do after the first run
character(len=11) :: user_species !< The species designed by the user
integer :: user_species_id !< corresponding id of the desired species of the user
integer :: output_id !< designed output id by the user
integer :: user_1D_id !< designed spatial id by the user
integer :: change_threshold !< change threshold (1 by default)

! For outputs
real(double_precision), allocatable, dimension(:) :: destructions !< production rates for one given species. Reactions where this species is not involved have 0 value
integer, allocatable, dimension(:) :: destructions_id !< corresponding ID sorted from the lest important to the most importants at the end
real(double_precision), allocatable, dimension(:) :: productions !< destruction rates for one given species. Reactions where this species is not involved have 0 value
integer, allocatable, dimension(:) :: productions_id !< corresponding ID sorted from the lest important to the most importants at the end
real(double_precision) :: destructions_sum
real(double_precision) :: productions_sum
real(double_precision) :: percentage
character(len=35) :: reaction_line ! the longest reaction so far is the line 239 of gas_reactions.in, involving CH2CHCHCH2. 
!! This short value is to align reactions correctly, but limiting the number of extra spaces

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

allocate(reaction_rates_out(spatial_resolution,nb_outputs, nb_reactions))

! Outputs
allocate(destructions(nb_reactions))
allocate(destructions_id(nb_reactions))
allocate(productions(nb_reactions))
allocate(productions_id(nb_reactions))

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

! For non existing reactants (whose index is 'nb_species+1') in reactions, we create a new species whose abundance is always 1, so that we can calculate the fluxes 
!! more easily.
abundances_out(1:nb_outputs, nb_species+1, 1:spatial_resolution) = 1.d0

!######################################################
! User asked section
!######################################################

! Read the threshold
PERCENTAGE_THRESHOLD = 1.0
write(*,*) 'Choose the percentage threshold?'
write(*,*) '   Yes = 1 : All reactions'
write(*,*) '   No  = 0 : Default'
read(*,*) change_threshold
if(change_threshold==1) PERCENTAGE_THRESHOLD = 0.0

! What time ?
20 if (change_time) then
  wrong_output = .true.
  do while (wrong_output)
    write(*,'(a,i0,a)') 'Select the output time (from 1 to ', nb_outputs, '):'
    read(*,*) output_ID
    
    
    if (output_ID.gt.nb_outputs) then
      write(*,'(a,i0)') "Choose output between 1 and ", nb_outputs
    else
      wrong_output = .false.
    endif
  enddo
  change_time = .false.
endif

! Which species ?
30 if (change_species) then
  wrong_species = .true.
  do while (wrong_species)
    write(*,*) 'Select one species:'
    read(*,*) user_species
    
    user_species_id = -1
    do i=1, nb_species
      if (user_species.eq.species_name(i)) then
      user_species_id = i
      endif
    enddo
    
    if (user_species_id.eq.-1) then
      write(*,*) "'", trim(user_species), "' doesn't exist"
    else
      wrong_species = .false.
    endif
  enddo
  change_species = .false.
endif

! What spatial point ?
40 if (change_space) then
  if (spatial_resolution.ne.1) then
    wrong_1D = .true.
    do while (wrong_1D)
      write(*,'(a,i0,a)') 'Select the spatial point (from 1 to ', spatial_resolution, '):'
      read(*,*) user_1D_id
      
      
      if ((user_1D_id.gt.spatial_resolution).or.(user_1D_id.lt.0)) then
      write(*,'(a,i0)') "Choose spatial point between 1 and ", spatial_resolution
      else
      wrong_1D = .false.
      endif
    enddo
  else
    ! Default value if only one point
    write(*,*) "We are in 0D, skipping spatial point choosing"
    user_1D_id = 1
  endif
  change_space = .false.
endif

!######################################################
! End User asked section
!######################################################

destructions(1:nb_reactions) = 0.d0
productions(1:nb_reactions) = 0.d0
destructions_sum = 0.d0
productions_sum = 0.d0

do reaction=1, nb_reactions
  if (any(REACTION_COMPOUNDS_ID(1:MAX_REACTANTS, reaction).eq.user_species_id)) then

    tmp = 1.d0
    do i=1, MAX_REACTANTS
      ! Skip if blanck species
      if (REACTION_COMPOUNDS_ID(i, reaction).ne.nb_species+1) then
        tmp = tmp * abundances_out(output_id, REACTION_COMPOUNDS_ID(i, reaction), user_1D_id) * density(user_1D_id, output_ID)
      endif
    enddo
    destructions(reaction) = reaction_rates_out(user_1D_id,output_id, reaction) * tmp
  endif
  
  if (any(REACTION_COMPOUNDS_ID(MAX_REACTANTS+1:MAX_COMPOUNDS, reaction).eq.user_species_id)) then
    tmp = 1.d0
    do i=1, MAX_REACTANTS
      ! Skip if blanck species
      if (REACTION_COMPOUNDS_ID(i, reaction).ne.nb_species+1) then
        tmp = tmp * abundances_out(output_id, REACTION_COMPOUNDS_ID(i, reaction), user_1D_id) * density(user_1D_id, output_ID)
      endif
    enddo
    productions(reaction) = reaction_rates_out(user_1D_id,output_id, reaction) * tmp
  endif
enddo
  
call get_sorted_index(productions(1:nb_reactions), productions_id(1:nb_reactions))
call get_sorted_index(destructions(1:nb_reactions), destructions_id(1:nb_reactions))

productions_sum = sum(productions(1:nb_reactions))
destructions_sum = sum(destructions(1:nb_reactions))

if (print_types) then
  write(*,*) '0    : Gas phase reactions with GRAINS'
  write(*,*) '1    : Photodissoc/ionisation with cosmic rays'
  write(*,*) '2    : Gas phase photodissoc/ionisations by secondary UV photons generated by CR'
  write(*,*) '3    : Gas phase photodissociations/ionisations by UV'
  write(*,*) '4-8  : Bimolecular gas phase reactions - several possible formula '
  write(*,*) '10-11: H2 formation on the grains when IS_GRAIN_REACTIONS eq 0'
  write(*,*) '14   : Grain surface reactions'
  write(*,*) '15   : Thermal evaporation'
  write(*,*) '16   : Cosmic-ray evaporation'
  write(*,*) '17-18: Photodissociations by Cosmic rays on grain surfaces'
  write(*,*) '19-20: Photodissociations by UV photons on grain surfaces'
  write(*,*) '21   : Grain surface reactions'
  write(*,*) '66   : Photodesorption by external UV'
  write(*,*) '67   : Photodesorption by CR generated UV'
  write(*,*) '98   : storage of H2S under a refractory form'
  write(*,*) '99   : Adsorption on grains'
endif

write(output_format,*)'(i2,4x,a,2x,es13.7,4x,f5.1,"%")'

write(*,*) ''
write(*,'(a,a,a,i0,a,es8.2,a,i0,a,es8.2,a)') 'For ', trim(user_species), ' at output n°', output_ID, &
' (', time(output_ID)/YEAR, ' years) and spatial point n°',user_1D_id,' (',grid_sample(user_1D_id)/AU,' AU)'
write(*,'(2(a,es8.2),a,es9.1e3,a)') 'Gas density = ', density(user_1D_id, output_ID), ' [part/cm^3] ; Av = ', &
visual_extinction_out(user_1D_id, output_ID), ' [mag] ; X rate = ', zeta(output), ' [s-1]'
write(*,'(2(a,es8.2),a)') 'Gas temp = ', gas_temperature_out(user_1D_id, output_ID), ' [K] ; Dust temp = ', &
dust_temperature_out(user_1D_id, output_ID), ' [K]'

write(*,*) "--------------- PRODUCTION (cm-3 s-1) -----------------    ------"
percentage = 100.d0
i = nb_reactions

! For first element
reaction = productions_id(i)
percentage = productions(reaction) / productions_sum * 100.d0

do while(percentage.gt.PERCENTAGE_THRESHOLD)

  call display_reaction(REACTION_COMPOUNDS_NAMES(1:MAX_COMPOUNDS, reaction), reaction_line)
  
  write(*, output_format) REACTION_TYPE(reaction), reaction_line, productions(reaction), percentage
  
  ! Decrement index
  i = i - 1
  
  ! Just in case, to avoid negative index for arrays
  if (i.lt.0) then
    percentage = 0.d0
  else
    ! Done at the end of the loop to avoid printing one extra reaction below the threshold
    reaction = productions_id(i)
    percentage = productions(reaction) / productions_sum * 100.d0
  endif
enddo

write(*,*) "--------------- DESTRUCTION (cm-3 s-1) ----------------    ------"
percentage = 100.d0
i = nb_reactions

! For first element
reaction = destructions_id(i)
percentage = destructions(reaction) / destructions_sum * 100.d0

do while(percentage.gt.PERCENTAGE_THRESHOLD)

  call display_reaction(REACTION_COMPOUNDS_NAMES(1:MAX_COMPOUNDS, reaction), reaction_line)
  
  write(*,output_format) REACTION_TYPE(reaction), reaction_line, destructions(reaction), percentage
  
  ! Decrement index
  i = i - 1
  
  ! Just in case, to avoid negative index for arrays
  if (i.lt.0) then
    percentage = 0.d0
  else
      ! Done at the end of the loop to avoid printing one extra reaction below the threshold
      reaction = destructions_id(i)
      percentage = destructions(reaction) / destructions_sum * 100.d0
  endif
enddo

! Ask user if he want another run
10 write(*,*) ''
write(*,*) 'change: t(ime), s(pecies), p(oint 1D), a(ll)'
write(*,*) 'l(egend) ; h(elp) ; q(uit)'
write(*,"(a)", advance='no') 'Please enter your selection now:'

! Use a C function, in getkey.c, whose object file is included.
! we use gcc -c getkey.c
! then add getkey.o when compiling the fortran main program
read(*,*) user_action

write(*,*) ''

if ((user_action.eq.'p').and.(spatial_resolution.eq.1)) then
  write(*,*) "/!\ You are in 0D !"
  goto 10
endif

select case(user_action)
  case('q')
    call exit(0) ! Exiting normally

  case('a') ! change all
    change_time = .true.
    change_species = .true.
    change_space = .true.
    goto 20

  case('t') ! change time
    change_time = .true.
    goto 20

  case('s') ! change species
    change_species = .true.
    goto 30

  case('p') ! change spatial point
    change_space = .true.
    goto 40
    
  case('l') ! Switch legend boolean
    if (print_types) then
      print_types = .false.
      write(*,*) 'Info: Reaction types legend will now be hidden'
    else
      print_types = .true.
      write(*,*) 'Info: Reaction types legend will now be shown everytime'
    endif
    goto 10
  
  case('h') ! Show help
    write(*,*) 't: change output time'
    write(*,*) 's: change species'
    write(*,*) 'p: change spatial point (only in 1D)'
    write(*,*) 'a: change all (time, species and spatial point)'
    write(*,*) 'l: permanently show/hide the reaction types legend'
    write(*,*) 'h: help. Show the present message'
    write(*,*) 'q: quit the application'
    goto 10
  
end select


contains


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Christophe Cossou
!
!> @date 2014
!
! DESCRIPTION: 
!> @brief Return a string to display a reaction
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine display_reaction(names,reaction_line)
implicit none

! Input
character(len=11), dimension(MAX_COMPOUNDS), intent(in) :: names !<[in] REACTION_COMPOUNDS_NAMES for a given reaction

! Output
character(len=*), intent(out) :: reaction_line !<[out] String that represent the reaction

! Locals
character(len=11) :: tmp_name
integer :: compound

reaction_line = trim(names(1))
do compound=2,MAX_REACTANTS
tmp_name = names(compound)
  if (tmp_name.ne.'') then
    reaction_line = trim(reaction_line)//" + "//trim(tmp_name)
  endif
enddo

reaction_line = trim(reaction_line)//" -> "//trim(names(MAX_REACTANTS+1))

do compound=MAX_REACTANTS+2,MAX_COMPOUNDS
tmp_name = names(compound)
  if (tmp_name.ne.'') then
  reaction_line = trim(reaction_line)//" + "//trim(tmp_name)
  endif
enddo

end subroutine display_reaction

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Christophe Cossou
!
!> @date 2014
!
! DESCRIPTION: 
!> @brief subroutine coming from Numerical Recipes index_sp function
!! adapted for this code. The routine return the ordered indexes corresponding
!! to the element of the input arrays, from lowest to highest. \n\n
!! For a given array [3, 1, 2, 4, 5], the function will return [2, 3, 1, 4, 5]
!! because the lowest element is the 2nd, then come the 3rd, the 1st 
!! and then 4-th and 5-th
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine get_sorted_index(arr,index)
implicit none
real(double_precision), dimension(:), intent(in) :: arr !<[in] the input array for the sorted process
integer, dimension(:), intent(out) :: index !<[out] the index of the elements, from lowest to highest corresponding values
integer, parameter :: nn=15, nstack=50
real(double_precision) :: a
integer :: n,k,i,j,indext,jstack,l,r
integer, dimension(nstack) :: istack
integer :: tmp_index

if (size(index).ne.size(arr)) then
  write(error_unit,*) 'in get_sorted_index. size are different for the two arguments'
  call exit(24)
endif

n = size(index)

! initialize list of index
do i=1,n
  index(i) = i
enddo

jstack=0
l=1
r=n
do
  if (r-l < nn) then
    do j=l+1,r
      indext=index(j)
      a=arr(indext)
      do i=j-1,l,-1
        if (arr(index(i)) <= a) exit
        index(i+1)=index(i)
      end do
      index(i+1)=indext
    end do
    if (jstack == 0) return
    r=istack(jstack)
    l=istack(jstack-1)
    jstack=jstack-2
  else
    k=(l+r)/2
    
    ! swaping indexes
    tmp_index = index(k)
    index(k) = index(l+1)
    index(l+1) = tmp_index
    
    call icomp_xchg(arr,index(l),index(r))
    call icomp_xchg(arr,index(l+1),index(r))
    call icomp_xchg(arr,index(l),index(l+1))
    i=l+1
    j=r
    indext=index(l+1)
    a=arr(indext)
    do
      do
        i=i+1
        if (arr(index(i)) >= a) exit
      end do
      do
        j=j-1
        if (arr(index(j)) <= a) exit
      end do
      if (j < i) exit
      tmp_index = index(i)
      index(i) = index(j)
      index(j) = tmp_index
    end do
    index(l+1)=index(j)
    index(j)=indext
    jstack=jstack+2
    if (jstack > nstack) then
      write(error_unit, *) 'indexx: nstack too small'
      call exit(24)
    endif
    if (r-i+1 >= j-l) then
      istack(jstack)=r
      istack(jstack-1)=i
      r=j-1
    else
      istack(jstack)=j-1
      istack(jstack-1)=l
      l=i
    end if
  end if
end do
end subroutine get_sorted_index

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> christophe cossou
!
!> @date 2014
!
! DESCRIPTION: 
!> @brief subroutine coming from numerical recipes icomp_xchg function
!! adapted for this code. the routine will swap the two indexes i and j
!! if the corresponding value of the first (i) is bigger than 
!! the value of the second (j).
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine icomp_xchg(arr,i,j)
real(double_precision), dimension(:), intent(in) :: arr !<[in] the reference array
integer, intent(inout) :: i,j !<[in,out] index we will swap if necessary
integer :: swp
if (arr(j) < arr(i)) then
  swp=i
  i=j
  j=swp
end if
end subroutine icomp_xchg

end program nautilus_major_reactions
