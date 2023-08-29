program nautilus_trace_major

use numerical_types
use iso_fortran_env
use utilities
use nautilus_main

implicit none

! Parameters
real(double_precision), parameter :: PERCENTAGE_THRESHOLD = 1.d0 !< Percentage below which the reaction will not be displayed

! Locals
character(len=80) :: filename_output
integer :: output, reaction, i, i_back ! index for loops
logical :: isDefined
real(double_precision) :: tmp !< temporary variable

! /!\ Variable names with _out are variable that already exist in the nautilus code, but here they are arrays, one value per output.

real(double_precision), dimension(:,:,:), allocatable :: abundances_out !< abundances_out over time for each species. (nb_outputs, nb_species)

real(double_precision), dimension(:), allocatable :: time
real(double_precision), dimension(:,:), allocatable :: gas_temperature_out
real(double_precision), dimension(:,:), allocatable :: dust_temperature_out
real(double_precision), dimension(:,:), allocatable :: density
real(double_precision), dimension(:,:), allocatable :: visual_extinction_out !< visual extinction
real(double_precision), dimension(:), allocatable :: zeta

! For rates
real(double_precision), allocatable, dimension(:,:,:) :: reaction_rates_out ! (spatial_reesolution,nb_outputs, nb_reactions)

! For storing info about production and destruction reaction of a given species through time
integer :: major_prod_size = 512 ! Default minimum size of some arrays
integer :: major_dest_size = 512 ! Default minimum size of some arrays
integer, allocatable, dimension(:) :: trace_prod_reaction !< The list of reaction ID's that are relevant to produce a given species
integer, allocatable, dimension(:) :: trace_dest_reaction !< The list of reaction ID's that are relevant to destruct a given species
real(double_precision), allocatable, dimension(:,:) :: trace_prod_percentage !< dim(nb_outputs, major_prod_size) the percentage of 
!! total production for the n-th species, the index being the same as in trace_prod_reaction
real(double_precision), allocatable, dimension(:,:) :: trace_dest_percentage !< dim(nb_outputs, major_prod_size) the percentage of 
!! total destruction for the n-th species, the index being the same as in trace_dest_reaction
integer :: nb_major_prod !< Actual number of relevant production reaction (less than 'major_prod_size')
integer :: nb_major_dest !< Actual number of relevant destruction reaction (less than 'major_dest_size')
integer :: store_id !< Corresponding index of the current reaction ID in the arrays

! Temp values to increase array size if needed
integer :: error !< In case of errors while allocating
integer, allocatable, dimension(:) :: reaction_temp
real(double_precision), allocatable, dimension(:,:) :: percentage_temp

! User asked values
logical :: wrong_species, wrong_1D !< Flags for while loops when asking the user something
character(len=11) :: user_species !< The species designed by the user
integer :: user_species_id !< corresponding id of the desired species of the user
integer :: output_id !< for loops
integer :: user_1D_id !< designed spatial id by the user

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

! Which species ?
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


! What spatial point ?

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


!######################################################
! End User asked section
!######################################################

allocate(trace_prod_reaction(major_prod_size))
allocate(trace_dest_reaction(major_dest_size))
trace_prod_reaction(major_prod_size) = 0
trace_dest_reaction(major_dest_size) = 0

allocate(trace_prod_percentage(nb_outputs, major_prod_size))
allocate(trace_dest_percentage(nb_outputs, major_dest_size))

! Initialize percentage to 0, because for reaction 
!! that do not exist from the beginning, we need to have a complete line of values. 
trace_prod_percentage(1:nb_outputs, 1:major_prod_size) = 0.d0
trace_dest_percentage(1:nb_outputs, 1:major_dest_size) = 0.d0

nb_major_prod = 0 ! To count the number of relevant reactions to make sure we do not cross the array size
nb_major_dest = 0 ! To count the number of relevant reactions to make sure we do not cross the array size

do output_id=1, nb_outputs

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
          tmp = tmp * abundances_out(output_id, REACTION_COMPOUNDS_ID(i, reaction), 1) * density(user_1D_id, output_ID)
        endif
      enddo
      destructions(reaction) = reaction_rates_out(user_1D_id,output_id, reaction) * tmp
    endif
    
    if (any(REACTION_COMPOUNDS_ID(MAX_REACTANTS+1:MAX_COMPOUNDS, reaction).eq.user_species_id)) then
      tmp = 1.d0
      do i=1, MAX_REACTANTS
        ! Skip if blanck species
        if (REACTION_COMPOUNDS_ID(i, reaction).ne.nb_species+1) then
          tmp = tmp * abundances_out(output_id, REACTION_COMPOUNDS_ID(i, reaction), 1) * density(user_1D_id, output_ID)
        endif
      enddo
      productions(reaction) = reaction_rates_out(user_1D_id,output_id, reaction) * tmp
    endif
  enddo
    
  call get_sorted_index(productions(1:nb_reactions), productions_id(1:nb_reactions))
  call get_sorted_index(destructions(1:nb_reactions), destructions_id(1:nb_reactions))

  productions_sum = sum(productions(1:nb_reactions))
  destructions_sum = sum(destructions(1:nb_reactions))



  !--------------- PRODUCTION (cm-3 s-1) -----------------    ------
  percentage = 100.d0
  i_back = nb_reactions ! We're counting backward
  
  ! For first element
  reaction = productions_id(i_back)
  percentage = productions(reaction) / productions_sum * 100.d0
  
  do while(percentage.gt.PERCENTAGE_THRESHOLD)
    
    
    if (nb_major_prod.eq.major_prod_size) then
      ! If the limit of the array is reach, we copy the values in a temporary array, allocate with a double size, et paste the 
      ! old values in the new bigger array
      allocate(reaction_temp(major_prod_size), stat=error)
      allocate(percentage_temp(nb_outputs,major_prod_size), stat=error)
      reaction_temp(1:major_prod_size) = trace_prod_reaction(1:major_prod_size)
      percentage_temp(1:nb_outputs, 1:major_prod_size) = trace_prod_percentage(1:nb_outputs, 1:major_prod_size)
      
      deallocate(trace_prod_reaction, stat=error)
      deallocate(trace_prod_percentage, stat=error)
      
      major_prod_size = major_prod_size * 2
      
      allocate(trace_prod_reaction(major_prod_size), stat=error)
      allocate(trace_prod_percentage(nb_outputs, major_prod_size), stat=error)
      ! Initialize percentage to 0, because for reaction 
      !! that do not exist from the beginning, we need to have a complete line of values. 
      trace_prod_reaction(1:major_prod_size) = 0
      trace_prod_percentage(1:nb_outputs,1:major_prod_size) = 0.d0
      trace_prod_reaction(1:major_prod_size/2) = reaction_temp(1:major_prod_size/2)
      trace_prod_percentage(1:nb_outputs,1:major_prod_size/2) = percentage_temp(1:nb_outputs,1:major_prod_size/2)
      
      deallocate(reaction_temp, stat=error)
      deallocate(percentage_temp, stat=error)
    end if
    
    store_id = nb_major_prod + 1 ! By default, we consider the current reaction do not exist in the list. 
    ! We test if it exist
    do i = 1, nb_major_prod
      if (trace_prod_reaction(i).eq.reaction) then
        store_id = i
        exit
      endif
    end do
    
    ! If we add a new species, we increment the counter
    if (store_id.eq.nb_major_prod + 1) then
      nb_major_prod = nb_major_prod + 1
    endif
    
    trace_prod_reaction(store_id) = reaction
    trace_prod_percentage(output_id, store_id) = percentage
    
    
    ! Decrement index
    i_back = i_back - 1
    
    ! Just in case, to avoid negative index for arrays
    if (i_back.lt.0) then
      percentage = 0.d0
    else
      ! Done at the end of the loop to avoid printing one extra reaction below the threshold
      reaction = productions_id(i_back)
      percentage = productions(reaction) / productions_sum * 100.d0
    endif
  enddo

  !--------------- DESTRUCTION (cm-3 s-1) ----------------    ------
  percentage = 100.d0
  i_back = nb_reactions ! We're counting backward
  
  ! For first element
  reaction = destructions_id(i_back)
  percentage = destructions(reaction) / destructions_sum * 100.d0
  
  do while(percentage.gt.PERCENTAGE_THRESHOLD)
    
    if (nb_major_dest.eq.major_dest_size) then
      ! If the limit of the array is reach, we copy the values in a temporary array, allocate with a double size, et paste the 
      ! old values in the new bigger array
      allocate(reaction_temp(major_dest_size), stat=error)
      allocate(percentage_temp(nb_outputs,major_dest_size), stat=error)
      reaction_temp(1:major_dest_size) = trace_dest_reaction(1:major_dest_size)
      percentage_temp(1:nb_outputs, 1:major_dest_size) = trace_dest_percentage(1:nb_outputs, 1:major_dest_size)
      
      deallocate(trace_dest_reaction, stat=error)
      deallocate(trace_dest_percentage, stat=error)
      
      major_dest_size = major_dest_size * 2
      
      allocate(trace_dest_reaction(major_dest_size), stat=error)
      allocate(trace_dest_percentage(nb_outputs, major_dest_size), stat=error)
      ! Initialize percentage to 0, because for reaction 
      !! that do not exist from the beginning, we need to have a complete line of values. 
      trace_dest_reaction(1:major_dest_size) = 0
      trace_dest_percentage(1:nb_outputs,1:major_dest_size) = 0.d0
      trace_dest_reaction(1:major_dest_size/2) = reaction_temp(1:major_dest_size/2)
      trace_dest_percentage(1:nb_outputs,1:major_dest_size/2) = percentage_temp(1:nb_outputs,1:major_dest_size/2)
      
      deallocate(reaction_temp, stat=error)
      deallocate(percentage_temp, stat=error)
    end if
    
    store_id = nb_major_dest + 1 ! By default, we consider the current reaction do not exist in the list. 
    ! We test if it exist
    do i = 1, nb_major_dest
      if (trace_dest_reaction(i).eq.reaction) then
        store_id = i
        exit
      endif
    end do
    
    ! If we add a new species, we increment the counter
    if (store_id.eq.nb_major_dest + 1) then
      nb_major_dest = nb_major_dest + 1
    endif
    
    trace_dest_reaction(store_id) = reaction
    trace_dest_percentage(output_id, store_id) = percentage

    
    ! Decrement index
    i_back = i_back - 1
    
    ! Just in case, to avoid negative index for arrays
    if (i_back.lt.0) then
      percentage = 0.d0
    else
      ! Done at the end of the loop to avoid printing one extra reaction below the threshold
      reaction = destructions_id(i_back)
      percentage = destructions(reaction) / destructions_sum * 100.d0
    endif
  enddo

enddo

! We store percentages for each reaction involved
write(filename_output, '(a,a,a)') 'trace_prod_',trim(user_species), '.percentage'
open(11, file=filename_output)
write(11,*) '! first line is the list of reaction ID involved.'
write(11,*) trace_prod_reaction(1:nb_major_prod)
write(11,*) '! time [year] ; one column per reaction, the percentage of production via this reaction.'
do output=1,nb_outputs
  write(11,*) time(output)/YEAR, trace_prod_percentage(output, 1:nb_major_prod)
enddo
close(11)

! We store percentages for each reaction involved
write(filename_output, '(a,a,a)') 'trace_dest_',trim(user_species), '.percentage'
open(11, file=filename_output)
write(11,*) '! first line is the list of reaction ID involved.'
write(11,*) trace_dest_reaction(1:nb_major_dest)
write(11,*) '! time [year] ; one column per reaction, the percentage of destruction via this reaction.'
do output=1,nb_outputs
  write(11,*) time(output)/YEAR, trace_dest_percentage(output, 1:nb_major_dest)
enddo
close(11)

write(filename_output, '(a,a,a)') 'trace_prod_',trim(user_species), '.reaction'
open(11, file=filename_output)
write(11,*) '! For each reaction involved, the reaction corresponding to a given ID.'
do i=1,nb_major_prod
  reaction = trace_prod_reaction(i)
  call display_reaction(REACTION_COMPOUNDS_NAMES(1:MAX_COMPOUNDS, reaction), reaction_line)
  write(11,'(i0,": ",a)') reaction, trim(reaction_line)
enddo
close(11)

write(filename_output, '(a,a,a)') 'trace_dest_',trim(user_species), '.reaction'
open(11, file=filename_output)
write(11,*) '! For each reaction involved, the reaction corresponding to a given ID.'
do i=1,nb_major_dest
  reaction = trace_dest_reaction(i)
  call display_reaction(REACTION_COMPOUNDS_NAMES(1:MAX_COMPOUNDS, reaction), reaction_line)
  write(11,'(i0,": ",a)') reaction, trim(reaction_line)
enddo
close(11)

write(*,*) 'ASCII Files have been generated. '
write(*,*) 'Please run nautilus-trace-species.py for the same species you asked here.'

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

end program nautilus_trace_major
