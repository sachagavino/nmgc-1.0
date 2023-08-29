program nautilus_select_outputs

use numerical_types
use iso_fortran_env
use utilities
use nautilus_main

implicit none

! Locals
character(len=80) :: filename_output
integer :: species, output, idx_1D,j ! index for loops
logical :: isDefined
logical :: wrong_species
integer :: use_ans
integer :: use_spe_nb
integer :: ans, i, resol
integer :: x_evol
character(len=11) :: user_species
character(len=11), dimension(:), allocatable :: user_name
integer, dimension(:), allocatable :: user_species_id
real(double_precision), dimension(:,:,:), allocatable :: user_ab

character(len=80) :: output_format !< format used to output data

real(double_precision), dimension(:,:,:), allocatable :: abundances_out !< abundances over time for each species. (nb_outputs, nb_species)

real(double_precision), dimension(:), allocatable :: time !< Simulation time [s]
real(double_precision), dimension(:,:), allocatable :: gas_temperature_out !< [K]
real(double_precision), dimension(:,:), allocatable :: dust_temperature_out !< [K]
real(double_precision), dimension(:,:), allocatable :: density !< [part/cm^3] 
real(double_precision), dimension(:,:), allocatable :: visual_extinction_out !< visual extinction [mag]
real(double_precision), dimension(:), allocatable :: x_rate !< X ionisation rate [s-1]

integer, dimension(:), allocatable :: compt

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

write(*,'(a,i0)') 'Spatial resolution: ', spatial_resolution
write(*,'(a,i0)') 'Number of time outputs: ', nb_outputs
write(*,'(a,i0)') 'Number of species: ', nb_species

! We allocate the output arrays
allocate(time(nb_outputs))
allocate(gas_temperature_out(spatial_resolution, nb_outputs))
allocate(dust_temperature_out(spatial_resolution, nb_outputs))
allocate(density(spatial_resolution, nb_outputs))
allocate(visual_extinction_out(spatial_resolution, nb_outputs))
allocate(x_rate(nb_outputs))

allocate(abundances_out(nb_outputs, nb_species, spatial_resolution))

allocate(compt(nb_species))
allocate(user_species_id(nb_species))

! The next write will be written in the same line
write(*,'(a)', advance='no') 'Reading unformatted outputs...'
! We read output files
do output=1,nb_outputs
  write(filename_output, '(a,i0.6,a)') 'abundances.',output,'.out'

  open(10, file=filename_output, status='old', form='unformatted')
  read(10) time(output)
  read(10) gas_temperature_out(1:spatial_resolution, output), dust_temperature_out(1:spatial_resolution, output), &
           density(1:spatial_resolution, output), &
           visual_extinction_out(1:spatial_resolution, output), x_rate(output)
  read(10) abundances_out(output,1:nb_species, 1:spatial_resolution)
  close(10)
enddo
! achar(13) is carriage return '\r'. Allow to go back to the beginning of the line
write(*,'(a,a)') achar(13), 'Reading unformatted outputs... Done'

!####################################################@@
! This part is to write one file per species, each line being one output time
!####################################################@@
write(*,*) " "
write(*,'(a)') 'Do you want to keep track of your commands?'
write(*,'(a)') 'No: 0   Yes: 1'
read(*,*) ans

if(ans .ne. 0) then
  open (10,file="selected_output.in",status='unknown')
  write(10,*) ans
endif

user_species_id = -1
if (spatial_resolution.eq.1) then
  use_ans = 1
  use_spe_nb = 0
  do while (use_ans==1)
     wrong_species = .true.
     use_spe_nb = use_spe_nb + 1
     do while (wrong_species)
        write(*,'(a)') 'Which species you want?'
        read(*,*) user_species
        if(ans .ne. 0) then
           write(10,*) user_species
        endif

        do i=1, nb_species
           if (user_species.eq.species_name(i)) then
               user_species_id(use_spe_nb) = i
           endif
        enddo

        if (user_species_id(use_spe_nb).eq.-1) then
            write(*,*) "'", trim(user_species), "' doesn't exist"
        else
            wrong_species = .false.
        endif
     enddo
     write(*,'(a)') "Continue?"
     write(*,'(a)') 'No: -1   Yes: 1'
     read(*,*) use_ans
     if(ans .ne. 0) then
       write(10,*) use_ans
     endif
  enddo
  allocate(user_ab(nb_outputs,use_spe_nb,spatial_resolution))
  allocate(user_name(use_spe_nb))
  do i=1,use_spe_nb
     user_ab(:,i,spatial_resolution) = abundances_out(:,user_species_id(i), spatial_resolution)
     user_name(i) = species_name(user_species_id(i))
  enddo

  open(20,file="selected_output.dat", status="unknown")
  write(output_format, *) '(a16,',use_spe_nb,'(a16))'
  write(20,output_format) "Time [yr]", user_name(:)
  write(20,*) ''
  write(output_format, *) '(es16.8,',use_spe_nb,'(es16.8))'
  do output=1,nb_outputs
      write(20,output_format) time(output)/year, user_ab(output, :, spatial_resolution)
  enddo
  close(20)
else
  write(*,'(a)') 'Evolution with:'
  write(*,'(a)') '1: Av   2: Time'
  read(*,*) x_evol
  if(ans .ne. 0) then
     write(10,*) x_evol
  endif
  use_ans = 1
  use_spe_nb = 0

  do while (use_ans==1)
     wrong_species = .true.
     use_spe_nb = use_spe_nb + 1
     do while (wrong_species)
        write(*,'(a)') 'Which species you want?'
        read(*,*) user_species
        if(ans .ne. 0) then
           write(10,*) user_species
        endif

        do i=1, nb_species
           if (user_species.eq.species_name(i)) then
               user_species_id(use_spe_nb) = i
           endif
        enddo

        if (user_species_id(use_spe_nb).eq.-1) then
            write(*,*) "'", trim(user_species), "' doesn't exist"
        else
            wrong_species = .false.
        endif
     enddo
     write(*,'(a)') "Continue?"
     write(*,'(a)') 'No: -1   Yes: 1'
     read(*,*) use_ans
     if(ans .ne. 0) then
        write(10,*) use_ans
     endif
  enddo

  allocate(user_ab(nb_outputs,use_spe_nb,spatial_resolution))
  allocate(user_name(use_spe_nb))

  do i=1,use_spe_nb
     user_ab(:,i,:) = abundances_out(:,user_species_id(i),:)
     user_name(i) = species_name(user_species_id(i))
  enddo

  open(20,file="selected_output.dat", status="unknown")
  if(x_evol==1) then

     write(output_format, *) '(2a16,',use_spe_nb,'(a16))'
     write(20,output_format) "Av [mag]","Time [yr]", user_name(:)
     write(20,*) ''
     write(output_format, *) '(2es16.8,',use_spe_nb,'(es16.8))'

     do output = 1, nb_outputs
        do resol=1,spatial_resolution
           write(20,output_format) visual_extinction_out(resol, output), time(output)/year, user_ab(output, :, resol)
        enddo
        write(20,*) ''
     enddo

     close(20)

  elseif(x_evol==2) then

     write(output_format, *) '(2a16,',use_spe_nb,'(a16))'
     write(20,output_format) "Time [yr]", "Av [mag]", user_name(:)
     write(20,*) ''
     write(output_format, *) '(2es16.8,',use_spe_nb,'(es16.8))'


     do resol=1,spatial_resolution
        do output = 1, nb_outputs
           write(20,output_format) time(output)/year,visual_extinction_out(resol, output), user_ab(output, :, resol)
        enddo
        write(20,*) ''
     enddo

     close(20)

   endif
endif

if(ans .ne. 0) then
   close(10)
endif


end program nautilus_select_outputs
