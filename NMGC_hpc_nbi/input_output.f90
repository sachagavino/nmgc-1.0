!******************************************************************************
! MODULE: input_output
!******************************************************************************
!
! DESCRIPTION: 
!> @brief Module that contains all the routines linked to reading input files
!! or writing output files. \n\n
!!
!! Input files tends to be named *.in
!! Output files tends to be named *.out
!! Temporary files that are overwritten at each timestep are named *.tmp
!
!******************************************************************************

module input_output

use iso_fortran_env


implicit none

contains


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Christophe Cossou
!
!> @date 2014
!
! DESCRIPTION: 
!> @brief Write information in an output file, notably the commit and branch of the compiled binary
!! Also show if the current version had uncommitted modification that can't
!! be traced.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine write_general_infos()
use global_variables
use git_infos

  implicit none
  
  character(len=80) :: filename = 'info.out'
  
  open(10, file=filename)
  write(10,'(a)') '!----------------------------'
  write(10,'(a)') '!     Nautilus Version      |'
  write(10,'(a)') '!----------------------------'
  write(10,'(a,a)') 'branch = ', trim(branch)
  write(10,'(a,a)') 'commit = ', trim(commit)
  write(10,'(a,a)') '!', trim(modifs)
  write(10,'(a)') ""
  write(10,'(a)') '!----------------------------'
  write(10,'(a)') '!      General infos        |'
  write(10,'(a)') '!----------------------------'
  write(10,'(a,i0)') 'Maximum number of non-zeros values in jacobian = ', nb_nonzeros_values
  close(10)
  
end subroutine write_general_infos

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Franck Hersant
!
!> @date 2000
!
! DESCRIPTION: 
!> @brief Write the list of species and the corresponding index in an
!! output file 'species.out'.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine write_species()
use global_variables

implicit none

! Locals
integer :: i

open(10, file='species.out')
! Write 'ggo_spec.d': 5 columns of numbered species=====================
write(10,'(5(I4,")",1X,A11,1X))') (I,species_name(I),I=1,nb_species)
close(10)

return
end subroutine write_species

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Christophe Cossou
!
!> @date 2014
!
! DESCRIPTION: 
!> @brief Write the list of elemental species with their abundances and mass.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine write_elemental_abundances(filename, el_abundances)

use global_variables

implicit none

! Input
character(len=*), intent(in) :: filename !< [in] the name of the output file
real(double_precision), intent(in), dimension(NB_PRIME_ELEMENTS) :: el_abundances !< [in] Elemental abundances, either initial or current ones

! Locals
integer :: i

open(10, file=filename)
write(10, '(a)') '! Species name ; abundance (relative to H) ; mass (AMU)'
do i=1, NB_PRIME_ELEMENTS
  
  write(10, '(a,es10.4e2, f8.3)') species_name(PRIME_ELEMENT_IDX(i)), el_abundances(i), elemental_mass(i)
  
enddo
close(10)

return
end subroutine write_elemental_abundances




!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Franck Hersant & Christophe Cossou
!
!> @date 2000
!
! DESCRIPTION: 
!> @brief Write all abundances for all species in an output file at each
!! output time. The total number of output files related to abundances 
!! will be equal to the number of timestep, not equally spaced in time.\n\n
!! Output filename is of the form : abundances.000001.out
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine write_current_output(index)
! Writes 1D outputs
use global_variables

implicit none

! Input
integer, intent(in) :: index !<[in] The reference index of the current output

! Locals
character(len=80) :: filename_output

write(filename_output, '(a,i0.6,a)') 'abundances.',index,'.out'


open(UNIT=35, file=filename_output, form='unformatted')

write(35) current_time
write(35) gas_temperature(1:spatial_resolution), dust_temperature(1,1:spatial_resolution), &
          H_number_density(1:spatial_resolution), visual_extinction(1:spatial_resolution), X_IONISATION_RATE
write(35) abundances(1:nb_species, 1:spatial_resolution)
close(35)

return
end subroutine write_current_output

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Franck Hersant
!
!> @date 2000
!
! DESCRIPTION: 
!> @brief Write rates of all chemical reactions for the current timestep.
!! The total number of files will be equal to the total number of timesteps, the
!! routine being called at the end of each timestep.\n\n
!! Output filename is of the form : rates.000001.out
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine write_current_rates(index)

use global_variables

implicit none

! Input
integer, intent(in) :: index !<[in] The reference index of the current output

! Locals
character(len=80) :: filename_output

write(filename_output, '(a,i0.6,a)') 'rates.',index,'.out'

open(45, file=filename_output, form='unformatted')

write(45) reaction_rates_1D 

close(45)

return 
end subroutine write_current_rates

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author
!> Valentine Wakelam
!
!> @date 2014
!
! DESCRIPTION:
!> @brief Write the H2 and CO column density computed by the model - used 
!   for the self-shielding of H2 and CO from the UV photons
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine write_H2_CO_col_dens(index)

use global_variables

implicit none

! Input
integer, intent(in) :: index !<[in] The reference index of the current output

! Locals
character(len=80) :: filename_output
integer :: i

write(filename_output, '(a,i0.6,a)') 'col_dens.',index,'.out'

open(55, file=filename_output)

! Header
write(55,'(50a)') 'H column density (cm-2) H2 column density (cm-2)  CO column density (cm-2)  N2 column density (cm-2)'

do i=1,spatial_resolution
    write(55,*) NH_z(i),NH2_z(i),NCO_z(i), NN2_z(i)
enddo

close(55)

return
end subroutine write_H2_CO_col_dens


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Franck Hersant
!
!> @date 2000
!
! DESCRIPTION: 
!> @brief Write the total chemical composition of the actual timestep in 
!! a file whose name is given as an input parameter. This allow to use the
!! same routine to write input, temporary and output files
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine write_abundances(filename)

use global_variables

implicit none

! Input
character(len=*), intent(in) :: filename !<[in] the name of the output file

! Locals
integer :: i
character(len=80) :: line_format

open(13, file=filename)

! Header
write(13,'(5(a,es10.3e2),a)') '!time =', current_time, ' s ; density = ', H_number_density, &
' part/cm^3 ; temperature=', gas_temperature,' K ; visual extinction = ', visual_extinction, &
' [mag] ; CR ionisation rate = ',CR_IONISATION_RATE,' s-1'

! To adapt the format in function of the 1D number of points
write(line_format,'(a,i0,a)') '(a," = ",', spatial_resolution, '(es12.5e2))'

do i=1,nb_species
  write(13,line_format) trim(species_name(i)), abundances(i,1:spatial_resolution)
enddo

close(13)

return
end subroutine write_abundances



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Christophe Cossou
!
!> @date 2014
!
! DESCRIPTION: 
!> @brief subroutine that write the simulation parameters into the file 'parameters.out'
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine write_parameters()
use global_variables

  implicit none
  
  character(len=80) :: filename = 'parameters.in'
  character(len=80) :: cp_command
  
  ! Copy the old file into a *.bak version just in case
  write(cp_command, '(5a)') 'cp ', trim(filename), ' ', trim(filename), '.bak'
  call system(cp_command)
  
  open(10, file=filename)
  write(10,'(a)') "!# ------------------------------------------------"
  write(10,'(a)') "!# Parameter file for various properties of the disk."
  write(10,'(a)') "!# ------------------------------------------------"
  write(10,'(a)') "!# blanck line or with spaces will be skipped."
  write(10,'(a)') "!# In fact, the only lines that matter are non commented lines with a"
  write(10,'(a)') "!# '=' character to distinguish the identificator and the value(s)"
  write(10,'(a)') "!# (each value must be separated with at least one space."
  write(10,'(a)') "!# Line must not be longer than 80 character, but comments can be far"
  write(10,'(a)') "!# bigger than that, even on line with a parameter to read."
  write(10,'(a)') ""
  write(10,'(a)') "!*****************************"
  write(10,'(a)') "!*   Switch 2/3 phase model  *"
  write(10,'(a)') "!*****************************"
  write(10,'(a)') ""
  write(10,'(a,i0,a)') 'is_3_phase = ', is_3_phase, ' ! 0: 2 phase, 1: 3 phase'
  write(10,'(a)') ""
  write(10,'(a)') "!*****************************"
  write(10,'(a)') "!*          Switches         *"
  write(10,'(a)') "!*****************************"
  write(10,'(a)') ""
  write(10,'(a,i0,a)') 'preliminary_test = ', IS_TEST, ' ! Will or will not do comprehensive tests &
  &before the simulation. Switch it off when lauching thousands or simulations'
  write(10,'(a,i0,a)') 'is_structure_evolution = ', IS_STRUCTURE_EVOLUTION, ' ! If 1, physical structure properties evolve with &
                        &time, values come from structure_evolution.dat file that must exists'
  write(10,'(a,a,a)') 'grain_temperature_type = ', trim(GRAIN_TEMPERATURE_TYPE), ' ! fixed, fixed_to_dust_size, gas,&
                      & table_evolv, table_1D or computed'
  write(10,'(a)') '! fixed: Tgrain = initial_dust_temperature. All dust grains have same temperature;'
  write(10,'(a)') '! fixed_to_dust_size = each grain have a fixed temperature defined in1D_grain_sizes.in or 0D_grain_sizes.in;' 
  write(10,'(a)') '! gas: Tgrain = Tgas ; '
  write(10,'(a)') '! table_evolv: Tgrain is interpolated from structure_evolution.dat data file (5th optional column) ; '
  write(10,'(a)') '! table_1D: Tgrain is read in the 1D_static. dat file (5th column) ; '
  write(10,'(a)') '! computed: calculated from uv_flux and visual extinction by radiative equilibrium'
  ! write(10,'(a,i0,a)') 'is_dust_1D = ', is_dust_1D, ' ! Reading the grain abundance and the NH/AV factor in the 1D_static.dat file & 
                     ! &(mostly for disks)'
  write(10,'(a,i0,a)') 'photo_disk = ', photo_disk, ' ! Computation of photodissociation rates for protoplanetary disks.'
  write(10,'(a,i0,a)') 'is_grain_reactions = ', IS_GRAIN_REACTIONS, ' ! Accretion, grain surface reactions'
  write(10,'(a,i0,a)') 'is_h2_adhoc_form = ', IS_H2_ADHOC_FORM, ' ! Ad hoc formation of H2 on grain surfaces (1=activated)'
  write(10,'(a,i0,a)') 'is_h2_formation_rate = ', is_h2_formation_rate, ' ! h2 formation rates on surfaces from Bron et al: (2014)'
  write(10,'(a,i0,a)') 'height_h2formation = ', height_h2formation, ' ! Spatial point above which B14s method is used. If 0 then &
                      &B14 is not used at all.'
  write(10,'(a,i0,a)') 'is_absorption_h2 = ', is_absorption_h2, ' ! H2 self-shielding from Lee & Herbst (1996) (1=activated)'
  write(10,'(a,i0,a)') 'is_absorption_co = ', is_absorption_co, ' ! CO self-shielding. (1: Lee & Herbst (1996), &
                                                                 & 2: Visser et al. (2009)'
  write(10,'(a,i0,a)') 'is_absorption_n2 = ', is_absorption_n2, ' ! N2 self-shielding from Li et al. (2013) (1=activated)'
  write(10,'(a,i0,a)') 'is_photodesorb = ', is_photodesorb, &
' ! Switch to turn on the photodesorption of ices (default yield is 1e-4)'
  write(10,'(a,i0,a)') 'is_crid = ', is_crid, &
' ! Switch to turn on the CRID (cosmic rays induced diffusion) mechanism'
  write(10,'(a,i0,a)') 'is_er_cir = ', is_er_cir, &
' ! Switch to turn on Eley-Rideal and Complex Induced Reaction mechanisms (default=0: desactivated)'
  write(10,'(a,i0,a)') 'grain_tunneling_diffusion = ', GRAIN_TUNNELING_DIFFUSION, &
  ' ! 0=thermal; For H,H2: 1=QM1; 2=QM2; 3=choose fastest'
  write(10,'(a,i0,a)') 'modify_rate_flag = ', MODIFY_RATE_FLAG, ' ! 1=modify H; 2=modify H,H2, 3=modify all, -1=H+H only'
  write(10,'(a,i0,a)') 'conservation_type = ', CONSERVATION_TYPE, ' ! 0=only e- conserved; 1=elem #1 conserved, 2=elem #1 & #2, etc'
  !write(10,'(a,i0,a)') 'is_dust_MRN = ', is_dust_MRN, '!  0 = custom; 1 = MRN distribution; 2 = WD distribution. 0D mode only.'
  write(10,'(a)') ""
  write(10,'(a)') "!*****************************"
  write(10,'(a)') "!* Number of active layers   *"
  write(10,'(a)') "!*****************************"
  write(10,'(a)') ""
  write(10,'(a,es10.3e2,a)') 'nb_active_lay = ',nb_active_lay, ' ! Number of active layers'
  write(10,'(a)') ""
  write(10,'(a)') "!******************************"
  write(10,'(a)') "!*          0D or 1D          *"
  write(10,'(a)') "!******************************"
  write(10,'(a)') "!(diffusion is for species, not the structure)"
  write(10,'(a)') ""
  write(10,'(a,a,a)') 'structure_type = ', trim(STRUCTURE_TYPE), ' ! 0D, 1D_diff, 1D_no_diff'
  write(10,'(a,i0,a)') 'spatial_resolution = ', spatial_resolution, &
    ' ! Number of lines in 1D. If 1, we are in 0D, else, we are in 1D, with diffusion between gas boxes.'
  write(10,'(a)') ""
  write(10,'(a)') '!******************************************************'
  write(10,'(a)') '!*        Use single-grain or multi-grain mode        *'
  write(10,'(a)') '!******************************************************'
  write(10,'(a)') ""
  write(10,'(a,i0,a)') 'multi_grain = ', multi_grain, ' ! 1 = multi-grain; 0 = single-grain. & 
                                               & If 1, then the grain parameters are read in 0D/1D_grain_sizes.in files.'
  write(10,'(a)') ""
  write(10,'(a)') "!*****************************"
  write(10,'(a)') "!*    Gas phase parameters   *"
  write(10,'(a)') "!*****************************"
  write(10,'(a)') ""
  write(10,'(a,es10.3e2,a)') 'initial_gas_density = ', initial_gas_density, ' ! initial gas density [part/cm-3]'
  write(10,'(a,es10.3e2,a)') 'initial_gas_temperature = ', initial_gas_temperature, ' ! initial gas temperature [K]'
  write(10,'(a,es10.3e2,a)') 'initial_visual_extinction = ', INITIAL_VISUAL_EXTINCTION, ' ! initial visual extinction'
  write(10,'(a,es10.3e2,a)') 'cr_ionisation_rate = ', CR_IONISATION_RATE, ' ! cosmic ray ionisation rate [s-1] (standard=1.3e-17)'
  write(10,'(a,es10.3e2,a)') 'x_ionisation_rate = ', X_IONISATION_RATE, ' ! Ionisation rate due to X-rays [s-1]'
  write(10,'(a,es10.3e2,a)') 'uv_flux = ', UV_FLUX, ' ! Scale factor for the UV flux, in unit of the reference flux (1.=nominal)'
  write(10,'(a)') ""
  write(10,'(a)') "!*****************************"
  write(10,'(a)') "!*      Grain parameters     *"
  write(10,'(a)') "!*****************************"
  write(10,'(a)') ""
  write(10,'(a,es10.3e2,a)') 'initial_dust_temperature = ', initial_dust_temperature, ' ! initial dust temperature [K]&
                              & when grain_temperature_type=fixed'
  write(10,'(a,es10.3e2,a)') 'initial_dtg_mass_ratio = ', initial_dtg_mass_ratio, ' ! dust-to-gas ratio by mass'
  write(10,'(a,es10.3e2,a)') 'sticking_coeff_neutral = ', sticking_coeff_neutral, ' ! sticking coeff for neutral species'
  write(10,'(a,es10.3e2,a)') 'sticking_coeff_positive = ', sticking_coeff_positive, ' ! sticking coeff for positive species'
  write(10,'(a,es10.3e2,a)') 'sticking_coeff_negative = ', sticking_coeff_negative, ' ! sticking coeff for negative species'
  write(10,'(a,es10.3e2,a)') 'grain_density = ', GRAIN_DENSITY, ' ! mass density of grain material'
  write(10,'(a,es10.3e2,a)') 'grain_radius = ', grain_radius, ' ! grain radius [cm]'
  write(10,'(a,es10.3e2,a)') 'diffusion_barrier_thickness = ', DIFFUSION_BARRIER_THICKNESS, ' ! Barrier thickness [cm]'
  write(10,'(a,es10.3e2,a)') 'surface_site_density = ', SURFACE_SITE_DENSITY, ' ! site density on one grain [cm-2]'
  write(10,'(a,es10.3e2,a)') 'diff_binding_ratio_surf = ', DIFF_BINDING_RATIO_SURF, &
                             ' ! Ratio used to compute the DIFFUSION_BARRIER from the BINDING_ENERGY if not known (surface species)'
  write(10,'(a,es10.3e2,a)') 'diff_binding_ratio_mant = ', DIFF_BINDING_RATIO_MANT, &
                             ' ! Ratio used to compute the DIFFUSION_BARRIER from the BINDING_ENERGY if not known  (mantle species)'
  write(10,'(a,es10.3e2,a)') 'chemical_barrier_thickness = ', CHEMICAL_BARRIER_THICKNESS, &
                             ' ! grain reaction activation energy barrier width. [cm]'
  write(10,'(a,es10.3e2,a)') 'cr_peak_grain_temp = ', CR_PEAK_GRAIN_TEMP, ' ! peak grain temperature [K] (CR heating)'
  write(10,'(a,es10.3e2,a)') 'cr_peak_duration = ', CR_PEAK_DURATION, ' ! duration [s] of peak grain temperature'
  write(10,'(a,es10.3e2,a)') 'Fe_ionisation_rate = ', FE_IONISATION_RATE, ' ! (cosmic) Fe-ion--grain encounter [s-1 grain-1] '
  write(10,'(a,es10.3e2,a)') '!! (for 0.1 micron grain) For cosmic photo desorptions, only Fe-ions are efficient to heat grains. '
  write(10,'(a,es10.3e2,a)') 'vib_to_dissip_freq_ratio = ', VIB_TO_DISSIP_FREQ_RATIO, &
                             ' ! [no unit] The ratio of the surface-molecule bond frequency to the frequency at'
  write(10,'(a)') '!! which energy is lost to the grain surface. Used for the RRK (Rice Ramsperger-Kessel) desorption mechanism'
  write(10,'(a)') '!! (see Garrod el al. 2007 for more). Assumed to be 1% by default.'
  write(10,'(a,es10.3e2,a)') 'ED_H2 = ', ED_H2, &
                            ' ! H2 binding energy over itself. Used for the desorption encounter mechanism. in K. '
  write(10,'(a)') ""
  write(10,'(a)') "!*****************************"
  write(10,'(a)') "!*  Integration and Outputs  *"
  write(10,'(a)') "!*****************************"
  write(10,'(a)') ""
  write(10,'(a,es10.3e2,a)') 'start_time = ', START_TIME/YEAR, ' ! [yrs] first output time'
  write(10,'(a,es10.3e2,a)') 'stop_time = ', STOP_TIME/YEAR, ' ! [yrs] last output time'
  write(10,'(a,i0,a)') 'nb_outputs = ', NB_OUTPUTS, ' ! Total number of outputs (used for linear or log spaced outputs)'
  write(10,'(a,a,a)') 'output_type = ', trim(OUTPUT_TYPE), ' ! linear, log'
  write(10, '(a)') '! linear: Output times are linearly spaced'
  write(10, '(a)') '! log   : Outputs times are log-spaced'
  write(10,'(a,es10.3e2, a)') 'relative_tolerance = ',RELATIVE_TOLERANCE, ' ! Relative tolerance of the solver'
  write(10,'(a,es10.3e2,a)') 'minimum_initial_abundance = ', MINIMUM_INITIAL_ABUNDANCE, ' ! default minimum initial &
                             &fraction abundance'
  close(10)
  
end subroutine write_parameters

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Wasim Iqbal
!
!> @date 2017
!
! DESCRIPTION: 
!> @brief reads grain_radii
!! from datafile 0D_grain_sizes.in
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine get_grain_radii() 
use global_variables
implicit none
! Locals
character(len=80) :: filename
character(len=80) :: filename_static
integer           :: i,j,nb_lines, nb_lines_static
character(len=1000) :: line
character(len=1), parameter :: comment_character = '!' ! character that will indicate that the rest of the line is a comment
integer :: comment_position ! the index of the comment character on the line. if zero, there is none on the current string
integer :: error ! to store the state of a read instruction
logical :: isDefined

!if (is_dust_1D.eq.0 ) then
if (STRUCTURE_TYPE.eq.'0D') then
  if (multi_grain.eq.0) then

    nb_grains = 1

    allocate(grain_radii(nb_grains))
    grain_radii(1:nb_grains) = 0.d0
    
    allocate(grain_temp(nb_grains))
    grain_temp(1:nb_grains) = 0.d0
    
    allocate(CR_PEAK_GRAIN_TEMP_all(nb_grains))
    CR_PEAK_GRAIN_TEMP_all(1:nb_grains) = 0.d0
    
    allocate(actual_dust_temp(nb_grains))
    actual_dust_temp(1:nb_grains) = 0.d0
    
    allocate(INDGRAIN(nb_grains))
    INDGRAIN(1:nb_grains) = 0
    
    allocate(INDGRAIN_MINUS(nb_grains))
    INDGRAIN_MINUS(1:nb_grains) = 0
    
    allocate(sumlaysurfsave(nb_grains))
    sumlaysurfsave(1:nb_grains) = 0.d0
    
    allocate(sumlaymantsave(nb_grains))
    sumlaymantsave(1:nb_grains) = 0.d0
    
    allocate(EVAPORATION_RATES_H2(nb_grains))
    EVAPORATION_RATES_H2(1:nb_grains) = 0.d0
    
    allocate(EVAPORATION_RATES_TEMPO_H2(nb_grains))
    EVAPORATION_RATES_TEMPO_H2(1:nb_grains) = 0.d0
    
    allocate(YGRAIN(nb_grains))
    YGRAIN(1:nb_grains) = ""
    
    allocate(YGRAIN_MINUS(nb_grains))
    YGRAIN_MINUS(1:nb_grains) = ""
    
    allocate(GTODN(nb_grains))
    GTODN(1:nb_grains) = 0.d0
    
    allocate(nb_sites_per_grain(nb_grains))
    nb_sites_per_grain(1:nb_grains) = 0.d0

    grain_radii(:)=grain_radius                  
    CR_PEAK_GRAIN_TEMP_all(:)=cr_peak_grain_temp  
    grain_temp(:)= initial_dust_temperature       

  elseif (multi_grain.eq.1) then
    filename='0D_grain_sizes.in'
  
    inquire(file=filename, exist=isDefined)
    
    if (isDefined) then
    
    call get_linenumber(filename, nb_lines)
    
    ! if (nb_grains.ne.nb_lines) then
    !   write(Error_unit,'(a)') 'Please check: nb_grains in parameters.in is different than the number of grains in 0D_grains_sizes.in'
    !   call exit(22)
    ! endif

    nb_grains = nb_lines

    allocate(grain_radii(nb_grains))
    grain_radii(1:nb_grains) = 0.d0
    
    allocate(grain_temp(nb_grains))
    grain_temp(1:nb_grains) = 0.d0

    allocate(GTODN_0D_temp(nb_grains))
    GTODN_0D_temp(1:nb_grains) = 0.d0
    
    allocate(CR_PEAK_GRAIN_TEMP_all(nb_grains))
    CR_PEAK_GRAIN_TEMP_all(1:nb_grains) = 0.d0
    
    allocate(actual_dust_temp(nb_grains))
    actual_dust_temp(1:nb_grains) = 0.d0
    
    allocate(INDGRAIN(nb_grains))
    INDGRAIN(1:nb_grains) = 0
    
    allocate(INDGRAIN_MINUS(nb_grains))
    INDGRAIN_MINUS(1:nb_grains) = 0
    
    allocate(sumlaysurfsave(nb_grains))
    sumlaysurfsave(1:nb_grains) = 0.d0
    
    allocate(sumlaymantsave(nb_grains))
    sumlaymantsave(1:nb_grains) = 0.d0
    
    allocate(EVAPORATION_RATES_H2(nb_grains))
    EVAPORATION_RATES_H2(1:nb_grains) = 0.d0
    
    allocate(EVAPORATION_RATES_TEMPO_H2(nb_grains))
    EVAPORATION_RATES_TEMPO_H2(1:nb_grains) = 0.d0
    
    allocate(YGRAIN(nb_grains))
    YGRAIN(1:nb_grains) = ""
    
    allocate(YGRAIN_MINUS(nb_grains))
    YGRAIN_MINUS(1:nb_grains) = ""
    
    allocate(GTODN(nb_grains))
    GTODN(1:nb_grains) = 0.d0
    
    allocate(nb_sites_per_grain(nb_grains))
    nb_sites_per_grain(1:nb_grains) = 0.d0
    
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
        read(line,*) grain_radii(i), GTODN_0D_temp(i), grain_temp(i), CR_PEAK_GRAIN_TEMP_all(i)
        i = i + 1
      endif
    enddo
    else   
      write(Error_unit,*) 'Error: The file ', trim(filename),' does not exist.'  
      call exit
    endif
    close(10)
  endif
endif
! do i=1,nb_grains
! write(*, *) i, grain_radii(i),grain_temp(i)
! enddo
! pause
!if (is_dust_1D.eq.1 ) then
if ((STRUCTURE_TYPE.eq.'1D_no_diff').or.(STRUCTURE_TYPE.eq.'1D_diff')) then  
  if (multi_grain == 0) then
    nb_grains = 1
    allocate(grain_radii(nb_grains))
    grain_radii(1:nb_grains) = 0.d0
    
    allocate(grain_temp(nb_grains))
    grain_temp(1:nb_grains) = 0.d0
    
    allocate(CR_PEAK_GRAIN_TEMP_all(nb_grains))
     CR_PEAK_GRAIN_TEMP_all(1:nb_grains) = 0.d0
    
    allocate(actual_dust_temp(nb_grains))
    actual_dust_temp(1:nb_grains) = 0.d0
    
    allocate(INDGRAIN(nb_grains))
    INDGRAIN(1:nb_grains) = 0
    
    allocate(INDGRAIN_MINUS(nb_grains))
    INDGRAIN_MINUS(1:nb_grains) = 0
    
    allocate(sumlaysurfsave(nb_grains))
    sumlaysurfsave(1:nb_grains) = 0.d0
    
    allocate(sumlaymantsave(nb_grains))
    sumlaymantsave(1:nb_grains) = 0.d0
    
    allocate(EVAPORATION_RATES_H2(nb_grains))
    EVAPORATION_RATES_H2(1:nb_grains) = 0.d0
    
    allocate(EVAPORATION_RATES_TEMPO_H2(nb_grains))
    EVAPORATION_RATES_TEMPO_H2(1:nb_grains) = 0.d0
    
    allocate(YGRAIN(nb_grains))
    YGRAIN(1:nb_grains) = ""
    
    allocate(YGRAIN_MINUS(nb_grains))
    YGRAIN_MINUS(1:nb_grains) = ""
    
    allocate(GTODN(nb_grains))
    GTODN(1:nb_grains) = 0.d0
    
    allocate(nb_sites_per_grain(nb_grains))
    nb_sites_per_grain(1:nb_grains) = 0.d0
    
    ! we do not read 1D_grain_sizes.in or 0D_grain_sizes.in, values are read from 1D_static.dat file
    grain_radii(:)=grain_radius                   !assigning value from parameter.in file, just to be sure that variable has some value 
    CR_PEAK_GRAIN_TEMP_all(:)=cr_peak_grain_temp  !assigning value from parameter.in file
    grain_temp(:)= initial_dust_temperature       !assigning value from parameter.in file
  elseif (multi_grain == 1) then
    filename='1D_grain_sizes.in'
    filename_static='1D_static.dat'
    inquire(file=filename, exist=isDefined)
    
    if (isDefined) then
  
      nb_grains = get_nb_columns(filename)/4.d0 !4 is the number of parameters in the file.

      call get_linenumber(filename, nb_lines)
      call get_linenumber(filename_static, nb_lines_static)

      if (nb_lines.ne.nb_lines_static) then
        write(Error_unit,'(a)') 'Please check: number of spatial points in 1D_static.dat different than that in 1D_grain_sizes.in'
        call exit(22)
      endif

      allocate(grain_radii(nb_grains))
      grain_radii(1:nb_grains) = 0.d0
      
      allocate(grain_radii_1D(nb_grains,spatial_resolution))
      grain_radii_1D(1:nb_grains,1:spatial_resolution) = 0.d0
      
      allocate(CR_PEAK_GRAIN_TEMP_all_1D(nb_grains,spatial_resolution))
      CR_PEAK_GRAIN_TEMP_all_1D(1:nb_grains,1:spatial_resolution) = 0.d0
      
      allocate(grain_temp_1D(nb_grains,spatial_resolution))
      grain_temp_1D(1:nb_grains,1:spatial_resolution) = 0.d0
      
      allocate(GTODN_1D_temp(nb_grains,spatial_resolution))
      GTODN_1D_temp(1:nb_grains,1:spatial_resolution) = 0.d0
      
      
      allocate(grain_temp(nb_grains))
      grain_temp(1:nb_grains) = 0.d0
      
      allocate(CR_PEAK_GRAIN_TEMP_all(nb_grains))
       CR_PEAK_GRAIN_TEMP_all(1:nb_grains) = 0.d0
      
      allocate(actual_dust_temp(nb_grains))
      actual_dust_temp(1:nb_grains) = 0.d0
      
      allocate(INDGRAIN(nb_grains))
      INDGRAIN(1:nb_grains) = 0
      
      allocate(INDGRAIN_MINUS(nb_grains))
      INDGRAIN_MINUS(1:nb_grains) = 0
      
      allocate(sumlaysurfsave(nb_grains))
      sumlaysurfsave(1:nb_grains) = 0.d0
      
      allocate(sumlaymantsave(nb_grains))
      sumlaymantsave(1:nb_grains) = 0.d0
      
      allocate(EVAPORATION_RATES_H2(nb_grains))
      EVAPORATION_RATES_H2(1:nb_grains) = 0.d0
      
      allocate(EVAPORATION_RATES_TEMPO_H2(nb_grains))
      EVAPORATION_RATES_TEMPO_H2(1:nb_grains) = 0.d0
      
      allocate(YGRAIN(nb_grains))
      YGRAIN(1:nb_grains) = ""
      
      allocate(YGRAIN_MINUS(nb_grains))
      YGRAIN_MINUS(1:nb_grains) = ""
      
      allocate(GTODN(nb_grains))
      GTODN(1:nb_grains) = 0.d0
      
      allocate(nb_sites_per_grain(nb_grains))
      nb_sites_per_grain(1:nb_grains) = 0.d0
      
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
          read(line,*) (grain_radii_1D(j,i),j=1,nb_grains),(GTODN_1D_temp(j,i),j=1,nb_grains),&
                       (grain_temp_1D(j,i),j=1,nb_grains),(CR_PEAK_GRAIN_TEMP_all_1D(j,i),j=1,nb_grains)
          !  write(*,'(10ES14.5)') (grain_radii_1D(j,i),j=1,nb_grains)
          !  write(*,'(40ES14.5)') (grain_radii_1D(j,i),GTODN_1D_temp(j,i),grain_temp_1D(j,i),CR_PEAK_GRAIN_TEMP_all_1D(j,i),j=1,nb_grains)
          !           pause
          i = i + 1
        endif
      enddo
      grain_radii(1:nb_grains)=grain_radii_1D(1:nb_grains,1)
      grain_temp(:)=grain_temp_1D(:,1)
      CR_PEAK_GRAIN_TEMP_all(:)=CR_PEAK_GRAIN_TEMP_all_1D(:,1)
      
!        write(*,*)'----',i,(grain_radii_1D(j,i),j=1,nb_grains)
!        write(*,*)grain_radii
!        write(*,*)grain_radii_1D(:,1)
!        pause
      else   
        write(Error_unit,*) 'Error: The file ', trim(filename),' does not exist.'  
        call exit
      endif
      close(10)        
      
     
  endif
endif
return
end subroutine get_grain_radii

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Wasim Iqbal
!
!> @date 2017
!
! DESCRIPTION: 
!> @brief assigning values to character variable YGRAIN and YGRAIN_MINUS
!! 
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

subroutine get_YGRAIN()
use global_variables
implicit none
integer :: i
character(2) :: i_c !character variable to convert integer value of j to character j

YGRAIN = ""
YGRAIN_MINUS = ""

do i=1,nb_grains
  write(i_c,'(I2.2)')i
  YGRAIN(i)      = "GRAIN"//trim(i_c)
  YGRAIN_MINUS(i)= "GRAIN"//trim(i_c)//"-"
enddo
!  write(*,*)(YGRAIN(I),I=1,nb_grains)
!   write(*,*)(YGRAIN_MINUS(I),I=1,nb_grains)
end subroutine get_YGRAIN


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Christophe Cossou
!
!> @date 2014
!
! DESCRIPTION: 
!> @brief Read name and mass of all 'basic' elements (basic brick for molecules
!! such as H, He and so on).
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine read_element_in()
! Writes 1D outputs
use global_variables

implicit none

! Locals
character(len=80), parameter :: filename='element.in'
integer :: i

character(len=80) :: line
character(len=1), parameter :: comment_character = '!' ! character that will indicate that the rest of the line is a comment
integer :: comment_position ! the index of the comment character on the line. if zero, there is none on the current string
integer :: error ! to store the state of a read instruction

logical :: isDefined

character(len=80) :: filename_gas = 'gas_species.in'
character(len=80) :: filename_grain = 'grain_species.in'

integer :: nb_columns_gas
integer :: nb_columns_grain

inquire(file=filename, exist=isDefined)

call get_linenumber(filename, NB_PRIME_ELEMENTS)

! We get the number of prime elements
nb_columns_gas = get_nb_columns(filename_gas)

if ((nb_columns_gas - 2).ne.NB_PRIME_ELEMENTS) then
  write (Error_unit,'(a,i0,a,a,a,i0,a)') 'The number of prime elements is different in "element.in" (', NB_PRIME_ELEMENTS, &
  ') and "', trim(filename_gas), '" (', nb_columns_gas-2, ') .'
  call exit(6)
endif

nb_columns_grain = get_nb_columns(filename_grain)

if ((nb_columns_grain - 2).ne.NB_PRIME_ELEMENTS) then
  write (Error_unit,'(a,i0,a,a,a,i0,a)') 'The number of prime elements is different in "element.in" (', NB_PRIME_ELEMENTS, &
  ') and "', trim(filename_grain), '" (', nb_columns_grain-2, ') .'
  call exit(6)
endif

! We allocate global variables
allocate(element_name(NB_PRIME_ELEMENTS))
allocate(elemental_mass(NB_PRIME_ELEMENTS))

element_name(1:NB_PRIME_ELEMENTS) = ''
elemental_mass(1:NB_PRIME_ELEMENTS) = 0.d0

if (isDefined) then

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
      read(line, '(a, f8.3)') element_name(i), elemental_mass(i)
      i = i + 1
    endif
  enddo
  close(10)
endif

! Other operations are done once species from gas and grain are read, because we need the indexes of all prime element in those arrays

return
end subroutine read_element_in

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Christophe Cossou
!
!> @date 2014
!
! DESCRIPTION: 
!> @brief Read simulation parameters from the file parameters.in
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine read_parameters_in()
use global_variables

implicit none

character(len=80) :: filename = 'parameters.in' !< name of the file in which parameters are stored
character(len=80) :: line
character(len=1), parameter :: comment_character = '!' !< character that will indicate that the rest of the line is a comment
integer :: comment_position !< the index of the comment character on the line. if zero, there is none on the current string
integer :: error !< to store the state of a read instruction

logical :: isParameter, isDefined
character(len=80) :: identificator, value
!------------------------------------------------------------------------------

inquire(file=filename, exist=isDefined)
if (isDefined) then

  open(10, file=filename, status='old')
  
  do
    read(10, '(a)', iostat=error) line
    if (error /= 0) exit
      
    ! We get only what is on the left of an eventual comment parameter
      comment_position = index(line, comment_character)
    
    ! if there are comments on the current line, we get rid of them
    if (comment_position.ne.0) then
      line = line(1:comment_position - 1)
    end if
    
    call get_parameter_value(line, isParameter, identificator, value)
      
    if (isParameter) then
      select case(identificator)
      ! Solver
      case('relative_tolerance', 'RTOL') ! The old name is kept for compatibility reasons
        read(value, '(e12.6)') RELATIVE_TOLERANCE
      
      !Switches      
      case('is_3_phase')
        read(value, '(i2)') is_3_phase

      !case('is_dust_1D')
      !  read(value, '(i2)') is_dust_1D

    case('photo_disk')
        read(value, '(i2)') photo_disk

      !case('nb_grains')
      !  read(value, '(i4)') nb_grains

      case('multi_grain')
        read(value, '(i2)') multi_grain

      case('is_grain_reactions', 'IDUST') ! The old name is kept for compatibility reasons
        read(value, '(i2)') IS_GRAIN_REACTIONS

      case('is_h2_adhoc_form')
        read(value, '(i2)') IS_H2_ADHOC_FORM

      case('is_h2_formation_rate')
        read(value, '(i2)') is_h2_formation_rate

      case('height_h2formation')
        read(value, '(i2)') height_h2formation
      
      case('preliminary_test')
        read(value, '(i2)') IS_TEST
      
      case('is_absorption_h2')
        read(value, '(i2)') is_absorption_h2

      case('is_absorption_co')
        read(value, '(i2)') is_absorption_co

      case('is_absorption_n2')
        read(value, '(i2)') is_absorption_n2

      case('is_photodesorb')
      read(value, '(i2)') is_photodesorb

      case('is_crid')
      read(value, '(i2)') is_crid

      case('is_er_cir')
      read(value, '(i2)') is_er_cir

      case('grain_tunneling_diffusion', 'IGRQM') ! The old name is kept for compatibility reasons
        read(value, '(i2)') GRAIN_TUNNELING_DIFFUSION
      
      case('modify_rate_flag', 'IMODH') ! The old name is kept for compatibility reasons
        read(value, '(i2)') MODIFY_RATE_FLAG
      
      case('conservation_type', 'ICONS') ! The old name is kept for compatibility reasons
        read(value, '(i2)') CONSERVATION_TYPE
        
      case('structure_type')
        read(value, *) STRUCTURE_TYPE
      
      ! 1D definitions
      case('spatial_resolution')
        read(value, '(i4)') spatial_resolution
              
      ! Gas phase
      case('initial_gas_density', 'XNT0') ! The old name is kept for compatibility reasons
        read(value, '(e12.6)') initial_gas_density
      
      case('initial_gas_temperature', 'TEMP0') ! The old name is kept for compatibility reasons
        read(value, '(e12.6)') initial_gas_temperature
      
      case('initial_visual_extinction', 'TAU0') ! The old name is kept for compatibility reasons
        read(value, '(e12.6)') INITIAL_VISUAL_EXTINCTION
      
      case('cr_ionisation_rate', 'ZETA0') ! The old name is kept for compatibility reasons
        read(value, '(e12.6)') CR_IONISATION_RATE
      
      case('x_ionisation_rate', 'ZETAX') ! The old name is kept for compatibility reasons
        read(value, '(e12.6)') X_IONISATION_RATE
      
      case('uv_flux', 'UVGAS') ! The old name is kept for compatibility reasons
        read(value, '(e12.6)') UV_FLUX
      
      ! Grain
      case('initial_dust_temperature', 'DTEMP0') ! The old name is kept for compatibility reasons
        read(value, '(e12.6)') initial_dust_temperature
      
      case('initial_dtg_mass_ratio', 'DTOGM') ! The old name is kept for compatibility reasons
        read(value, '(e12.6)') initial_dtg_mass_ratio
      
      case('sticking_coeff_neutral', 'STICK0') ! The old name is kept for compatibility reasons
        read(value, '(e12.6)') sticking_coeff_neutral
      
      case('sticking_coeff_positive', 'STICKP') ! The old name is kept for compatibility reasons
        read(value, '(e12.6)') sticking_coeff_positive
      
      case('sticking_coeff_negative', 'STICKN') ! The old name is kept for compatibility reasons
        read(value, '(e12.6)') sticking_coeff_negative
      
      case('grain_density', 'RHOD') ! The old name is kept for compatibility reasons
        read(value, '(e12.6)') GRAIN_DENSITY
      
      case('grain_radius', 'RD') ! The old name is kept for compatibility reasons
        read(value, '(e12.6)') GRAIN_RADIUS
        
      case('diffusion_barrier_thickness', 'ACM') ! The old name is kept for compatibility reasons
        read(value, '(e12.6)') DIFFUSION_BARRIER_THICKNESS
      
      case('surface_site_density', 'SNS') ! The old name is kept for compatibility reasons
        read(value, '(e12.6)') SURFACE_SITE_DENSITY
        
      case('diff_binding_ratio_surf')
        read(value, '(e12.6)') DIFF_BINDING_RATIO_SURF
      
      case('diff_binding_ratio_mant')
        read(value, '(e12.6)') DIFF_BINDING_RATIO_MANT
      
      case('chemical_barrier_thickness', 'ACT') ! The old name is kept for compatibility reasons
        read(value, '(e12.6)') CHEMICAL_BARRIER_THICKNESS
      
      case('cr_peak_grain_temp', 'TSMAX') ! The old name is kept for compatibility reasons
        read(value, '(e12.6)') CR_PEAK_GRAIN_TEMP
      
      case('cr_peak_duration', 'CRT') ! The old name is kept for compatibility reasons
        read(value, '(e12.6)') CR_PEAK_DURATION
      
      case('Fe_ionisation_rate', 'CRFE') ! The old name is kept for compatibility reasons
        read(value, '(e12.6)') FE_IONISATION_RATE
      
      case('vib_to_dissip_freq_ratio', 'ARRK') ! The old name is kept for compatibility reasons
        read(value, '(e12.6)') VIB_TO_DISSIP_FREQ_RATIO

      case('ED_H2')
      read(value, '(e12.6)') ED_H2
    
      !case('is_dust_MRN')
      !ead(value, '(i2)') is_dust_MRN
     
      case('nb_active_lay')
      read(value, '(e12.6)') nb_active_lay
      
      ! Outputs
      case('nb_outputs')
        read(value, '(i4)') NB_OUTPUTS
      
      case('start_time')
        read(value, '(e12.6)') START_TIME
      
      case('stop_time')
        read(value, '(e12.6)') STOP_TIME
      
      case('output_type')
        read(value, *) OUTPUT_TYPE
      
      ! Initial abundances
      case('minimum_initial_abundance')
        read(value, '(e12.6)') MINIMUM_INITIAL_ABUNDANCE
      
      case('is_structure_evolution')
        read(value, '(i2)') IS_STRUCTURE_EVOLUTION
      
      case('grain_temperature_type')
        read(value, *) GRAIN_TEMPERATURE_TYPE
         
      case default
        write(*,*) 'Warning: An unknown parameter has been found'
        write(*,*) "identificator='", trim(identificator), "' ; value(s)='", trim(value),"'"
      end select
    end if
  end do
  close(10)
  
else
  write (*,*) 'Warning: The file "parameters.in" does not exist. Default values have been used'
end if

START_TIME = START_TIME * YEAR
STOP_TIME = STOP_TIME * YEAR

if ((STRUCTURE_TYPE.eq.'0D').and.(spatial_resolution.ne.1)) then
  write(Error_unit,'(a,i0,a)') 'Error: In 0D, we must have one point (spatial_resolution=',spatial_resolution,')'
  call exit(22)
endif

if ((STRUCTURE_TYPE.ne.'0D').and.(spatial_resolution.eq.1)) then
  write(Error_unit,'(a,i0,a)') 'Error: In 1D, we must have more than one point (spatial_resolution=',spatial_resolution, ')'
  call exit(22)
endif

if ((STRUCTURE_TYPE.ne.'0D').and.(IS_STRUCTURE_EVOLUTION.ne.0)) then
  write(Error_unit,*) 'Error: In 1D, structure evolution is currently not supported. If is_structure_evolution=1, you must have'
  write(Error_unit,*) '       structure_type="0D"'
  
  call exit(23)
endif

return
end subroutine read_parameters_in


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Christophe Cossou
!
!> @date 2014
!
! DESCRIPTION: 
!> @brief Read all species information from gas_species.in and grain_species.in
!! files.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
!
!  modified by Wasim Iqbal on 23-02-2017 to add grain size distribution
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

subroutine read_species()

use global_variables
use iso_fortran_env

implicit none

! Locals
integer           :: i,j,k
character(2)      :: i_c !character variable to convert integer value of j to character j
character(len=80) :: filename !< name of the file to be read
character(len=80) :: line
character(len=1), parameter :: comment_character = '!' !< character that will indicate that the rest of the line is a comment
integer :: comment_position !< the index of the comment character on the line. if zero, there is none on the current string
integer :: error !< to store the state of a read instruction
logical :: isDefined

! Variables for the unordered reaction file
character(len=11), dimension(nb_species_for_gas)           :: gas_species_label 
integer, dimension(nb_species_for_gas)                     :: ICG1 
integer, dimension(NB_PRIME_ELEMENTS, nb_species_for_gas)  :: gas_species_composition 

character(len=11), dimension(nb_species_for_grain)         :: surface_species_label 
integer, dimension(nb_species_for_grain)                   :: ICG2 
integer, dimension(NB_PRIME_ELEMENTS, nb_species_for_grain):: grain_species_composition 


! Reading list of species for gas phase
filename = 'gas_species.in'
inquire(file=filename, exist=isDefined)

if (isDefined) then
  open(10, file=filename, status='old')
  i = 0
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
      i = i + 1
      read(line, '(A11,i3,13(I3))')  gas_species_label(I),ICG1(I),(gas_species_composition(K,I),K=1,NB_PRIME_ELEMENTS) 
      
      if ((gas_species_label(i)(:5) == 'GRAIN')) then
        if ((gas_species_label(i) == 'GRAIN0     ')) then
          gas_species_label(i)='GRAIN01    '           ! channging old name to new
          do j=1,nb_grains-1
            i=i+1
            write(i_c,'(I2.2)')j+1
            gas_species_label(i)=gas_species_label(i-j)(:5)//trim(i_c)//gas_species_label(i-j)(8:)
            ICG1(I)=  ICG1(I-j)
            do k=1,NB_PRIME_ELEMENTS
              gas_species_composition(K,I)=gas_species_composition(K,I-j)
            enddo  
          enddo   
        elseif ((gas_species_label(i) == 'GRAIN-     ')) then 
          gas_species_label(i)='GRAIN01-   '
          do j=1,nb_grains-1
            i=i+1
            write(i_c,'(I2.2)')j+1
            gas_species_label(i)=gas_species_label(i-j)(:5)//trim(i_c)//gas_species_label(i-j)(8:)
            ICG1(I)=  ICG1(I-j)
            do k=1,NB_PRIME_ELEMENTS
              gas_species_composition(K,I)=gas_species_composition(K,I-j)
            enddo  
          enddo  
        endif
      endif        
    end if
  end do
  close(10)
else
  write(Error_unit,*) 'Error: The file ', trim(filename),' does not exist.'
  call exit(1)
end if


! Reading list of species for grain surface
filename = 'grain_species.in'
inquire(file=filename, exist=isDefined)

if (isDefined) then
  open(10, file=filename, status='old')
  i = 0
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
      i = i + 1
      read(line, '(A11,i3,13(I3))')  surface_species_label(I),ICG2(I),(grain_species_composition(K,I),K=1,NB_PRIME_ELEMENTS) 
      
!       write(*,*)surface_species_label(I),ICG2(I),(grain_species_composition(K,I),K=1,NB_PRIME_ELEMENTS) 
      
      if ((surface_species_label(i)(1:1).eq.'J').or.(surface_species_label(i)(1:1).eq.'K')) then
        write(i_c,'(I2.2)')1
        surface_species_label(i)=surface_species_label(i)(:1)//trim(i_c)//surface_species_label(i)(2:)

        do j=1,nb_grains-1
          i=i+1
          write(i_c,'(I2.2)')j+1
          surface_species_label(i)=surface_species_label(i-j)(:1)//trim(i_c)//surface_species_label(i-j)(4:)
          ICG2(I)=  ICG2(I-j)
          do k=1,NB_PRIME_ELEMENTS
            grain_species_composition(K,I)=grain_species_composition(K,I-j)
          enddo  
          
!           write(*,*)surface_species_label(I),ICG2(I),(grain_species_composition(K,I),K=1,NB_PRIME_ELEMENTS) 
        enddo   
      endif    
    end if
  enddo
!   do i=1,nb_species_for_grain
!     write(*,*)surface_species_label(I),ICG2(I),(grain_species_composition(K,I),K=1,NB_PRIME_ELEMENTS)
!   enddo
!      pause
  close(10)
  
else
  write(Error_unit,*) 'Error: The file ', trim(filename),' does not exist.'
  call exit(1)
end if


! putting everything back into the big tables
do I=1,nb_species_for_gas 
  species_name(I)=gas_species_label(I)
  SPECIES_CHARGE(I)=ICG1(I)
  do k=1,NB_PRIME_ELEMENTS
    species_composition(K,I)=gas_species_composition(K,I)
  enddo
enddo
do I=1,nb_species_for_grain 
  species_name(nb_species_for_gas+I)=surface_species_label(I)
  SPECIES_CHARGE(nb_species_for_gas+I)=ICG2(I)
  do k=1,NB_PRIME_ELEMENTS
    species_composition(K,nb_species_for_gas+I)=grain_species_composition(K,I)
  enddo
enddo

!  write(*,*)"--------------------------------------------------------------------------------------"
!     do i=1,nb_species
!     write(9990,*)species_name(I),SPECIES_CHARGE(I),(species_composition(K,I),k=1,NB_PRIME_ELEMENTS)
!      enddo
!      pause

return
end subroutine read_species


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Franck Hersant
!
!> @date 2000
!
! DESCRIPTION: 
!> @brief Read all reactions both for gas phase and grain surface and order them
! by ITYPE
! TODO: explain the reordering process
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
!
!  modified by Wasim Iqbal on 23-02-2017 to add grain size distribution
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine read_reactions()

use global_variables

implicit none

! Locals
integer :: i,j,k,k1,k2,jk
integer :: NUM1start,NUM1end,NUM2start,NUM2end ! integer to count increase in reaction number due to extra grains and accordingly assign NUM2 values 
integer :: oldid ! integer used to save id of reaction
logical :: isidsame
character(2)      :: i_c !character variable to convert integer value of j to character j
character(len=80) :: filename !< name of the file to be read
character(len=200):: line
character(len=2)  :: c_i   !variable to read grain number, i.e. to read a given reaction is for which grain_radii
character(len=1), parameter :: comment_character = '!' !< character that will indicate that the rest of the line is a comment
integer :: comment_position !< the index of the comment character on the line. if zero, there is none on the current string
integer :: error !< to store the state of a read instruction
logical :: isDefined

! Variables for the unordered reaction file
character(len=11), dimension(MAX_COMPOUNDS,nb_gas_phase_reactions) :: SYMBOLUO1 
real(double_precision), dimension(nb_gas_phase_reactions)          :: AUO1,BUO1,CUO1 
integer, dimension(nb_gas_phase_reactions) :: itypeUO1,Tmin1,Tmax1,FORMULA1,NUM1 

character (len=11), dimension(MAX_COMPOUNDS,nb_surface_reactions)  :: SYMBOLUO2 
real(double_precision), dimension(nb_surface_reactions)            :: AUO2,BUO2,CUO2 
integer, dimension(nb_surface_reactions)   :: itypeUO2,Tmin2,Tmax2,FORMULA2,NUM2 

character (len=11), dimension(MAX_COMPOUNDS,nb_reactions) :: SYMBOLUO
real(double_precision), dimension(nb_reactions)           :: AUO,BUO,CUO
integer, dimension(nb_reactions) :: itypeUO,TminUO,TmaxUO,FORMULAUO,NUMUO

character(len=80) :: line_format !< format of one line of gas_reaction.in or grain_reaction.in

! Definition of the line format, common to gas_reaction.in and grain_reaction.in
write(line_format, '(a,i0,a,i0,a)') '(', MAX_REACTANTS, 'A11,1x,', MAX_PRODUCTS, 'A11,3D11.3,23x,I3,2i7,i3,i6)'

! Reading list of reaction for gas phase
filename = 'gas_reactions.in'
inquire(file=filename, exist=isDefined)

if (isDefined) then
  open(10, file=filename, status='old')
  j = 0
  NUM1start = 0
  NUM1end   = NUM1start 
  oldid     = 0
  isidsame  = .FALSE.
  
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
      j = j + 1
      read(line, line_format) (SYMBOLUO1(I,J),I=1,MAX_COMPOUNDS),AUO1(J),BUO1(J),CUO1(J), &
                                  ITYPEUO1(J),Tmin1(j),Tmax1(j),FORMULA1(J),NUM1(J)

      if ((SYMBOLUO1(1,j)(1:5).eq.'GRAIN').or. (SYMBOLUO1(2,j)(1:5).eq.'GRAIN').or.&
          (SYMBOLUO1(3,j)(1:5).eq.'GRAIN').or. (SYMBOLUO1(4,j)(1:5).eq.'GRAIN').or.&
          (SYMBOLUO1(5,j)(1:5).eq.'GRAIN').or. (SYMBOLUO1(6,j)(1:5).eq.'GRAIN').or.&
          (SYMBOLUO1(7,j)(1:5).eq.'GRAIN').or. (SYMBOLUO1(8,j)(1:5).eq.'GRAIN')) then
        if (NUM1(j) == oldid) then
          NUM1(J)= NUM1start
          isidsame = .TRUE.
        else
          isidsame = .FALSE.
          NUM1end  = NUM1end+1
          NUM1start= NUM1end
          oldid    = NUM1(J)
          NUM1(J)  = NUM1start
        endif  
        
        
        do k2=1,8
          if(SYMBOLUO1(k2,j) == 'GRAIN0     ') then
            SYMBOLUO1(k2,j) = 'GRAIN01    '
          elseif (SYMBOLUO1(k2,j) == 'GRAIN-     ') then
            SYMBOLUO1(k2,j) = 'GRAIN01-   '
          endif             
        enddo
        
        do k1=1,nb_grains-1  
          j=j+1
          write(i_c,'(I2.2)')k1+1
          do k2=1,8
            if(SYMBOLUO1(k2,j-k1) == 'GRAIN01    ') then
              SYMBOLUO1(k2,j)=SYMBOLUO1(k2,j-k1)(:5)//trim(i_c)//SYMBOLUO1(k2,j-k1)(8:)
            elseif (SYMBOLUO1(k2,j-k1) == 'GRAIN01-   ') then
              SYMBOLUO1(k2,j)=SYMBOLUO1(k2,j-k1)(:5)//trim(i_c)//SYMBOLUO1(k2,j-k1)(8:) 
            else
              SYMBOLUO1(k2,j)=SYMBOLUO1(k2,j-k1)
            endif             
          enddo
          AUO1(J)=AUO1(J-k1)
          BUO1(J)=BUO1(J-k1)
          CUO1(J)=CUO1(J-k1)
          ITYPEUO1(J)=ITYPEUO1(J-k1)
          Tmin1(J)=  Tmin1(J-k1)
          Tmax1(J)=Tmax1(J-k1)
          FORMULA1(J)=FORMULA1(J-k1)
          if (isidsame) then       
            NUM1end  = NUM1start + k1
            NUM1(J)  = NUM1end 
          else
            NUM1end  = NUM1end+1
            NUM1(J)  = NUM1end 
          endif  
        enddo
      else
        if(NUM1start /= NUM1end) then
          NUM1end  = NUM1end+1
          NUM1start= NUM1end
          oldid    = NUM1(J)
          NUM1(J)  = NUM1end
        else
          if (NUM1(j) == oldid) then
            NUM1(J)= NUM1end          
          else  
            NUM1end  = NUM1end+1
            NUM1start= NUM1end
            oldid    = NUM1(J)
            NUM1(J)  = NUM1end          
          endif  
        endif         
      endif
    end if
  end do
  close(10)
  
else
  write(Error_unit,*) 'Error: The file ', trim(filename),' does not exist.'
  call exit(1)
end if

! Reading list of reaction for grains
NUM2start= NUM1end
NUM2end  = NUM1end
oldid    = 0
isidsame = .FALSE.

filename = 'grain_reactions.in'
inquire(file=filename, exist=isDefined)
if (isDefined) then

  open(10, file=filename, status='old')
  
  j = 0
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
      j = j + 1
      read(line, line_format)  (SYMBOLUO2(I,J),I=1,MAX_COMPOUNDS),AUO2(J),BUO2(J),CUO2(J), &
                                   ITYPEUO2(J),Tmin2(j),Tmax2(j),FORMULA2(J),NUM2(J)

      if ((SYMBOLUO2(1,j)(1:1).eq.'J').or.(SYMBOLUO2(1,j)(1:1).eq.'K') .or. &
          (SYMBOLUO2(2,j)(1:1).eq.'J').or.(SYMBOLUO2(2,j)(1:1).eq.'K') .or. &
          (SYMBOLUO2(3,j)(1:1).eq.'J').or.(SYMBOLUO2(3,j)(1:1).eq.'K') .or. &
          (SYMBOLUO2(4,j)(1:1).eq.'J').or.(SYMBOLUO2(4,j)(1:1).eq.'K') .or. &
          (SYMBOLUO2(5,j)(1:1).eq.'J').or.(SYMBOLUO2(5,j)(1:1).eq.'K') .or. &
          (SYMBOLUO2(6,j)(1:1).eq.'J').or.(SYMBOLUO2(6,j)(1:1).eq.'K') .or. &
          (SYMBOLUO2(7,j)(1:1).eq.'J').or.(SYMBOLUO2(7,j)(1:1).eq.'K') .or. &
          (SYMBOLUO2(8,j)(1:1).eq.'J').or.(SYMBOLUO2(8,j)(1:1).eq.'K') ) then
        write(i_c,'(I2.2)')1
        do k2=1,8
            if(SYMBOLUO2(k2,j)(:1)=='J' .or. SYMBOLUO2(k2,j)(:1)=='K') then
              SYMBOLUO2(k2,j)=SYMBOLUO2(k2,j)(:1)//trim(i_c)//SYMBOLUO2(k2,j)(2:)
            endif             
        enddo
                
        if (NUM2(j) == oldid) then
          NUM2(J)= NUM2start
          isidsame = .TRUE.
        else
          isidsame = .FALSE.
          NUM2end  = NUM2end+1
          NUM2start= NUM2end
          oldid    = NUM2(J)
          NUM2(J)  = NUM2start
        endif          
        do k1=1,nb_grains-1  
          j=j+1
          write(i_c,'(I2.2)')k1+1
          do k2=1,8
            if(SYMBOLUO2(k2,j-k1)(:1)=='J' .or. SYMBOLUO2(k2,j-k1)(:1)=='K') then
              SYMBOLUO2(k2,j)=SYMBOLUO2(k2,j-k1)(:1)//trim(i_c)//SYMBOLUO2(k2,j-k1)(4:)
            else
              SYMBOLUO2(k2,j)=SYMBOLUO2(k2,j-k1)
            endif             
          enddo
          AUO2(J)=AUO2(J-k1)
          BUO2(J)=BUO2(J-k1)
          CUO2(J)=CUO2(J-k1)
          ITYPEUO2(J)=ITYPEUO2(J-k1)
          Tmin2(J)=  Tmin2(J-k1)
          Tmax2(J)=Tmax2(J-k1)
          FORMULA2(J)=FORMULA2(J-k1)
          if (isidsame) then       
            NUM2end  = NUM2start + k1
            NUM2(J)  = NUM2end 
          else
            NUM2end  = NUM2end+1
            NUM2(J)  = NUM2end 
          endif  
        enddo
      else
        if(NUM2start /= NUM2end) then
          NUM2end  = NUM2end+1
          NUM2start= NUM2end
          oldid    = NUM2(J)
          NUM2(J)  = NUM2end
        else
          if (NUM2(j) == oldid) then
            NUM2(J)= NUM2end          
          else  
            NUM2end  = NUM2end+1
            NUM2start= NUM2end
            oldid    = NUM2(J)
            NUM2(J)  = NUM2end          
          endif  
        endif  
      endif  
    end if
  end do
  close(10)
  
else
  write(Error_unit,*) 'Error: The file ', trim(filename),' does not exist.'
  call exit(1)
end if

! putting everything back into the big tables
do I=1,nb_gas_phase_reactions 
  do k=1,MAX_COMPOUNDS
    SYMBOLUO(k,I)=SYMBOLUO1(k,I)
  enddo
  AUO(I)=AUO1(I)
  BUO(I)=BUO1(I)
  CUO(I)=CUO1(I)
  ITYPEUO(I)=ITYPEUO1(I)
  TminUO(I) = Tmin1(I)
  TmaxUO(I) = Tmax1(I)
  FORMULAUO(I) = FORMULA1(I)
  NUMUO(I) = NUM1(I)
enddo

do I=1,nb_surface_reactions 
  do k=1,MAX_COMPOUNDS
    SYMBOLUO(k,nb_gas_phase_reactions+I)=SYMBOLUO2(k,I)
  enddo
  AUO(nb_gas_phase_reactions+I)=AUO2(I)
  BUO(nb_gas_phase_reactions+I)=BUO2(I)
  CUO(nb_gas_phase_reactions+I)=CUO2(I)
  ITYPEUO(nb_gas_phase_reactions+I)=ITYPEUO2(I)
  TminUO(nb_gas_phase_reactions+I) = Tmin2(I)
  TmaxUO(nb_gas_phase_reactions+I) = Tmax2(I)
  FORMULAUO(nb_gas_phase_reactions+I) = FORMULA2(I)
  NUMUO(nb_gas_phase_reactions+I) = NUM2(I)
enddo

! write(*,*)"--------------------------------------------------------------------------------------"
!    do j=1, nb_reactions
!      write(9,line_format)(SYMBOLUO(I,J),I=1,MAX_COMPOUNDS),AUO(J),BUO(J),CUO(J), &
!                                       ITYPEUO(J),TminUO(j),TmaxUO(j),FORMULAUO(J),NUMUO(J)
!    enddo
!  write(*,*)"--------------------------------------------------------------------------------------"
!    pause

! Reorder reaction file entries with ITYPE
jk=1
do i=0,MAX_NUMBER_REACTION_TYPE
  do j=1,nb_reactions
    if (itypeuo(j).eq.i) then
      REACTION_COMPOUNDS_NAMES(:,jk)=SYMBOLUO(:,j)     
      RATE_A(jk)=AUO(j)
      RATE_B(jk)=BUO(j)
      RATE_C(jk)=CUO(j)
      REACTION_TYPE(jk)=itypeuo(j)
      REACTION_TMIN(jk) = dble(TminUO(j))
      REACTION_TMAX(jk) = dble(TmaxUO(j))
      RATE_FORMULA(jk) = FORMULAUO(J)
      REACTION_ID(jk) = NUMUO(j)
      jk=jk+1
    endif
  enddo
enddo

! write(*,*)"--------------------------------------------------------------------------------------"
!   do j=1, nb_reactions
!     write(*,line_format)(REACTION_COMPOUNDS_NAMES(I,J),I=1,MAX_COMPOUNDS),RATE_A(J),RATE_B(J),RATE_C(J), &
!                                      REACTION_TYPE(J),TminUO(j),TmaxUO(j),RATE_FORMULA(J),REACTION_ID(J)
!   enddo
!  write(*,*)"--------------------------------------------------------------------------------------"
!   pause

! Get grain size or grain rank for each reordered reactions

do i=1,nb_reactions
  do j=1,MAX_COMPOUNDS
    if (REACTION_COMPOUNDS_NAMES(j,i)(:1) == 'J' .or. REACTION_COMPOUNDS_NAMES(j,i)(:1) == 'K' ) then
      c_i=REACTION_COMPOUNDS_NAMES(j,i)(2:3)
      read(c_i,'(I2)')GRAIN_RANK(i)
!             write(*,*)c_i

    elseif (REACTION_COMPOUNDS_NAMES(j,i)(:5) == 'GRAIN') then
      c_i=REACTION_COMPOUNDS_NAMES(j,i)(6:7)
!       write(*,*)c_i
      read(c_i,'(I2)')GRAIN_RANK(i)
    endif
  enddo 
!  if (i < 20000)  then
!  write(798,*)i,GRAIN_RANK(i),grain_radii(GRAIN_RANK(i)),grain_radii(GRAIN_RANK(i)-1),grain_radii(GRAIN_RANK(i)+1),&
!              REACTION_COMPOUNDS_NAMES(:,i)
!  endif
enddo  
!  pause
if (jk.ne.nb_reactions+1) then
  write(Error_unit,*) 'Some reaction was not found by the reorder process'
  write(Error_unit,*) jk,'=/',nb_reactions+1 
  call exit(4)
endif

!       replace the species names by blanks for non chemical species                                                                        
do j=1,nb_reactions-1
  do i=1,MAX_COMPOUNDS
    select case(REACTION_COMPOUNDS_NAMES(i,j))
      case ('CR', 'CRP', 'Photon')
        REACTION_COMPOUNDS_NAMES(i,j) = '           '
    end select
  enddo
enddo 

! pause
return
end subroutine read_reactions


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Christophe Cossou
!
!> @date 2014
!
! DESCRIPTION: 
!> @brief Read abundances from abundances.in file. All abundances not defined
!! here will have the default value MINIMUM_INITIAL_ABUNDANCE
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine read_abundances()
! Writes 1D outputs
use global_variables

implicit none

! Locals
character(len=80), parameter :: filename='abundances.in'
integer :: i, j

real(double_precision), allocatable, dimension(:) :: temp_abundances
character(len=11), allocatable, dimension(:) :: temp_names
integer :: nb_lines

character(len=80) :: line
character(len=1), parameter :: comment_character = '!' ! character that will indicate that the rest of the line is a comment
integer :: comment_position ! the index of the comment character on the line. if zero, there is none on the current string
integer :: error ! to store the state of a read instruction

logical :: isParameter, isDefined
character(len=80) :: identificator, value

call get_linenumber(filename, nb_lines)

allocate(temp_abundances(nb_lines))
allocate(temp_names(nb_lines))

temp_abundances(1:nb_lines) = 0.d0
temp_names(1:nb_lines) = ''

  !------------------------------------------------------------------------------
  
inquire(file=filename, exist=isDefined)
if (isDefined) then

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
    
    call get_parameter_value(line, isParameter, identificator, value)
      
    if (isParameter) then
      read(value, '(e12.6)') temp_abundances(I)
      read(identificator, *) temp_names(I)
      i = i + 1
    end if
  enddo
  close(10)
endif

! We check if all species in abundances.in exists in the simulation
do j=1,nb_lines
  error = 1
  do i=1,nb_species
    if (temp_names(j).eq.species_name(i)) then
      error = 0 ! The species exist
    endif
  enddo
  
  if (error.eq.1) then
    write(*,*) j
    write(*,*) temp_names
    write(Error_Unit,*) 'Input species "', trim(temp_names(j)), '" in "', trim(filename), '" do not match those in reaction file'
    call exit(2)
  endif
enddo

! Set initial abundances================================================
abundances(1:nb_species, 1:spatial_resolution) = MINIMUM_INITIAL_ABUNDANCE

! Initial abundance for one species is assumed to be the same throughout the 1D structure initially
do i=1,nb_species
  do j=1,nb_lines
    if ((species_name(i).EQ.temp_names(j)).AND.(temp_abundances(j).NE.0.D0)) then
      abundances(i,1:spatial_resolution)=temp_abundances(j)
    endif
  enddo
enddo

deallocate(temp_abundances)
deallocate(temp_names)

return
end subroutine read_abundances


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Christophe Cossou
!
!> @date 2014
!
! DESCRIPTION: 
!> @brief From the label of each species, determine if this is a gas
!! phase or a surface species. Set the values of nb_gaseous_species and
!! nb_surface_species global parameters that count the total number of species
!! in each category.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine get_gas_surface_species()

use global_variables

implicit none

! Locals
integer :: i

! We retrieve the total number of gas and surface species (not the ones that are involved in reactions, but the actual position of the species)
nb_gaseous_species = 0
nb_surface_species = 0
do i=1,nb_species
  if ((species_name(i)(1:1).eq.'J').or.(species_name(i)(1:1).eq.'K')) then
    nb_surface_species = nb_surface_species + 1
  else
    nb_gaseous_species = nb_gaseous_species + 1
  endif
end do
!  write(*,*)nb_surface_species , nb_gaseous_species
end subroutine get_gas_surface_species



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Wasim Iqbal
!
!> @date 2017
!
! DESCRIPTION: 
!> @brief calculates mass of different grains sizes
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

subroutine get_grain_mass() 
use global_variables

implicit none
character(2)      :: i_c !character variable to convert integer value of j to character j
character(len=11) :: temp_name1,temp_name2
integer           :: i,j

  ! ------ Calculate mass of grains

temp_name1 =''
temp_name2 =''
do i=1,nb_species 
  if (species_name(I)(:5).EQ.'GRAIN') then
    do j=1,nb_grains
      write(i_c,'(I2.2)')j
      temp_name1 ="GRAIN"//trim(i_c)
      temp_name2 ="GRAIN"//trim(i_c)//"-"
      if(species_name(I).EQ. temp_name1 .or. species_name(I).EQ.temp_name2) then 
        SPECIES_MASS(I)=4.0*PI*grain_radii(j)*grain_radii(j)*grain_radii(j)*GRAIN_DENSITY/3.0/AMU
!            write(*,'(I0,2ES15.6)')i,SPECIES_MASS(I),grain_radii(j)
      endif  
    enddo  
  endif  
enddo   
      
!  pause
return
end subroutine get_grain_mass



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Wasim Iqbal
!
!> @date 2017
!
! DESCRIPTION: 
!> @brief calculates dust to gas ratio using MRN distribution
!!   (1/n_H)*(dn_gr/dr)= C * r^(-3.5) for r_min < r < r_max, 
!!!    where, n_H is umber density of H, n_gr is number density of grain,
!!!!   C is the grain constant and r is the grain radius. 
!!!!!  This relation is valid between r_min = 50 Å (5 nm) and r_max = 0.25 μm.
!!!!!! The grain constant C is 10^(−25.13) for carbon and 10^(−25.11) cm^2.5 for silicates.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine get_MRN_distribution()
use global_variables
implicit none

real(double_precision), parameter :: r_min= 5.d-7 , r_max=1.d-4  ! in cm ....minimum and maximum values of radii in MRN model  !5.d-7 3.11914d-6 arbitrary value
real(double_precision), parameter :: const_C=-25.11d0 !The grain constant C is 10^(−25.13) for carbon and 10^(−25.11) cm^2.5 for silicates
integer          :: i
real(double_precision)           :: a0,a1
   do i=1,nb_grains
    if (i==1) then
      a0=r_min
    else  
      a0=(grain_radii(i) + grain_radii(i-1))/2.0
    endif
    if (i==nb_grains) then
      a1=r_max
    else
      a1=(grain_radii(i) + grain_radii(i+1))/2.0
    endif
        GTODN(i)= ((a0**(-2.5)-a1**(-2.5)) * 10**const_C)/2.5d0
        GTODN(i)= 1.d0/GTODN(i)
   enddo
end subroutine get_MRN_distribution

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Wasim Iqbal
!
!> @date 2017
!
! DESCRIPTION: 
!> @brief calculates dust to gas ratio using WD distribution
!!  This subroutine uses subroutine GRAIN_DIST, I have put 
!!! subroutine GRAIN_DIST just below this subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
SUBROUTINE  get_WD_distribution()
use global_variables
IMPLICIT NONE
INTEGER             :: INDEX,i,dtype
real(double_precision)             :: a0,a1,da

INDEX =16
dtype = 1
 
DO i=1,nb_grains
  if (i==1 ) then
    a0=grain_radii(i)-(grain_radii(i+1)-grain_radii(i))/2.d0
  ELSE
    a0= grain_radii(i-1)
  endif
  a1=grain_radii(i)
  da=a1-a0
  CALL GRAIN_DIST(INDEX,DTYPE,grain_radii(i),GTODN(i))
  GTODN(i)= 1.d0/(GTODN(i)*da)
ENDDO
END SUBROUTINE  get_WD_distribution



SUBROUTINE GRAIN_DIST(INDEX,DTYPE,A,DNDA)
use global_variables
IMPLICIT NONE

!C  History:  Feb 12, 2002:  Corrected bug--previous version did not
!C  multiply
!C            very small carbonaceous grain size dist by b_C.

!C  Output:   DNDA = (1/n_H) dn_gr/da, where n_gr(a) is the number
!C  density
!C                   of grains with radius < a and n_H is the H nucleus
!C                   number density   (cm^-1)
!C  Input:   DTYPE = 1 for silicate dust  
!C                   2 for carbonaceous dust
!C               A = grain radius (cm)

!C            INDEX  is an integer indicating which set of conditions to 
!C            adopt

!C          INDEX    R_V    10^5 b_C   case  

!C            1      3.1      0.0       A        
!C            2      3.1      1.0       A       
!C            3      3.1      2.0       A        
!C            4      3.1      3.0       A        
!C            5      3.1      4.0       A        
!C            6      3.1      5.0       A       
!C            7      3.1      6.0       A        
!C            8      4.0      0.0       A
!C            9      4.0      1.0       A
!C            10     4.0      2.0       A
!C            11     4.0      3.0       A
!C            12     4.0      4.0       A
!C            13     5.5      0.0       A
!C            14     5.5      1.0       A
!C            15     5.5      2.0       A
!C            16     5.5      3.0       A
!C            17     4.0      0.0       B
!C            18     4.0      1.0       B
!C            19     4.0      2.0       B
!C            20     4.0      3.0       B
!C            21     4.0      4.0       B
!C            22     5.5      0.0       B
!C            23     5.5      1.0       B
!C            24     5.5      2.0       B
!C            25     5.5      3.0       B

INTEGER ::INDEX,DTYPE
real(double_precision)  ::A,DNDA

real(double_precision) :: ALPHAGARR(25),BETAGARR(25),ATGARR(25),ACGARR(25),CGARR(25),&
&      ALPHASARR(25),BETASARR(25),ATSARR(25),CSARR(25),ACS,       &
&      BC5ARR(25),ALPHAG,BETAG,ATG,ACG,CG,ALPHAS,BETAS,ATS,CS,BC5,&
&      DNDAVSG,A01,A02,SIG,B1,B2  
 
DATA ALPHAGARR/&
 &-2.25,-2.17,-2.04,-1.91,-1.84,-1.72,-1.54,-2.26,-2.16,-2.01,      &
 &-1.83,-1.64,-2.35,-2.12,-1.94,-1.61,-2.62,-2.52,-2.36,-2.09,      &
 &-1.96,-2.80,-2.67,-2.45,-1.90/
DATA BETAGARR/&
 &-0.0648,-0.0382,-0.111,-0.125,-0.132,-0.322,-0.165,-0.199,-0.0862,&
 &-0.0973,-0.175,-0.247,-0.668,-0.67,-0.853,-0.722,-0.0144,-0.0541, &
 &-0.0957,-0.193,-0.813,0.0356,0.0129,-0.00132,-0.0517/
DATA ATGARR/&
 &0.00745,0.00373,0.00828,0.00837,0.00898,0.0254,.0107,0.0241,      &
 &0.00867,0.00811,0.0117,0.0152,0.148,0.0686,0.0786,0.0418,0.0187,  &
 &0.0366,0.0305,0.0199,0.0693,0.0203,0.0134,0.0275,0.012/
DATA ACGARR/&
 &0.606,0.586,0.543,0.499,0.489,0.438,0.428,0.861,0.803,0.696,      &
 &0.604,0.536,1.96,1.35,0.921,0.72,5.74,6.65,6.44,4.6,3.48,3.43,    &
 &3.44,5.14,7.28/
DATA CGARR/&
 &9.94e-11,3.79e-10,5.57e-11,4.15e-11,2.90e-11,3.20e-12,9.99e-12,   &
 &5.47e-12,4.58e-11,3.96e-11,1.42e-11,5.83e-12,4.82e-14,3.65e-13,   &
 &2.57e-13,7.58e-13,6.46e-12,1.08e-12,1.62e-12,4.21e-12,2.95e-13,   &
 &2.74e-12,7.25e-12,8.79e-13,2.86e-12/
DATA ALPHASARR/&
 &-1.48,-1.46,-1.43,-1.41,-2.1,-2.1,-2.21,-2.03,-2.05,-2.06,-2.08,  &
 &-2.09,-1.57,-1.57,-1.55,-1.59,-2.01,-2.11,-2.05,-2.1,-2.11,-1.09, &
 &-1.14,-1.08,-1.13/
DATA BETASARR/&
 &-9.34,-10.3,-11.7,-11.5,-0.114,-0.0407,0.3,0.668,0.832,0.995,     &
 &1.29,1.58,1.1,1.25,1.33,2.12,0.894,1.58,1.19,1.64,2.1,-0.37,      &
 &-0.195,-0.336,-0.109/
DATA ATSARR/&
 &0.172,0.174,0.173,0.171,0.169,0.166,0.164,0.189,0.188,0.185,      &
 &0.184,0.183,0.198,0.197,0.195,0.193,0.198,0.197,0.197,0.198,0.198,&
 &0.218,0.216,0.216,0.211/
DATA CSARR/&
 &1.02e-12,1.09e-12,1.27e-12,1.33e-12,1.26e-13,1.27e-13,1.e-13,     &
 &5.2e-14,4.81e-14,4.7e-14,4.26e-14,3.94e-14,4.24e-14,4.e-14,       &
 &4.05e-14,3.2e-14,4.95e-14,3.69e-14,4.37e-14,3.63e-14,3.13e-14,    &
 &1.17e-13,1.05e-13,1.17e-13,1.04e-13/
DATA BC5ARR/&
 &0.,1.,2.,3.,4.,5.,6.,0.,1.,2.,3.,4.,0.,1.,2.,3.,0.,1.,2.,3.,4.,   &
 &0.,1.,2.,3./

ALPHAG=ALPHAGARR(INDEX)
BETAG=BETAGARR(INDEX)
ATG=ATGARR(INDEX)*1.E-4
ACG=ACGARR(INDEX)*1.E-4
CG=CGARR(INDEX)
ALPHAS=ALPHASARR(INDEX)
BETAS=BETASARR(INDEX)
ATS=ATSARR(INDEX)*1.E-4
ACS=1.E-5
CS=CSARR(INDEX)
BC5=BC5ARR(INDEX)
IF (DTYPE .EQ. 1) THEN
  DNDA=(CS/A)*(A/ATS)**ALPHAS
  IF (BETAS .GE. 0.) THEN
   DNDA=DNDA*(1.+BETAS*A/ATS)
  ELSE
   DNDA=DNDA/(1.-BETAS*A/ATS)
  ENDIF
  IF (A .GT. ATS) DNDA=DNDA*EXP(((ATS-A)/ACS)**3)
ENDIF
IF (DTYPE .EQ. 2) THEN
  DNDA=(CG/A)*(A/ATG)**ALPHAG
  IF (BETAG .GE. 0.) THEN
    DNDA=DNDA*(1.+BETAG*A/ATG)
  ELSE
    DNDA=DNDA/(1.-BETAG*A/ATG)
  ENDIF
  IF (A .GT. ATG) DNDA=DNDA*EXP(((ATG-A)/ACG)**3)
  A01=3.5E-8
  A02=3.E-7
  SIG=0.4
  B1=2.0496E-7
  B2=9.6005E-11
  DNDAVSG=(B1/A)*EXP(-0.5*(LOG(A/A01)/SIG)**2)+ (B2/A)*EXP(-0.5*(LOG(A/A02)/SIG)**2)
  IF (DNDAVSG .GE. 0.0001*DNDA) DNDA=DNDA+BC5*DNDAVSG
ENDIF
END SUBROUTINE GRAIN_DIST


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Christophe Cossou
!
!> @date 2014
!
! DESCRIPTION: 
!> @brief subroutine that try to split the line in two part, given a 
!! separator value (set in parameter of the subroutine)
!
!> @warning The first character of the parameter value MUST NOT be a string. All spaces will be truncated
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine get_parameter_value(line, isParameter, id, value)

  implicit none
  
  ! Input
  character(len=80), intent(in) :: line !< [in] Input line in which we want to retrieve a parameter and its value
  
  ! Output
  logical, intent(out) :: isParameter !< [out] a boolean to say whether or not there is a parameter on this line. 
!! i.e if there is an occurence of the separator in the input line
  character(len=80), intent(out) :: id !< [out] the name of the parameter
  character(len=80), intent(out) :: value !< [out] a string that contains the value(s) associated with the parameter name. 
!!         Note that a special attention is given to the fact that the first character of 'value' must NOT be a 'space'
  
  ! Local
  character(len=1), parameter :: SEP = '=' ! the separator of a parameter line
  character(len=1) :: first_character
  integer :: id_first_char
  integer :: sep_position ! an integer to get the position of the separator

  !------------------------------------------------------------------------------

  sep_position = index(line, SEP)
  
  if (sep_position.ne.0) then
    isParameter = .true.
    id = line(1:sep_position-1)
    
    id_first_char = sep_position +1
    first_character = line(id_first_char:id_first_char)
    do while (first_character.eq.' ')
      id_first_char = id_first_char +1
      first_character = line(id_first_char:id_first_char)
    end do
    value = line(id_first_char:)
  else
    isParameter = .false.
  end if

end subroutine get_parameter_value


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author
!> Sacha Gavino
!
!> @date June 2019
!
! DESCRIPTION: read the file branching_ratios.in line by line and store the
!              value of the branching ratio in a REAL variable. This work is
!              necessary to be able to calculate the photodissociation rates
!              of species with several channels.
!
!              WARNING: make sure that the file branching_ratios.in is here !
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!subroutine read_br()
!  use utilities
!  use global_variables
!
!  implicit none
!
!  !----------------------- VARIABLES -----------------------
!  integer :: ios, nb_line, i,j
!  integer :: read_branching=9
!
!  ! Variables for the unordered reaction file
!  !character(len=11), dimension(MAX_COMPOUNDS,37) :: SYMBOLUO1
!  real(double_precision) :: br
!
!  character(len=80) :: branching
!  character(len=177) :: line_reactions
!  character(len=11) :: A1,A2,A3,A4,A5,A6,A7,A8
!  branching = "branching_ratios.in"
!  !---------------------------------------------------------
!
!
!  !------OPEN branching_ratios.in here------
!  open(unit=read_branching, file=branching, status="old", action="read", iostat=ios)
!  if (ios /= 0) stop "Error opening file branching_ratios.in. Check if branching_ratios.in is there."
!
!  !------GET the number of lines in branching_ratios.in------
!  !nb_line = 0
!  !DO
!  !READ (read_branching,*, END=10)
!  !nb_line = nb_line + 1
!  !END DO
!  !10 CLOSE (1)
!  !print *, "number of line in braching_ratios.in is:", nb_line
!
!  ! TO DO / READ THE NUMBER OF LINES IN ANOTHER PLACE AND ALLOCATE THE TABLE CORRECTLY. FOR THE
!  ! MOMENT, THE NUMBER IS GIVEn IN THE DEFINITION.
!
!  !------LOOP over the photodissociation reactions in branching_ratios.in------
!  !REWIND(read_branching)
!
!  do j = 1,nb_line_photorates
!      read(read_branching,100) A1,A2,A3,A4,A5,A6,A7,A8,br
!      reaction_compounds_names_photo(J,1) = A1
!      reaction_compounds_names_photo(J,2) = A2
!      reaction_compounds_names_photo(J,3) = A3
!      reaction_compounds_names_photo(J,4) = A4
!      reaction_compounds_names_photo(J,5) = A5
!      reaction_compounds_names_photo(J,6) = A6
!      reaction_compounds_names_photo(J,7) = A7
!      reaction_compounds_names_photo(J,8) = A8
!      br_photo(J) = br
!  enddo
!
!  100 FORMAT(3A11,1X,5A11,2X,F14.12)
!
!  !-----CLOSE open files-----
!  close(read_branching)
!end subroutine read_br
!
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!> @author
!!> Sacha Gavino
!!
!!> @date June 2019
!!
!! DESCRIPTION: Subroutine to read the species cross sections from folder
!!
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!subroutine read_spec_photo()
!  use global_variables
!  implicit none
!
!  integer :: i,j
!  character(len=11) :: A1
!  character(len=32) :: filename
!  real(kind=8), dimension (4) :: line
!
!  open(unit=9, file='cross-sections/param.txt', status="old", action="read")
!
!  do i=1,nb_line_spec_photorates
!    read(9,'(A)') A1
!    spec_photo(i) = A1
!    if (A1.eq.YH2) id_photo_H2 = i
!    if (A1.eq.YCO) id_photo_CO = i
!    if (A1.eq.YN2) id_photo_N2 = i
!
!    if (A1.eq.YH2O) id_photo_H2O = i
!    if (A1.eq.YCO2) id_photo_CO2 = i
!    if (A1.eq.YN2O) id_photo_N2O = i
!
!    if (A1.eq.YCH) id_photo_CH = i
!    if (A1.eq.YCH3) id_photo_CH3 = i
!    if (A1.eq.YCH4) id_photo_CH4 = i
!
!    if (A1.eq.YOH) id_photo_OH = i
!    if (A1.eq.YHCO) id_photo_HCO = i
!    if (A1.eq.YH2CO) id_photo_H2CO = i
!
!    if (A1.eq.YCN) id_photo_CN = i
!    if (A1.eq.YHCN) id_photo_HCN = i
!    if (A1.eq.YHNC) id_photo_HNC = i
!    
!    if (A1.eq.YNH) id_photo_NH = i
!    if (A1.eq.YNH2) id_photo_NH2 = i
!    if (A1.eq.YNH3) id_photo_NH3 = i
!
!
!  enddo
!
!  close(9)
!
!  do i=1,nb_line_spec_photorates
!      filename = 'cross-sections/'//trim(spec_photo(i))//'.txt'
!      open(10, file=filename)
!      do j = 1,6
!          read(10,*)
!      enddo
!      do j = 1,nb_line_table_flux
!         ! read(10,'(e12.6,3(1x,e12.6))')  line
!          read(10,*)  line
!          cross_sections(j,i,:) = line(:)
!      enddo
!      close(10)
!  enddo
!
!end subroutine read_spec_photo
!
!
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!> @author
!!> Sacha Gavino
!!
!!> @date June 2019
!!
!! DESCRIPTION: subroutine to read the ISRF, local UV, star flux
!!
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!subroutine read_flux_tables()
!  use global_variables
!  implicit none
!
!  integer :: i
!  character(len=11) :: A1
!  character(len=35) :: AA
!  real(double_precision), dimension(2) :: line
!  real(double_precision), dimension(spatial_resolution) :: line_uv_disk ! different dimension because local_flux_dust.txt has a column for each spatial point.
!  real(double_precision), dimension(spatial_resolution) :: line_uv_local ! different dimension because local_flux_dust.txt has a column for each point.
!
!  open(unit=11, file='flux/local_flux_dust.txt', status="old", action="read") ! this table is from the radiative transfer code POLARIS.
!  open(unit=12, file='flux/flux_isrf_dust.txt', status="old", action="read") ! this table is from the radiative transfer code POLARIS.
!
!  do i=1,nb_line_table_flux
!    read(11,*) line_uv_local
!    local_flux_dust(i,:) = line_uv_local(:)
!    read(12,*) line_uv_disk
!    flux_isrf_dust(i,:) = line_uv_disk(:)
!  enddo
!
!  close(11)
!  close(12)
!
!end subroutine read_flux_tables

end module input_output
