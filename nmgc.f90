! -----------------------------------------------------------------------
!   #     #     #     #     #  #######  ###  #        #     #   #####   
!   ##    #    # #    #     #     #      #   #        #     #  #     #  
!   # #   #   #   #   #     #     #      #   #        #     #  #        
!   #  #  #  #     #  #     #     #      #   #        #     #   #####   
!   #   # #  #######  #     #     #      #   #        #     #        #  
!   #    ##  #     #  #     #     #      #   #        #     #  #     #  
!   #     #  #     #   #####      #     ###  #######   #####    #####   
! -----------------------------------------------------------------------
!                            ?
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~|^"~~~~~~~~~~~~~~~~~~~~~~~~~o~~~~~~~~~~~
!        o                   |                  o      __o
!         o                  |                 o     |X__>
!       ___o                 |                __o
!     (X___>--             __|__            |X__>     o
!                         |     \                   __o
!                         |      \                |X__>
!  _______________________|_______\________________
! <                                                \____________   _
!  \                                                            \ (_)
!   \    O       O       O                                       >=)
!    \__________________________________________________________/ (_)
!
!                            ___
!                           / o \
!                      __   \   /   _
!                        \__/ | \__/ \
!                       \___//|\\___/\
!                        ___/ | \___
!                             |     \
!                            /
! -----------------------------------------------------------------------
! 
! A fast 1D gas-grain chemical model by FH (2008)
! Based upon the OSU gas-grain chemical model
! Uses the OSU chemical network
! Updated from gg_osu_2006v1d by RTG/VW
! Rate equations from Hasegawa & Herbst (1992)
! Modified rates following Stantcheva et al. (2001)
! Stiff solver for sparse Jacobians: LSODES/ODEPACK (Hindmarsh 1983)
! Turbulent mixing implemented through first order operator splitting
! 
! April 2011 VW  An appoximative calculation of the X-ray ionization rate 
! have been added according to Glasgold at al. (1999) -> DOES NOT WORK YET
!
! Strong modifications of the code by C. Cossou from Jan. to July 2014
!
! INPUT FILES
! 
! All parameters file can have comments, either a full line or the end of a line. The comment character being the '!' character. 
! Blanck lines are ignored
! 
! parameters.in : parameter file of the code, with various flags
! 
! abundances.in : Give initial abundances for a set of species (a reduced number or all. Default minimum values are applied for those
! that do not exist in this file.
! 
! activation_energies.in : Activation energy for endothermic reactions
! 
! element.in : name and mass in AMU of all elemental species existing in the simulation
! 
! gas_reactions.in : Reaction that occurs in gas phase
! 
! gas_species.in : species that are involved in gas phase reactions
! 
! grain_reactions.in : Reactions that occurs on grain surface
! 
! grain_species.in : Species that are involved in grain surface reactions
! 
! surface_parameters.in : various energies and parameters for diffusion and movements on the grain surface
! structure_evolution.dat : providing the structure of the source evolving with time in the case the model is used this way.
! 1D_static.dat : providing the 1D static physical structure of a source when the model is used this way.
!
! OUTPUT FILES
! *.out files are output files. *.tmp files are file that are generated at each timestep, either to continue a 
! simulation or check if there is a problem
!
! abundances.*.out : writing in binary format the abundances of all species at each timestep (one file per timestep-output)
!
! abundances.tmp : writing in ASCII format the abundances of all species at the current timestep-output
!
! rates.*.out : writing in binary format the critical reactions rates
!
! info.out : writing information on the code, in particular about the code version (commit ID and branch)
!
! species.out : Writing the list of species and their corresponding index
!
! elemental_abundances.tmp/out : writing information about prime elements, their total abundances and mass
!
!
! 
! -----------------------------------------------------------------------

PROGRAM Gasgrain

use global_variables
use iso_fortran_env
use shielding
use utilities
use dust_temperature_module
use nautilus_main
use photorates

implicit none

! Parameters for DLSODES. RTOL is the RELATIVE_ABUNDANCE parameter in global_variables.f90
integer :: itol = 2 !< ITOL = 1 or 2 according as ATOL (below) is a scalar or array.
integer :: itask = 1 !< ITASK = 1 for normal computation of output values of Y at t = TOUT.
integer :: istate = 1 !< ISTATE = integer flag (input and output).  Set ISTATE = 1.
integer :: iopt = 1 !< IOPT = 0 to indicate no optional inputs used.
integer :: mf = 121 !< method flag.  Standard values are:
!!\n          10  for nonstiff (Adams) method, no Jacobian used
!!\n          121 for stiff (BDF) method, user-supplied sparse Jacobian
!!\n          222 for stiff method, internally generated sparse Jacobian
real(double_precision) :: atol = 1.d-99 !< absolute tolerance parameter (scalar or array).
!!\n          The estimated local error in Y(i) will be controlled so as
!!\n          to be roughly less (in magnitude) than
!!\n             EWT(i) = RTOL*ABS(Y(i)) + ATOL     if ITOL = 1, or
!!\n             EWT(i) = RTOL*ABS(Y(i)) + ATOL(i)  if ITOL = 2.
!!\n          Thus the local error test passes if, in each component,
!!\n          either the absolute error is less than ATOL (or ATOL(i)),
!!\n          or the relative error is less than RTOL.
!!\n          Use RTOL = 0.0 for pure absolute error control, and
!!\n          use ATOL = 0.0 (or ATOL(i) = 0.0) for pure relative error
!!\n          control.  Caution: actual (global) errors may exceed these
!!\n          local tolerances, so choose them conservatively.
real(double_precision) :: output_timestep !< Timestep to reach the next output time, starting from old output time [s]
real(double_precision) :: integration_timestep !< Timestep [s] to have an accurate integration. 
!! This is smaller than output_timestep in case of diffusion, or equal to output_timestep

real(double_precision), dimension(:), allocatable :: output_times !< [s] dim(NB_OUTPUTS or nb_times in data file) 
real(double_precision) :: output_step !< [s] used to compute the list of output time, depending if its linear or log spaced
integer :: output_idx !< integer for the output time loop

real(double_precision) :: code_start_time, code_current_time, code_elapsed_time !< cpu time in seconds, allowing to predict the expected ending time of the simulation
real(double_precision) :: remaining_time !< estimated remaining time [s]

integer :: i,j,k, ic_i !< For loops
character(2) :: c_i !character variable to convert integer value of j to character j


call initialisation()

call initialize_work_arrays()

if (is_structure_evolution.eq.0) then
select case(OUTPUT_TYPE)
  case('linear')! nb_output is used
    allocate(output_times(NB_OUTPUTS))
    output_step = (STOP_TIME - START_TIME) / dfloat(NB_OUTPUTS - 1)
    do i=1, NB_OUTPUTS - 1
      output_times(i) = START_TIME + output_step * dfloat(i)
    enddo
    output_times(NB_OUTPUTS) = STOP_TIME ! To ensure the exact same final value
    

  case('log')! nb_output is used
    allocate(output_times(NB_OUTPUTS))
    
    output_step = (STOP_TIME/START_TIME) ** (1.d0/dfloat(NB_OUTPUTS - 1))
    do i=1, NB_OUTPUTS -1
      output_times(i) = START_TIME * output_step ** (i - 1.d0)
    enddo
    output_times(NB_OUTPUTS) = STOP_TIME ! To ensure the exact same final value

  case default
    write(error_unit,*) 'The OUTPUT_TYPE="', OUTPUT_TYPE,'" cannot be found.'
    write(error_unit,*) 'Values possible : linear, log, table'
    write(error_unit, '(a)') 'Error in nautilus: main program gasgrain' 
    call exit(12)
end select
endif

    ! We do not want 0 as first output time.
    if (is_structure_evolution.eq.1) then
        if (structure_time(1).eq.0.d0) then
            START_TIME = structure_time(2)
            STOP_TIME = structure_time(structure_sample)
            NB_OUTPUTS = structure_sample - 1
            allocate(output_times(NB_OUTPUTS))
            output_times(1:NB_OUTPUTS) = structure_time(2:structure_sample)
        else
            START_TIME = structure_time(1)
            STOP_TIME = structure_time(structure_sample)
            NB_OUTPUTS = structure_sample
            allocate(output_times(NB_OUTPUTS))
            output_times(1:NB_OUTPUTS) = structure_time(1:NB_OUTPUTS)
        endif
    endif

current_time = 0.d0 ! In seconds
call cpu_time(code_start_time)

! loop on output times (sub-step can exists, especially for diffusion process occuring on shorted timescales)
do output_idx=1, NB_OUTPUTS

  output_timestep = output_times(output_idx) - current_time
  
  call get_timestep(current_time=current_time, final_time=output_times(output_idx), next_timestep=integration_timestep)
  
  ! Loop on sub-timesteps. If diffusion (1D) then sub-timesteps are constrained by diffusion timescale.
  !! Else, only one step is done equal to the output_timestep
  do while (current_time.lt.output_times(output_idx))

    ! So far, structure read from ASCII file is not compatible with 1D. 
    !! If this changes, this call must be included in the loop on 1D sample
    !! We thus assume that everything is stored in the first element of the array, since this can't be valid for nb > 1
    if (spatial_resolution.eq.1) then
      call get_structure_properties(time=current_time, & ! Inputs
                              av=visual_extinction(1), density=H_number_density(1), & ! Outputs
                              gas_temperature=gas_temperature(1)) ! Outputs
    endif

    ! Column densities are set to zero before each time evolution computation.

    NH_z(1:spatial_resolution) = 0.d0
    NH2_z(1:spatial_resolution) = 0.d0
    NN2_z(1:spatial_resolution) = 0.d0
    NCO_z(1:spatial_resolution) = 0.d0
    NH2O_z(1:spatial_resolution) = 0.d0
    NCO2_z(1:spatial_resolution) = 0.d0
    NN2O_z(1:spatial_resolution) = 0.d0
    NCH_z(1:spatial_resolution) = 0.d0
    NCH3_z(1:spatial_resolution) = 0.d0
    NCH4_z(1:spatial_resolution) = 0.d0
    NOH_z(1:spatial_resolution) = 0.d0
    NHCO_z(1:spatial_resolution) = 0.d0
    NH2CO_z(1:spatial_resolution) = 0.d0
    NCN_z(1:spatial_resolution) = 0.d0
    NHCN_z(1:spatial_resolution) = 0.d0
    NHNC_z(1:spatial_resolution) = 0.d0
    NNH_z(1:spatial_resolution) = 0.d0
    NNH2_z(1:spatial_resolution) = 0.d0
    NNH3_z(1:spatial_resolution) = 0.d0

    do x_i=1,spatial_resolution

        ! For the case of 1D, some grain parameters are changed according to the spatial point:
        !IF (is_dust_1D.eq.1) THEN
        IF (STRUCTURE_TYPE.eq.'1D_no_diff') THEN  
            GTODN(:)=GTODN_1D(:,x_i)
            UV_flux = UV_flux_1D(x_i)
            grain_radii(:)=grain_radius_1D(:,x_i)
            if(nb_grains /=1)   CR_PEAK_GRAIN_TEMP_all(:)=CR_PEAK_GRAIN_TEMP_all_1D(:,x_i)
!             grain_radius=grain_radius_1D(x_i)
              do i=1,nb_grains
  !               SPECIES_MASS(INDGRAIN(i))=4.0*PI*grain_radius*grain_radius*grain_radius*GRAIN_DENSITY/3.0/AMU
                SPECIES_MASS(INDGRAIN(i))=4.0*PI*grain_radii(i)*grain_radii(i)*grain_radii(i)*GRAIN_DENSITY/3.0/AMU
  !               SPECIES_MASS(INDGRAIN_MINUS(i))=4.0*PI*grain_radius*grain_radius*grain_radius*GRAIN_DENSITY/3.0/AMU
                SPECIES_MASS(INDGRAIN_MINUS(i))=4.0*PI*grain_radii(i)*grain_radii(i)*grain_radii(i)*GRAIN_DENSITY/3.0/AMU
                nb_sites_per_grain(i) = SURFACE_SITE_DENSITY*4.d0*PI*grain_radii(i)**2
              enddo  
! VW Changed on september 2018. the NH is directly computed from the physical structure 
! so we do not need this. The most outer part is anyway set to Av = 0. If not, then 
! the AV_NH_ratio will be the typical ISM one.  
!            AV_NH_ratio=AV_NH_1D(x_i)
            endif

          call get_grain_temperature(space=x_i,time=current_time, gas_temperature=gas_temperature(x_i), &
                av=visual_extinction(x_i), & ! inputs
                grain_temperature=dust_temperature(:,x_i)) ! Outputs
          
          ! We can't add input variables in dlsodes called routines, so we must store the values as global variables
          actual_gas_temp = gas_temperature(x_i)
          actual_dust_temp(:) = dust_temperature(:,x_i)
          actual_av = visual_extinction(x_i)
          actual_gas_density = H_number_density(x_i)
      
!--------------------------------------------------------------
! below few lines of code scales computed grain temp for different grain radius 
! computed grain temperature is only for .1 micron grain so a scaling is done here
      if (trim(GRAIN_TEMPERATURE_TYPE)=='computed') then
!   write(*,*)grain_radius,  actual_dust_temp(1)
        write(*,*)
        do i=1,nb_grains
          actual_dust_temp(i) = actual_dust_temp(i)*(grain_radii(i)/grain_radius)**(-1.d0/6.d0)
!     write(*,*)grain_radii(i),actual_dust_temp(i)
        enddo
      endif
! pause
!--------------------------------------------------------------


      
      DO i=1,nb_grains
        IF ((actual_dust_temp(i).gt.15.0d+00).and.(is_er_cir.eq.1)) THEN
          PRINT*, "Warning, for grain",i,"with radius =",grain_radii(i)
          PRINT*, "the temperature range is above the validity domain of the complex induced reaction mechanism..."
        ENDIF
      ENDDO
      if ((x_i.eq.1)) then

          NH_z(x_i)= actual_av/AV_NH_ratio * abundances(indH, x_i)
          NH2_z(x_i)= actual_av/AV_NH_ratio * abundances(indH2, x_i)
          NN2_z(x_i)= actual_av/AV_NH_ratio * abundances(indN2, x_i)
          NCO_z(x_i)= actual_av/AV_NH_ratio * abundances(indCO, x_i)
          NH2O_z(x_i)= actual_av/AV_NH_ratio * abundances(indH2O, x_i)
          NCH_z(x_i)= actual_av/AV_NH_ratio * abundances(indCH, x_i)
          NCH3_z(x_i)= actual_av/AV_NH_ratio * abundances(indCH3, x_i)
          NH2CO_z(x_i)= actual_av/AV_NH_ratio * abundances(indH2CO, x_i)
          NCO2_z(x_i)= actual_av/AV_NH_ratio * abundances(indCO2, x_i)
          NN2O_z(x_i)= actual_av/AV_NH_ratio * abundances(indN2O, x_i)
          NCH4_z(x_i)= actual_av/AV_NH_ratio * abundances(indCH4, x_i)
          NOH_z(x_i)= actual_av/AV_NH_ratio * abundances(indOH, x_i)
          NHCO_z(x_i)= actual_av/AV_NH_ratio * abundances(indHCO, x_i)
          NCN_z(x_i)= actual_av/AV_NH_ratio * abundances(indCN, x_i)
          NHCN_z(x_i)= actual_av/AV_NH_ratio * abundances(indHCN, x_i)
          NHNC_z(x_i)= actual_av/AV_NH_ratio * abundances(indHNC, x_i)
          NNH_z(x_i)= actual_av/AV_NH_ratio * abundances(indNH, x_i)
          NNH2_z(x_i)= actual_av/AV_NH_ratio * abundances(indNH2, x_i)
          NNH3_z(x_i)= actual_av/AV_NH_ratio * abundances(indNH3, x_i)

      else

          grid_cell_size = abs(grid_sample(x_i) - grid_sample(x_i-1))

          NH_z(x_i) = NH_z(x_i-1) + actual_gas_density * grid_cell_size * abundances(indH, x_i)
          NH2_z(x_i) = NH2_z(x_i-1) + actual_gas_density * grid_cell_size * abundances(indH2, x_i)
          NN2_z(x_i) = NN2_z(x_i-1) + actual_gas_density * grid_cell_size * abundances(indN2, x_i)
          NCO_z(x_i) = NCO_z(x_i-1) + actual_gas_density * grid_cell_size * abundances(indCO, x_i)
          NH2O_z(x_i) = NH2O_z(x_i-1) + actual_gas_density * grid_cell_size * abundances(indH2O, x_i)
          NCH_z(x_i) = NCH_z(x_i-1) + actual_gas_density * grid_cell_size * abundances(indCH, x_i)
          NCH3_z(x_i) = NCH3_z(x_i-1) + actual_gas_density * grid_cell_size * abundances(indCH3, x_i)
          NH2CO_z(x_i) = NH2CO_z(x_i-1) + actual_gas_density * grid_cell_size * abundances(indH2CO, x_i)
          NCO2_z(x_i) = NCO2_z(x_i-1) + actual_gas_density * grid_cell_size * abundances(indCO2, x_i)
          NN2O_z(x_i) = NN2O_z(x_i-1) + actual_gas_density * grid_cell_size * abundances(indN2O, x_i)
          NCH4_z(x_i) = NCH4_z(x_i-1) + actual_gas_density * grid_cell_size * abundances(indCH4, x_i)
          NOH_z(x_i) = NOH_z(x_i-1) + actual_gas_density * grid_cell_size * abundances(indOH, x_i)
          NHCO_z(x_i) = NHCO_z(x_i-1) + actual_gas_density * grid_cell_size * abundances(indHCO, x_i)
          NCN_z(x_i) = NCN_z(x_i-1) + actual_gas_density * grid_cell_size * abundances(indCN, x_i)
          NHCN_z(x_i) = NHCN_z(x_i-1) + actual_gas_density * grid_cell_size * abundances(indHCN, x_i)
          NHNC_z(x_i) = NHNC_z(x_i-1) + actual_gas_density * grid_cell_size * abundances(indHNC, x_i)
          NNH_z(x_i) = NNH_z(x_i-1) + actual_gas_density * grid_cell_size * abundances(indNH, x_i)
          NNH2_z(x_i) = NNH2_z(x_i-1) + actual_gas_density * grid_cell_size * abundances(indNH2, x_i)
          NNH3_z(x_i) = NNH3_z(x_i-1) + actual_gas_density * grid_cell_size * abundances(indNH3, x_i)

      endif

      NH = NH_z(x_i)
      NH2 = NH2_z(x_i)
      NN2 = NN2_z(x_i)
      NCO = NCO_z(x_i)
      NH2O = NH2O_z(x_i)
      NCH = NCH_z(x_i)
      NCH3 = NCH3_z(x_i)
      NH2CO = NH2CO_z(x_i)
      NCO2 = NCO2_z(x_i)
      NN2O = NN2O_z(x_i)
      NCH4 = NCH4_z(x_i)
      NOH = NOH_z(x_i)
      NHCO = NHCO_z(x_i)
      NCN = NCN_z(x_i)
      NHCN = NHCN_z(x_i)
      NHNC = NHNC_z(x_i)
      NNH = NNH_z(x_i)
      NNH2 = NNH2_z(x_i)
      NNH3 = NNH3_z(x_i)


      call integrate_chemical_scheme(delta_t=output_timestep, temp_abundances=abundances(1:nb_species, x_i),& ! Inputs
      itol=itol, atol=atol, itask=itask, iopt=iopt, mf=mf, & ! Inputs
      istate=istate) ! Output

      reaction_rates_1D(x_i,1:nb_reactions)=reaction_rates(1:nb_reactions)

    if (istate.eq.-3) stop

    ! Prevent too low abundances
    sumlaysurfsave = 0.0d0
    sumlaymantsave = 0.0d0
    do i=1,nb_species
      if (abundances(i,x_i).le.1.d-99) then
        abundances(i,x_i) = 1.d-99
      endif
      if(species_name(i)(1:1).eq."J") then 
        c_i = species_name(i)(2:3)  
        read(c_i,'(I2)')ic_i
        sumlaysurfsave(ic_i) = sumlaysurfsave(ic_i) + abundances(i,x_i)
!         do j=1,nb_grains
!           write(i_c,'(I2.2)')j
!           if(species_name(i)(2:3) == i_c) sumlaysurfsave(j) = sumlaysurfsave(j) + abundances(i,x_i)
!         enddo
      endif  
      if(species_name(i)(1:1).eq."K") then
        c_i = species_name(i)(2:3)  
        read(c_i,'(I2)')ic_i
        sumlaymantsave(ic_i) = sumlaymantsave(ic_i) + abundances(i,x_i)
!         do j=1,nb_grains
!           write(i_c,'(I2.2)')j
!           if(species_name(i)(2:3) == i_c)       sumlaymantsave(j) = sumlaymantsave(j) + abundances(i,x_i)
!         enddo
      endif  
    enddo

    
    do i=1,nb_grains
       sumlaysurfsave(i) = sumlaysurfsave(i) * GTODN(i) / nb_sites_per_grain(i)
       sumlaymantsave(i) = sumlaymantsave(i) * GTODN(i) / nb_sites_per_grain(i)
    enddo  
    !    where(abundances <= 1.d-99) abundances=1.d-99
    WRITE(*,'(a)')'----------------------------------------------------------------------------------------------'
    if ((structure_type.eq.'1D_no_diff').or.(structure_type.eq.'1D_diff')) then
      WRITE(*,'(a,4x,i4,4x,a)') " ---- SPATIAL POINT = ",x_i,"----"
    endif
    WRITE(*,'(2x,a)') "Physical parameters :"
    WRITE(*,'(4x,a19,ES13.6,1x,a6,4x,a19,ES13.6,1x,a6)') "Density = ",actual_gas_density,"[cm-3]","Av = ",actual_av,"[mag]"
    WRITE(*,'(4x,a16,ES13.6,1x,a6,4x,a27,ES13.6,1x,a10)') "Tgas = ",actual_gas_temp,"[K]", "UV flux = ", UV_flux,"[standard]"
    IF(is_3_phase.eq.0) THEN
      WRITE(*,'(2x,a)') "2 phase model (Rq. Nb. mant. layer should be equal to 0.00):"
    ELSE
      WRITE(*,'(2x,a59,f5.2,a)') "3 phase model (Nb. surf. layer should never be >",nb_active_lay,"):"
    ENDIF

    !WRITE(*,'(a8, ES13.6)') "GTODN = ", GTODN(1)
    !WRITE(*,'(a8, ES13.6)') "GTODN = ", GTODN(2)
    WRITE(*,'(2x, a)')'Grain parameters :'
    do i=1,nb_grains
!       WRITE(*,'(4x,a18,x,i3)')'Results for grain',i
      WRITE(*,'(a25,i2, a2)')'Grain rank =',i,': '
      WRITE(*,91) "GTODN = ", GTODN(i),"Grain radius =",grain_radii(i),"[cm]","Grain temperature =",actual_dust_temp(i),"[K]"
! WRITE(*,*)"Nb. surf. layer = ",sumlaysurfsave(i),"Nb. tot.  layer = ",sumlaysurfsave(i)+sumlaymantsave(i)
      WRITE(*,92) "Nb. surf. layer = ",sumlaysurfsave(i),"Nb. mant. layer = ", sumlaymantsave(i),&
                  &"Nb. total layer = ",sumlaysurfsave(i)+sumlaymantsave(i)
! WRITE(*,*) "Nb. mant. layer = ", sumlaymantsave(i) 
      WRITE(*,'(a)') " "
    enddo
    enddo
    WRITE(*,'(a)')'=============================================================================================='
!91  FORMAT(2x,a12,x,i3,4x,a14,ES10.2,x,a4,4x,a19,F8.2,x,a4) 
91  FORMAT(2x,a19,x,ES10.3,4x,a14,ES10.2,x,a4,4x,a19,F8.2,x,a4) 
92  FORMAT(2x,a29,x,F10.6,2x,a18,x,F10.6,2x,a18,x,F10.6)
    call write_H2_CO_col_dens(index=output_idx)


    if (spatial_resolution.eq.1) then
      call check_conservation(abundances(1:nb_species, 1))
    endif
    
    ! We did one timestep. We update the diffusion for this timestep (mainly species diffusion)
    call structure_diffusion(timestep=integration_timestep, temp_abundances=abundances(1:nb_species, 1:spatial_resolution))

    current_time = current_time + integration_timestep ! New current time at which abundances are valid
    
    ! We compute the next integration timestep
    call get_timestep(current_time=current_time, final_time=output_times(output_idx), next_timestep=integration_timestep)

  enddo
  
  call cpu_time(code_current_time)
  code_elapsed_time =  code_current_time - code_start_time
  remaining_time = code_elapsed_time*NB_OUTPUTS/output_idx-code_elapsed_time
  
  if (remaining_time.lt.60.d0) then
    write(Output_Unit,'(a,en11.2e2,a,f5.1,a,f5.1,a)') 'T =',current_time/YEAR,&
    ' years [', 100.d0 * output_idx/NB_OUTPUTS, ' %] ; Estimated time remaining: ', remaining_time, ' s'
  else if (remaining_time.lt.3600.d0) then
    write(Output_Unit,'(a,en11.2e2,a,f5.1,a,f5.1,a)') 'T =',current_time/YEAR,&
    ' years [', 100.d0 * output_idx/NB_OUTPUTS, ' %] ; Estimated time remaining: ', remaining_time / 60.d0, ' min.'
  else
    write(Output_Unit,'(a,en11.2e2,a,f5.1,a,f5.1,a)') 'T =',current_time/YEAR,&
    ' years [', 100.d0 * output_idx/NB_OUTPUTS, ' %] ; Estimated time remaining: ', remaining_time / 3600.d0, ' h'
  endif

    ! Output of the rates
    call write_current_rates(index=output_idx)

    ! Output of the abundances
    call write_current_output(index=output_idx)




  first_step_done = .true.
enddo

if (spatial_resolution.eq.1) call write_abundances('abundances.tmp')

contains



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Franck Hersant
!
!> @date 2003
!
! DESCRIPTION: 
!> @brief Chemically evolve for a given time delta_t
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine integrate_chemical_scheme(delta_t,temp_abundances,itol,atol,itask,iopt,mf,istate)

  use global_variables
  
  implicit none

  ! Inputs
  real(double_precision), intent(in) :: delta_t !<[in] time during which we must integrate
  integer, intent(in) :: itol !<[in] ITOL   = 1 or 2 according as ATOL (below) is a scalar or array.
  integer, intent(in) :: itask !<[in] ITASK  = 1 for normal computation of output values of Y at t = TOUT.
  integer, intent(in) :: iopt !<[in] IOPT   = 0 to indicate no optional inputs used.
  integer, intent(in) :: mf !<[in] method flag.  Standard values are:
!!\n          10  for nonstiff (Adams) method, no Jacobian used
!!\n          121 for stiff (BDF) method, user-supplied sparse Jacobian
!!\n          222 for stiff method, internally generated sparse Jacobian
  real(double_precision), intent(in) :: atol !<[in] integrator tolerance

  ! Outputs
  integer, intent(out) :: istate !<[out] ISTATE = 2  if DLSODES was successful, negative otherwise.
!!\n          -1 means excess work done on this call (perhaps wrong MF).
!!\n          -2 means excess accuracy requested (tolerances too small).
!!\n          -3 means illegal input detected (see printed message).
!!\n          -4 means repeated error test failures (check all inputs).
!!\n          -5 means repeated convergence failures (perhaps bad Jacobian
!!\n             supplied or wrong choice of MF or tolerances).
!!\n          -6 means error weight became zero during problem. (Solution
!!\n             component i vanished, and ATOL or ATOL(i) = 0.)
!!\n          -7 means a fatal error return flag came from sparse solver
!!\n             CDRV by way of DPRJS or DSOLSS.  Should never happen.
!!\n          A return with ISTATE = -1, -4, or -5 may result from using
!!\n          an inappropriate sparsity structure, one that is quite
!!\n          different from the initial structure.  Consider calling
!!\n          DLSODES again with ISTATE = 3 to force the structure to be
!!\n          reevaluated. 
 
  ! Input/Output
  real(double_precision), dimension(nb_species), intent(inout) :: temp_abundances !<[in,out] Temporary array that contain the 
  !! abundances for the duration of the timestep, a sort of buffer.

  ! Locals
  integer :: i
  real(double_precision) :: t !< The local time, starting from 0 to delta_t
  real(double_precision), dimension(nb_species) :: satol !< Array that contain the absolute tolerance, 
  !! one value per species, either a minimum value or a value related to its abundance.


  real(double_precision) :: t_stop_step

  t_stop_step = delta_t
  t = 0.d0

  do while (t.lt.t_stop_step)

    istate = 1

    ! Adaptive absolute tolerance to avoid too high precision on very abundant species,
    ! H2 for instance. Helps running a bit faster

    do i=1,nb_species
      satol(i) = max(atol, 1.d-16 * temp_abundances(i))
    enddo

    ! Added on Dec 2019 to prevent from NaN issues.
    if (temp_abundances(i).le.1.d-60) then
        temp_abundances(i) = 1.d-60
    end if

    ! Feed IWORK with IA and JA

    call set_work_arrays(Y=temp_abundances)
    

    call dlsodes(get_temporal_derivatives,nb_species,temp_abundances,t,t_stop_step,itol,RELATIVE_TOLERANCE,&
    satol,itask,istate,iopt,rwork,lrw,iwork,liw,get_jacobian,mf)       



    ! Whenever the solver fails converging, print the reason.
    ! cf odpkdmain.f for translation
    if (istate.ne.2) then
      write(*,*)  'ISTATE = ', ISTATE
    endif

  enddo

  return 
  end subroutine integrate_chemical_scheme

! ======================================================================
! ======================================================================
END program gasgrain
