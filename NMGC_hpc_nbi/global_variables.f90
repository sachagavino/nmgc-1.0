!******************************************************************************
! MODULE: global_variables
!******************************************************************************
!
! DESCRIPTION: 
!> @brief Contains global variables and parameters
!
!******************************************************************************

module global_variables
use numerical_types
use iso_fortran_env
use utilities

implicit none

!Integer parameter to define number of grains to be used in simulation.
! integer, parameter     :: nb_grains=1
integer                  :: nb_grains
integer                  :: multi_grain
real(double_precision),allocatable, dimension(:):: grain_radii,grain_temp
! real(double_precision), dimension(1:):: grain_radii
! Theses 3 parameters are only intnb_line_table_fluxended to easy the transition when one want to add a reactant or a
!! product (product are fairly easy, reactant are not). But some things needed to be modified in the 
!! code anyway when you want to add either a reactant or a product
integer, parameter :: MAX_REACTANTS = 3 !< The maximum number of reactants for one reaction. 
!! Do not think that changing this parameter alone is sufficient to allow the code to handle it directly !!
integer, parameter :: MAX_PRODUCTS = 5 !< The maximum number of products for one reaction. 
!! Do not think that changing this parameter alone is sufficient to allow the code to handle it directly !!
integer, parameter :: MAX_COMPOUNDS = MAX_REACTANTS + MAX_PRODUCTS !< Total maximum number of compounds for one reaction (reactants + products)
!! Warning: If this number change, get_jacobian(N, T, Y, J, IAN, JAN, PDJ) must be actualised, since each reactant and product has
!! its own variable, a new one must be created for the new column possible. 

integer :: nb_reactions !< total number of reactions
integer :: nb_species !< total number of species. Species are ordered. First, there are the nb_gaseous_species ones, and after, the nb_surface_species ones.
integer :: nb_gaseous_species !< number of species that are gaseous
integer :: nb_surface_species !< number of species that are on the surface of grains
integer :: NB_PRIME_ELEMENTS !< Number of prime element that compose molecules, such as H, He, C and so on.
integer :: nb_species_for_grain !< number of species involved in grain surface reactions (can be gas or grain phase elements)
integer :: nb_surface_reactions !< number of reactions on the grain surface
integer :: nb_species_for_gas !< number of species involved in gas phase reactions (can be gas or grain phase elements)
integer :: nb_gas_phase_reactions !< number of reactions in gas phase
real(double_precision), parameter :: MINIMUM_RATE_COEFFICIENT=1.0D-99 !< Minimum rate coefficient (Below, coefficients are forced to 0)
real(double_precision), parameter :: K_B = 1.3806488d-16 !< Boltzmann constant in CGS (cm^2 g s^⁻2 K-1)
real(double_precision), parameter :: PLANCK_CONSTANT = 6.62565d-27 !< Planck constant in CGS (erg.s or g cm2 s-1)
real(double_precision), parameter :: SPEED_OF_LIGHT = 2.99792458d10 !< speed of light in CGS (cm/s)
real(double_precision), parameter :: PI = 3.1415926535898d0 !< The number Pi
real(double_precision), parameter :: H_BARRE = 1.054571628d-27 !< Reduced Planck constant h/2*pi in CGS (g cm2 s-1)
real(double_precision), parameter :: AMU = 1.66053892d-24 !< Atomic mass unit in g
real(double_precision), parameter :: ELECTRON_MASS = 0.000548579909d0 !< Electron mass in AMU (close to 1/1836)
real(double_precision), parameter :: AVOGADRO = 6.02214129d23 !< avogadro number : number of atom in 1 mol
real(double_precision), parameter :: YEAR = 3.15576d7 !< one year in seconds
real(double_precision), parameter :: AU = 1.49597871d13 !< Astronomical unit in cm (mean earth-sun distance)

!-----Variables for protoplanetary disk models and photorates computation-----
integer(double_precision), parameter :: nb_line_photorates = 134 ! number of lines of the files for br branching ratios
integer(double_precision), parameter :: nb_line_spec_photorates = 80 ! number of species that have a cross-sections file.
integer(double_precision), parameter :: wv_max_factor = 309 ! number of first lines to integrate over for UV factor calculation because we don't take the full wv range.
integer(double_precision), parameter :: nb_line_table_flux = 309 ! number of lines in the fluxes files. Must find another way to get this value.
real(double_precision) :: INT_LOCAL_FLUX ! UV factor.
!----------------------------------------------------------------------------

real(double_precision) :: RELATIVE_TOLERANCE !< relative tolerance parameter (scalar) of DLSODES (of the ODEPACK package to solve ODE's)
real(double_precision), allocatable, dimension(:) :: EVAPORATION_RATES_H2 !< evaporation rate for H2 [s-1]. We sum all
! real(double_precision), dimension(1:nb_grains) :: EVAPORATION_RATES_H2 !< evaporation rate for H2 [s-1]. We sum all
!! the processes that are not time dependent. This is used for the encounter desorption process
real(double_precision), allocatable, dimension(:) :: EVAPORATION_RATES_TEMPO_H2 !< evaporation rate.Here we add to EVAPORATION_RATES_H2
! real(double_precision), dimension(1:) :: EVAPORATION_RATES_TEMPO_H2 !< evaporation rate.Here we add to EVAPORATION_RATES_H2
!! the time depedent processes.
real(double_precision) :: ED_H2 !< Binding energy of H2 over itself

! Name of key species
character(len=11) :: YH     = 'H          ' !< Gas phase Hydrogen
character(len=11) :: YN2    = 'N2         ' !< Gas phase N2
character(len=11) :: YJH    = 'J01H         ' !< Hydrogen on grains
character(len=11) :: YH2    = 'H2         ' !< Gas phase Dihydrogen
character(len=11) :: YJH2   = 'J01H2        ' !< Dihydrogen on grains
character(len=11) :: YJO    = 'J01O         ' !< Oxygen on grains
character(len=11) :: YHE    = 'He         ' !< Gas phase Helium
character(len=11) :: YHEP   = 'He+        ' !< Gas phase Helium+
character(len=11) :: YE     = 'e-         ' !< Gas phase electrons
! character(len=11) :: YGRAIN = 'GRAIN0     ' !< Grain
! character(len=11), dimension(1:nb_grains) :: YGRAIN      !< modified version of above line to handle many dust grains
character(len=11), allocatable, dimension(:) :: YGRAIN      !< modified version of above line to handle many dust grains
character(len=11), allocatable, dimension(:) :: YGRAIN_MINUS     !< modified version of above line to handle many dust grains
! character(len=11), dimension(1:nb_grains) :: YGRAIN_MINUS     !< modified version of above line to handle many dust grains

character(len=11) :: YCO    = 'CO         ' !< Gas phase CO
character(len=11) :: YKCO   = 'K01CO      ' !< CO in the mantle of the first grain
character(len=11) :: YJCO   = 'J01CO      ' !< CO on the surface of the first grain
character(len=11) :: YH2O   = 'H2O        ' !< Gas phase H2O
character(len=11) :: YKH2O  = 'K01H2O     ' !< H2O in the mantle of the first grain
character(len=11) :: YJH2O  = 'J01H2O     ' !< H2O on the surface of the first grain
character(len=11) :: YKH    = 'K01H       ' !< Hydrogen on first grain
character(len=11) :: YKH2   = 'K01H2      ' !< Dihydrogen on first grain
character(len=11) :: YCH   = 'CH         ' !< Gas phase CH
character(len=11) :: YCH3   = 'CH3        ' !< Gas phase CH3
character(len=11) :: YCO2   = 'CO2        ' !< Gas phase CO2
character(len=11) :: YH2CO   = 'H2CO       ' !< Gas phase H2CO
character(len=11) :: YN2O   = 'N2O       ' !< Gas phase N2O
character(len=11) :: YCH4   = 'CH4       ' !< Gas phase CH4
character(len=11) :: YOH   = 'OH       ' !< Gas phase OH
character(len=11) :: YHCO   = 'HCO       ' !< Gas phase HCO
character(len=11) :: YCN   = 'CN       ' !< Gas phase CN
character(len=11) :: YHCN   = 'HCN       ' !< Gas phase HCN
character(len=11) :: YHNC   = 'HNC       ' !< Gas phase HNC
character(len=11) :: YNH   = 'NH       ' !< Gas phase NH
character(len=11) :: YNH2   = 'NH2       ' !< Gas phase NH2
character(len=11) :: YNH3   = 'NH3       ' !< Gas phase NH3



integer :: INDCO   !< Index corresponding to CO in nb_species length arrays
integer :: INDKCO  !< Index corresponding to CO on the first grain in nb_species length arrays
integer :: INDJCO  !< Index corresponding to CO on the first grain in nb_species length arrays
integer :: INDN2   !< Index corresponding to N2 in nb_species length arrays
integer :: INDH2   !< Index corresponding to H2 in nb_species length arrays
integer :: INDH2O  !< Index corresponding to H2O in nb_species length arrays
integer :: INDKH2O !< Index corresponding to H2O on the first grain in nb_species length arrays
integer :: INDJH2O !< Index corresponding to H2O on the first grain in nb_species length arrays
integer :: INDCH  !< Index corresponding to CH in nb_species length arrays
integer :: INDCH3  !< Index corresponding to CH3 in nb_species length arrays
integer :: INDH2CO  !< Index corresponding to H2CO in nb_species length arrays
integer :: INDCO2  !< Index corresponding to CO2 in nb_species length arrays
integer :: INDN2O  !< Index corresponding to N2O in nb_species length arrays
integer :: INDCH4  !< Index corresponding to CH4 in nb_species length arrays
integer :: INDOH  !< Index corresponding to OH in nb_species length arrays
integer :: INDHCO  !< Index corresponding to HCO in nb_species length arrays
integer :: INDCN  !< Index corresponding to CN in nb_species length arrays
integer :: INDHCN  !< Index corresponding to HCN in nb_species length arrays
integer :: INDHNC  !< Index corresponding to HNC in nb_species length arrays
integer :: INDNH  !< Index corresponding to NH in nb_species length arrays
integer :: INDNH2  !< Index corresponding to NH2 in nb_species length arrays
integer :: INDNH3  !< Index corresponding to NH3 in nb_species length arrays
integer :: INDH    !< Index corresponding to H in nb_species length arrays
integer :: INDHE   !< Index corresponding to He in nb_species length arrays
integer :: INDEL   !< Index corresponding to e- in nb_species length arrays
integer, allocatable, dimension(:) :: INDGRAIN !< Index corresponding to GRAIN0 in nb_species length arrays
! integer, dimension(1:nb_grains) :: INDGRAIN !< Index corresponding to GRAIN0 in nb_species length arrays
integer, allocatable, dimension(:) :: INDGRAIN_MINUS !< Index corresponding to GRAIN- in nb_species length arrays
! integer, dimension(1:nb_grains) :: INDGRAIN_MINUS !< Index corresponding to GRAIN- in nb_species length arrays
integer, allocatable, dimension(:) :: GRAIN_RANK ! dim(nb_reactions) it says which size of grain is involved in a particular reaction

! 3 phase model
REAL(double_precision), allocatable, dimension(:) :: sumlaysurfsave !< Total number of layers on the grain surface
! REAL(double_precision), dimension(1:nb_grains) :: sumlaysurfsave !< Total number of layers on the grain surface
REAL(double_precision), allocatable, dimension(:) :: sumlaymantsave !< Total number of layers on the grain mantle
! REAL(double_precision), dimension(1:nb_grains) :: sumlaymantsave !< Total number of layers on the grain mantle
real(double_precision), allocatable, dimension(:) :: rate_tot_acc,rate_tot_des, rate_tot
real(double_precision) :: nb_active_lay

! Arrays about prime elements
real(double_precision), allocatable, dimension(:) :: INITIAL_ELEMENTAL_ABUNDANCE !< dim(NB_PRIME_ELEMENTS) Store abundances 
!! (relative to H) for all elemental species before running the simulation [number ratio]
real(double_precision), allocatable, dimension(:) :: elemental_mass !< dim(NB_PRIME_ELEMENTS) elemental mass [a.m.u]
character(len=11), allocatable, dimension(:) :: element_name !< dim(NB_PRIME_ELEMENTS) name of the prime elements
integer, allocatable, dimension(:) :: PRIME_ELEMENT_IDX ! < dim(NB_PRIME_ELEMENTS) Tell for each prime element its index in the global array of all elements.

! Arrays about species
character(len=11), allocatable, dimension(:) :: species_name !< dim(nb_species)
integer, allocatable, dimension(:,:) :: species_composition !< dim(NB_PRIME_ELEMENTS,nb_species) number of atom of each element composition the given species.
real(double_precision), allocatable, dimension(:,:) :: abundances !< dim(nb_species,spatial_resolution) Species abundances (relative to H) [number ratio]
real(double_precision), allocatable, dimension(:) :: SPECIES_MASS !< dim(nb_species) Species mass [a.m.u]
real(double_precision), allocatable, dimension(:) :: STICK_SPEC !< dim(nb_species) Species sticking probability
real(double_precision), allocatable, dimension(:) :: THERMAL_HOPING_RATE !< dim(nb_species) Diffusion rate by thermal hopping [s-1]:
!! 1/time required for an adsorbed particle to sweep over a number of sites equivalent to the whole grain surface 
real(double_precision), allocatable, dimension(:) :: CR_HOPING_RATE !< dim(nb_species) Diffusion rate by thermal hopping due to cosmic rays heating [s-1]:
!! 1/time required for an adsorbed particle to sweep over a number of sites equivalent to the whole grain surface
real(double_precision), allocatable, dimension(:) :: ACCRETION_RATES !< dim(nb_species) Accretion rate for a given species onto the grain surface [s-1]
real(double_precision), allocatable, dimension(:) :: EVAPORATION_RATES !< dim(nb_species) evaporation rate for a given species [s-1]. We sum all
!! the processes that are not time dependent.
real(double_precision), allocatable, dimension(:) :: EVAPORATION_RATES_TEMPO !< dim(nb_species) evaporation rate.Here we add to EVAPORATION_RATES
!! the time depedent processes.
real(double_precision), allocatable, dimension(:) :: BINDING_ENERGY !< dim(nb_species) [K] Binding energy of a species to the surface (specific to each species). Parameter read in the file surface_parameters.in
real(double_precision), allocatable, dimension(:) :: DIFFUSION_BARRIER !< dim(nb_species) [K] Potential energy barrier between adjacent surface potential energy wells. 
!! It is usually a fraction of the binding energy. This parameter is used to compute the rate for thermal hopping and the diffusion 
!! by quantum tunneling using Hasegawa’s formalism. Parameter read in the file surface_parameters.in for some species and computed 
!! in the model for the others. (specific of each species)
real(double_precision), allocatable, dimension(:) :: GAP_ENERGY_BANDS !< dim(nb_species) [K] gap between the energy bands corresponding 
!! to the fundamental and first excited state (Ricca et al., 1969). This parameter is used to compute the diffusion of species on 
!! the surface by tunneling using the formalism by Hollenbach & Salpeter (1970) and Watson (1976). It is the parameter dEb in the 
!! equation 20 of Reboussin et al. (2014). Values of this parameters from the literature are also given in the paper.
real(double_precision), allocatable, dimension(:) :: FORMATION_ENTHALPY !< dim(nb_species) Enthalpy of formation of the species 
!! (read in kcal/mol in the file surface_parameters.in and then converted into Kelvin/reaction via DHFSUM)
real(double_precision), allocatable, dimension(:) :: VIBRATION_FREQUENCY !< dim(nb_species) Characteristic vibration frequency [s-1] of the adsorbed species  as from a harmonic oscillator hypothesis (Hasegawa & Herbst 1992)
real(double_precision), allocatable, dimension(:) :: ACC_RATES_PREFACTOR !< dim(nb_reactions) Interrim calculation variable for ACCRETION_RATES
!! that contain cross section, probability and constant needed for thermal motion
real(double_precision), allocatable, dimension(:) :: TUNNELING_RATE_TYPE_1 !< dim(nb_species) Quantum tunneling diffusion rate [s-1] (Watson 1976) (dEB.BOLTZ) / (4.HBAR.nb_sites_per_grain)
real(double_precision), allocatable, dimension(:) :: TUNNELING_RATE_TYPE_2 !< dim(nb_species) Quantum tunneling diffusion rate [s-1] (Hasegawa & Herbst 1992) VIBRATION_FREQUENCY / nb_sites_per_grain.EXP(-2.DIFFUSION_BARRIER_THICKNESS / HBAR.(2.AMU.SMA.BOLTZ.EB)^1/2)
integer, allocatable, dimension(:) :: SPECIES_CHARGE !< dim(nb_species) !< electric charge [in e-] for each species, 0 if neutral, positive or negative if ions.

! Arrays about reactions
character(len=11), allocatable, dimension(:,:) :: REACTION_COMPOUNDS_NAMES !< dim(MAX_COMPOUNDS,nb_reactions). Empty string means no species
integer, allocatable, dimension(:,:) :: REACTION_COMPOUNDS_ID !< dim(MAX_COMPOUNDS, nb_reactions) for all reactions, 
!! list for reactants (first 3) and products (last 5). "nb_species+1" means no species
real(double_precision), allocatable, dimension(:) :: branching_ratio !< dim(nb_reactions) Branching ratio of each reaction
real(double_precision), allocatable, dimension(:) :: RATE_A !< dim(nb_reactions) Coefficient used to compute the reaction rate. Formula (and unit) is different in function of the reaction type.
real(double_precision), allocatable, dimension(:) :: RATE_B !< dim(nb_reactions) Coefficient used to compute the reaction rate. Formula (and unit) is different in function of the reaction type.
real(double_precision), allocatable, dimension(:) :: RATE_C !< dim(nb_reactions) Coefficient used to compute the reaction rate. Formula (and unit) is different in function of the reaction type.
real(double_precision), allocatable, dimension(:) :: reaction_rates !< dim(nb_reactions) reaction rate [unit depend on the reaction]
real(double_precision), allocatable, dimension(:,:) :: reaction_rates_1D  !< dim(spatial_resolution,nb_reactions) reaction rate for all spatial point at one time step [unit depend on the reaction]
real(double_precision), allocatable, dimension(:) :: ACTIVATION_ENERGY !< dim(nb_reactions) Activation energy for surface reactions [K]
real(double_precision), allocatable, dimension(:) :: REACTION_TMIN !< dim(nb_reactions) min temperature boundary of each reactions [K]
real(double_precision), allocatable, dimension(:) :: REACTION_TMAX !< dim(nb_reactions) max temperature boundary of each reactions [K]
real(double_precision), allocatable, dimension(:) :: SURF_REACT_PROBA !< dim(nb_reactions) probability for the surface reactions to 
!! happen upon an encounter. The probability is unity for exothermic reactions without activation energy and has to be computed 
!! otherwise. See equation 6 from Hasewaga, Herbst & Leung (1992).
integer, allocatable, dimension(:) :: REACTION_TYPE !< dim(nb_reactions) For each reaction, what is its type (cosmic ray evaporation, etc...)
integer, allocatable, dimension(:) :: RATE_FORMULA !< dim(nb_reactions) The index tracing the formula used for each specific 
!! reaction, defining its reaction rates in function of temperature and abundances.
integer, allocatable, dimension(:) :: REACTION_ID !< dim(nb_reactions) index of the reactions (one of the columns of the concerned file, 
!! declaring a given number for each reaction, like a hashtag. Some reactions may have the same ID though, when different 
!! temperature regimes are defined for instance

! Specific variables for first or second reactant of each reactions
integer, allocatable, dimension(:) :: reactant_1_idx !< dim(nb_reactions) Index of the first reactant species involved in the reaction
integer, allocatable, dimension(:) :: reactant_2_idx !< dim(nb_reactions) Index of the second reactant species involved in the reaction
real(double_precision), allocatable, dimension(:) :: DIFFUSION_RATE_1 !< dim(nb_reactions) [s-1] Diffusion rates used to compute the
!! grain reaction rate for reactant 1. It is equal to either the diffusion rate by thermal hopping or the diffusion rate by quantum 
!! tunneling.
real(double_precision), allocatable, dimension(:) :: DIFFUSION_RATE_2 !< dim(nb_reactions) [s-1] Diffusion rates used to compute the
!! grain reaction rate for reactant 2. It is equal to either the diffusion rate by thermal hopping or the diffusion rate by quantum 
!! tunneling.
real(double_precision), allocatable, dimension(:) :: EVAP_OVER_ACC_RATIO_1 !< dim(nb_reactions) EVAPORATION_RATES/ACCRETION_RATES
!! for reactant 1. This parameter is used in the modified rate subroutine.
real(double_precision), allocatable, dimension(:) :: EVAP_OVER_ACC_RATIO_2 !< dim(nb_reactions) EVAPORATION_RATES/ACCRETION_RATES 
!! for reactant 2. This parameter is used in the modified rate subroutine.

real(double_precision) :: initial_dtg_mass_ratio !< [no unit] initial dust to gas mass ratio
real(double_precision), allocatable, dimension(:) :: GTODN !< Gas to dust number ratio. 1/GTODN is equivalent to the grain abundance [no unit]
! real(double_precision), dimension(1:nb_grains) :: GTODN !< Gas to dust number ratio. 1/GTODN is equivalent to the grain abundance [no unit]
real(double_precision)                         :: GTODN_FIXED !< it is used for cases when we dont want to use grain size distribution
!integer                                        :: is_dust_MRN ! 1 =  distribution; 0 = WD distribution. (for specific case of WD see subroutine GRAIN_DIST(INDEX,DTYPE,A,DNDA))
real(double_precision) :: AV_NH_ratio !< Extinction over total hydrogen column density [mag/cm-2]
real(double_precision) :: grain_radius !< Grain radius [cm] 
real(double_precision) :: GRAIN_DENSITY !< grain density [g/cm^3]
real(double_precision) :: sticking_coeff_neutral  !< sticking coefficient for neutral  species on grain surface [no unit]
real(double_precision) :: sticking_coeff_positive !< sticking coefficient for positive species on grain surface [no unit]
real(double_precision) :: sticking_coeff_negative !< sticking coefficient for negative species on grain surface [no unit]
real(double_precision) :: MINIMUM_INITIAL_ABUNDANCE !< minimum value of the abundance (relative to H) [number ratio]
real(double_precision) :: initial_gas_density !< [part/cm^3] initial gas density of the structure
real(double_precision) :: initial_gas_temperature !< initial gas temperature [K], simulation parameter
real(double_precision) :: initial_dust_temperature !< initial dust temperature [K], simulation parameter
real(double_precision) :: INITIAL_VISUAL_EXTINCTION !< initial visual extinction [mag] 
real(double_precision) :: CR_IONISATION_RATE !< cosmic ray ionisation rate [s-1]
real(double_precision) :: UV_FLUX !< Scale factor for the UV flux, in unit of the reference flux (1.=nominal)
real(double_precision) :: DIFFUSION_BARRIER_THICKNESS !< [cm] thickness of the barrier that a surface species need to cross while 
!! undergoing quantum tunneling to diffuse from one surface site to another. This is used in the formalism by Hasegawa et al. (1992)
!! , see equation 10 of their paper (parameter a).
real(double_precision) :: SURFACE_SITE_DENSITY !< density of sites at the surface of each grain [cm-2]
real(double_precision), allocatable, dimension(:) :: nb_sites_per_grain !< Number of site per grain (site density * surface of the grain)
! real(double_precision), dimension(1:nb_grains) :: nb_sites_per_grain !< Number of site per grain (site density * surface of the grain)
real(double_precision) :: CHEMICAL_BARRIER_THICKNESS !< [cm] Parameter used to compute the probability for a surface reaction with 
!! activation energy to occur through quantum tunneling. This is the thickness of the energy barrier. See equation 6 from 
!! Hasegawa et al. (1992).
real(double_precision) :: CR_PEAK_GRAIN_TEMP !< Peak grain temperature when struck by a cosmic ray [K]
real(double_precision), dimension(:), allocatable  :: CR_PEAK_GRAIN_TEMP_all !< Peak grain temperature for all grains when struck by a cosmic ray [K]
real(double_precision) :: CR_PEAK_DURATION !< Peak duration [s] of CR_PEAK_GRAIN_TEMP
! real(double_precision) :: CR_PEAK_DURATION_r_dpnt !< Peak duration [s] of CR_PEAK_GRAIN_TEMP
real(double_precision) :: FE_IONISATION_RATE !< (cosmic) Fe-ion--grain encounter [s-1 grain-1] (for 0.1 micron grain) 
real(double_precision) :: FE_IONISATION_RATE_r_dpnt !< (cosmic) r dependent Fe-ion--grain encounter [s-1 grain-1] (for 0.1 micron grain) 
!! For cosmic photo desorptions, only Fe-ions are efficient to heat grains. 
real(double_precision) :: DIFF_BINDING_RATIO_SURF !< [no unit] Ratio used to compute the DIFFUSION_BARRIER from the BINDING_ENERGY if not known for surface species
real(double_precision) :: DIFF_BINDING_RATIO_MANT !< [no unit] Ratio used to compute the DIFFUSION_BARRIER from the BINDING_ENERGY if not known for mantle species
real(double_precision) :: START_TIME !< Start time of the simulation [s]
real(double_precision) :: STOP_TIME !< Stop time of the simulation [s]
real(double_precision) :: current_time !< Global current time of the simulation [s]
real(double_precision) :: VIB_TO_DISSIP_FREQ_RATIO !< [no unit] For the RRK (Rice Ramsperger-Kessel) desorption mechanism. Ratio of the vibration 
!! frequency (proper energy of a species when it is created on a grain) to the dissipation frequency (energy needed by the 
!! molecule to be evaporated from the grain surface). This ratio help to determine if a species evaporate after its formation 
!! on the grain surface. Since the dissipation frequency is usually unknown, this ratio is a free parameter. A common value is 1%.

! 1D variables
integer :: spatial_resolution = 1 !< sample in 1D dimension for the physical structure. If 1, then we are in 0D, else, we are in 1D
real(double_precision) :: grid_max_edge !< [AU] Maximum distance in 1D. We assume the minimum distance is always 0 AU
real(double_precision) :: grid_cell_size !< [cm] Grid cell size, or separation between two consecutive 1D points
real(double_precision), dimension(:), allocatable :: grid_sample !< dim(spatial_resolution) 1D sampling [cm] Must be linearly and equally spaced because of how the diffusion is treated.
real(double_precision), dimension(:), allocatable :: gas_temperature !< dim(spatial_resolution) current gas temperature [K]
real(double_precision), dimension(:), allocatable :: visual_extinction !< dim(spatial_resolution) visual extinction [mag] of the molecular cloud (or other astronomical object). 
!! It's the magnitude attenuation, difference from the absolute magnitude of the object and its apparent magnitude
real(double_precision), dimension(:,:), allocatable :: dust_temperature !< dim(nb_grains,spatial_resolution) current dust temperature [K]
real(double_precision), dimension(:), allocatable :: tmp_grain_temperature !< Grain temperature [K]
real(double_precision), dimension(:), allocatable :: H_number_density !< dim(spatial_resolution) [part/cm^3] Total H number density (both H and H2), representing the total gas density
real(double_precision), dimension(:), allocatable :: diffusion_coefficient !< dim(spatial_resolution) [cm^2/s] Diffusion coefficient for a 1D case
real(double_precision), dimension(:,:), allocatable :: GTODN_1D !< dim(spatial_resolution) [no unit] Gas to dust density number ratio read in the 1D_static.dat file, mostly for disk application
real(double_precision), dimension(:,:), allocatable :: GTODN_1D_temp  !< same as GTODN_1D but for temp storing data
real(double_precision), dimension(:), allocatable :: GTODN_0D !< dim(nb_grains) [no unit] Gas to dust density number ratio read in the 0D_grain_sizes.in file.
real(double_precision), dimension(:), allocatable :: GTODN_0D_temp  !< same as GTODN_0D but for temp storing data
real(double_precision), dimension(:), allocatable   :: AV_NH_1D !< dim(spatial_resolution) [??] conversion factor of AV to NH read the 1D_static.dat file, mostly for disk application (depends on the gas-to-dust mass ratio and the grain sizes)
real(double_precision), dimension(:,:), allocatable :: grain_radius_1D !< dim(spatial_resolution) [cm] grain sizes in 1D read the 1D_static.dat file, mostly for disk application (depends on the gas-to-dust mass ratio and the grain sizes)
real(double_precision), dimension(:,:), allocatable :: grain_radii_1D !< same as grain_radius_1D but for temp storing data
real(double_precision), dimension(:), allocatable   :: UV_flux_1D !< dim(spatial_resolution) UV flux factor in unit of reference flux for 1D structure. Usefull for YSOs.
real(double_precision), dimension(:,:), allocatable :: CR_PEAK_GRAIN_TEMP_all_1D !< for temp storing data
real(double_precision), dimension(:,:), allocatable :: grain_temp_1D          !< for temp storing data
! real(double_precision), dimension(:), allocatable :: grain_radius_1D !< dim(spatial_resolution) [cm] grain sizes in 1D read the 1D_static.dat file, mostly for disk application (depends on the gas-to-dust mass ratio and the grain sizes)
real(double_precision), dimension(:), allocatable :: NH_z !< dim(spatial_resolution) [cm^-2] H column density for a 1D case
real(double_precision), dimension(:), allocatable :: NH2_z !< dim(spatial_resolution) [cm^-2] H2 column density for a 1D case
real(double_precision), dimension(:), allocatable :: NCO_z !< dim(spatial_resolution) [cm^-2] CO column density for a 1D case
real(double_precision), dimension(:), allocatable :: NN2_z !< dim(spatial_resolution) [cm^-2] N2 column density for a 1D case
real(double_precision), dimension(:), allocatable :: NCH_z !< dim(spatial_resolution) [cm^-2] CH column density for a 1D case
real(double_precision), dimension(:), allocatable :: NCH3_z !< dim(spatial_resolution) [cm^-2] CH3 column density for a 1D case
real(double_precision), dimension(:), allocatable :: NCH4_z !< dim(spatial_resolution) [cm^-2] CH4 column density for a 1D case
real(double_precision), dimension(:), allocatable :: NH2O_z !< dim(spatial_resolution) [cm^-2] H2O column density for a 1D case
real(double_precision), dimension(:), allocatable :: NCO2_z !< dim(spatial_resolution) [cm^-2] CO2 column density for a 1D case
real(double_precision), dimension(:), allocatable :: NN2O_z !< dim(spatial_resolution) [cm^-2] N2O column density for a 1D case
real(double_precision), dimension(:), allocatable :: NOH_z !< dim(spatial_resolution) [cm^-2] OH column density for a 1D case
real(double_precision), dimension(:), allocatable :: NHCO_z !< dim(spatial_resolution) [cm^-2] HCO column density for a 1D case
real(double_precision), dimension(:), allocatable :: NH2CO_z !< dim(spatial_resolution) [cm^-2] H2 column density for a 1D case
real(double_precision), dimension(:), allocatable :: NCN_z !< dim(spatial_resolution) [cm^-2] CN column density for a 1D case
real(double_precision), dimension(:), allocatable :: NHCN_z !< dim(spatial_resolution) [cm^-2] HCN column density for a 1D case
real(double_precision), dimension(:), allocatable :: NHNC_z !< dim(spatial_resolution) [cm^-2] HNC column density for a 1D case
real(double_precision), dimension(:), allocatable :: NNH_z !< dim(spatial_resolution) [cm^-2] NH column density for a 1D case
real(double_precision), dimension(:), allocatable :: NNH2_z !< dim(spatial_resolution) [cm^-2] NH2 column density for a 1D case
real(double_precision), dimension(:), allocatable :: NNH3_z !< dim(spatial_resolution) [cm^-2] NH3 column density for a 1D case


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Sacha Gavino
!
!> @date 2020
!
! DESCRIPTION: 
!> Variables for photorates computation using cross-sections
!! some global size are set (nb_species_for_gas, nb_gas_phase_reactions, nb_species_for_gas, nb_surface_reactions, nb_species, nb_reactions)\n
!! This routine prepare allocation of global dynamical arrays
!
real(kind=8), dimension(nb_line_table_flux,2) :: table_ISRF  !
real(kind=8), dimension(nb_line_table_flux,2) :: table_Istar
real(kind=8), dimension(nb_line_table_flux,92) :: local_flux_dust ! 2 dimension table. (number of wavelength, number of spatial points)
real(kind=8), dimension(nb_line_table_flux,92) :: flux_isrf_dust ! 2 dimension table. (number of wavelength, number of spatial points)
real(kind=8), dimension(nb_line_table_flux) :: local_flux
real(kind=8), dimension(nb_line_table_flux) :: mol_opacity ! Table of opacity due to the main molecules as a function of wavelength
character(len=11) :: spec_photo_select
integer :: x_i
real(double_precision), dimension(nb_line_photorates) :: BR_photo !< dim(nbr lines in the file) Branching ratios for the photo reactions for which a detailed computation is done for disks
character(len=11), dimension(nb_line_photorates,8) :: REACTION_COMPOUNDS_NAMES_photo !< dim(MAX_COMPOUNDS,nbr lines in the file). Empty string means no species
character(len=11), dimension(nb_line_spec_photorates) :: spec_photo !<
real(kind=8), dimension(nb_line_photorates) :: table_photo_final !
real(kind=8), dimension(nb_line_table_flux,nb_line_spec_photorates,4) :: cross_sections
integer :: id_photo_H2, id_photo_CO, id_photo_N2 
integer :: id_photo_CO2, id_photo_H2O, id_photo_N2O
integer :: id_photo_CH, id_photo_CH3, id_photo_CH4 
integer :: id_photo_OH, id_photo_HCO, id_photo_H2CO
integer :: id_photo_NH, id_photo_NH2, id_photo_NH3
integer :: id_photo_CN, id_photo_HCN, id_photo_HNC

TYPE photodiss_type
    character(len=11) :: name
    real(8) :: k_diss
    real(8) :: k_ion
    ENDTYPE photodiss_type
TYPE (photodiss_type) :: photorates(nb_line_spec_photorates)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



integer, parameter :: MAX_NUMBER_REACTION_TYPE=100 !< Max number of various reaction type
! The following arrays start at 0 because the index correspond to the reaction type as indexed elsewhere, and there is a type 0 for reactions.
integer, dimension(0:MAX_NUMBER_REACTION_TYPE-1) :: type_id_start !< list of id start for each reaction type given their type number
integer, dimension(0:MAX_NUMBER_REACTION_TYPE-1) :: type_id_stop !< list of id stop for each reaction type given their type number

character(len=80) :: GRAIN_TEMPERATURE_TYPE = 'fixed' !< ('gas', 'computed', 'table_evolv', 'table_1D', 'fixed') How the grain temperature is computed in the code


procedure(get_grain_temperature_interface), pointer :: get_grain_temperature !< Pointer toward the routine that will calculate grain temperature

abstract interface 
  subroutine get_grain_temperature_interface(space,time, gas_temperature, av, grain_temperature)
  import
  
  implicit none

  ! Inputs
  integer, intent(in) :: space !<[in] current spatial point in 1D
  real(double_precision), intent(in) :: time !<[in] current time of the simulation [s]
  real(double_precision), intent(in) :: gas_temperature !<[in] gas temperature [K]
  real(double_precision), intent(in) :: av !<[in] visual extinction [mag]
  
  ! Outputs
  real(double_precision),dimension(:), intent(out) :: grain_temperature !<[out] grain temperature [K]
  !------------------------------------------------------------------------------
  
  end subroutine get_grain_temperature_interface
end interface

!--- Flags in file parameters.in ---
integer :: is_3_phase !< Flag for 3 phase model. 1:Activated, 0:2 phase model
integer :: IS_TEST = 1!< 1:True, 0:False. To thoroughly test reactions and stuff before running the code. This can take a few seconds
!! so you might switch it off when running several simulations with the exact same chemical network.
!integer :: is_dust_1D !< ! Reading the grain abundance and the NH/AV factor in the 1D_static.dat file (mostly for disks)
integer :: photo_disk !< ! photodissociation rates of molecules in the gas phase computation in the photorates_disks.f90 routines.
integer :: is_h2_formation_rate !< ! H2 formation rates on grain surfaces with the fluctuation effects. (Bron et al. (2014))
integer :: height_h2formation !< ! spatial point above which the H2 formation rates following B14 are computed (see documentation).
integer :: IS_GRAIN_REACTIONS !< Accretion, grain surface reactions
integer :: IS_H2_ADHOC_FORM !< Ad hoc formation of H2 on grain surfaces (1=activated)
integer :: GRAIN_TUNNELING_DIFFUSION !< How grain tunneling diffusion is handled
integer :: CONSERVATION_TYPE !< 0=only e- conserved; 1=elem #1 conserved, 2=elem #1 & #2, etc
integer :: MODIFY_RATE_FLAG !< Modify rates flag ; 1=modify H; 2=modify H,H2, 3=modify all, -1=H+H only
integer :: is_absorption_h2 !< H2 self-shielding from Lee & Herbst (1996) (1=activated)
integer :: is_absorption_co !< CO self-shielding. (1: Lee & Herbst (1996), 2: Visser et al. (2009)
integer :: is_absorption_n2 !< N2 self-shielding from Li et al. (2013) (1=activated)
integer :: is_photodesorb !< photodesorption processes flag (desactivated = 0)
integer :: is_crid !< CRID (cosmic rays induced diffusion) mechanism flag (desactivated = 0)
integer :: is_er_cir !< Eley-Rideal and complex induced reaction mechanisms flag (desactivated = 0)
integer :: is_reac_diff = 1 !< Flag for the reaction-diffusion competition
integer :: is_chem_des = 0 !< Flag for chemical desorption (0: Garrod 2007, 1: Minissale et al. 2016)

! About IS_STRUCTURE_EVOLUTION, describing the evolution of the physical structure properties with time
integer :: IS_STRUCTURE_EVOLUTION = 0 !< if 1, physical structure properties evolve with time. They come from structure_evolution.dat file, containing
!! {time [Myr], Av [mag], number density [part/cm3], gas temperature [K] and possibly grain temperature [K]} for the structure

procedure(get_structure_properties_interface), pointer :: get_structure_properties

abstract interface 
  subroutine get_structure_properties_interface(time, Av, density, gas_temperature)
  import 
  
  implicit none
  
  ! Inputs
  real(double_precision), intent(in) :: time !<[in] Current time of the simulation [s]
  
  ! Outputs
  real(double_precision), intent(out) :: Av !<[out] Visual extinction [mag]
  real(double_precision), intent(out) :: gas_temperature !<[out] gas temperature [K]
  real(double_precision), intent(out) :: density !<[out] gas density [part/cm^3]
  
  end subroutine get_structure_properties_interface
end interface


procedure(get_timestep_interface), pointer :: get_timestep !< Pointer that compute the next integration sub-step 
!! (withing an output step) depending on the situation (1D or not, and the type of 1D)

abstract interface 
  subroutine get_timestep_interface(current_time, final_time, next_timestep)
  import 
  
  implicit none
  
  ! Inputs
  real(double_precision), intent(in) :: current_time !<[in] current time [s]
  real(double_precision), intent(in) :: final_time !<[in] Final output time of the current 
  !! loop occurence. The last sub-step must lead exactly to this time [s]
  
  ! Outputs
  real(double_precision), intent(out) :: next_timestep !<[out] The next integration sub timestep withing an output integration step [s]
  
  end subroutine get_timestep_interface
end interface

character(len=80) :: STRUCTURE_TYPE = '0D' !< (0D, 1D_diff, 1D_no_diff)

procedure(a_structure_diffusion_interface), pointer :: structure_diffusion !< Pointer that diffuse the structure 
!! (mainly abundances) for a given timestep

abstract interface 
  subroutine a_structure_diffusion_interface(timestep, temp_abundances)
  import 
  
  implicit none
    
  ! Inputs
  real(double_precision), intent(in) :: timestep !<[in] timestep for the diffusion process [s]
  
  ! Inputs/Outputs
  real(double_precision), dimension(:,:), intent(inout) :: temp_abundances !<[in,out] dim(nb_species, spatial_resolution) 
  !! The abundances for all species, and 
  !! all 1D mesh points (relative to H) [number ratio]
  
  end subroutine a_structure_diffusion_interface
end interface

! About the optimization so that, for each species, we only check the reactions we know the species is involved.
integer :: max_reactions_same_species !< Maximum number of reactions in which any species will be involved. Used to set the array 'relevant_reactions'
integer, dimension(:,:), allocatable :: relevant_reactions !< dim(max_reactions_same_species, nb_species) For each species, store the list of reactions involving the species. 
!! When 0 are encountered for a given species, this means that no more reactions involve it. dimensions : (max_reactions_same_species, nb_species)
integer, dimension(:), allocatable :: nb_reactions_using_species !< dim(nb_species) For each species, the number of reactions in which it is used

! For LSODES
integer :: lrw !< declared length of RWORK (in user dimension).
integer :: liw !< declared length of IWORK (in user dimension).
integer, dimension(:), allocatable :: IWORK !< dim(liw) integer work array of length at least 30.
real(double_precision), dimension(:), allocatable :: RWORK !< dim(lrw) real work array of length at least:
!!\n             20 + 16*NEQ            for MF = 10,
!!\n             20 + (2 + 1./LENRAT)*NNZ + (11 + 9./LENRAT)*NEQ
!!\n                                    for MF = 121 or 222,
!!\n          where:
!!\n          NNZ    = the number of nonzero elements in the sparse
!!\n                   Jacobian (if this is unknown, use an estimate), and
!!\n          LENRAT = the real to integer wordlength ratio (usually 1 in
!!\n                   single precision and 2 in double precision).
!!\n          In any case, the required size of RWORK cannot generally
!!\n          be predicted in advance if MF = 121 or 222, and the value
!!\n          above is a rough estimate of a crude lower bound.  Some
!!\n          experimentation with this size may be necessary.
!!\n          (When known, the correct required length is an optional
!!\n          output, available in IWORK(17).)
integer :: nb_nonzeros_values !< number of non-zeros values in the jacobian. This is usefull for ODEPACK, to increase speed

! Diffusion and 1D variables
real(double_precision) :: X_IONISATION_RATE !< Ionisation rate due to X-rays [s-1]
real(double_precision) :: NH  ! column density [cm-2] (for the self shielding)
real(double_precision) :: NH2 ! column density [cm-2] (for the self shielding)
real(double_precision) :: NCO ! column density [cm-2] (for the self shielding)
real(double_precision) :: NN2 ! column density [cm-2] (for the self shielding)
real(double_precision) :: NH2O ! column density [cm-2] (for the self shielding)
real(double_precision) :: NCO2 ! column density [cm-2] (for the self shielding)
real(double_precision) :: NN2O ! column density [cm-2] (for the self shielding)
real(double_precision) :: NCH ! column density [cm-2] (for the self shielding)
real(double_precision) :: NCH3 ! column density [cm-2] (for the self shielding)
real(double_precision) :: NCH4 ! column density [cm-2] (for the self shielding)
real(double_precision) :: NOH ! column density [cm-2] (for the self shielding)
real(double_precision) :: NHCO ! column density [cm-2] (for the self shielding)
real(double_precision) :: NH2CO ! column density [cm-2] (for the self shielding)
real(double_precision) :: NCN ! column density [cm-2] (for the self shielding)
real(double_precision) :: NHCN ! column density [cm-2] (for the self shielding)
real(double_precision) :: NHNC ! column density [cm-2] (for the self shielding)
real(double_precision) :: NNH ! column density [cm-2] (for the self shielding)
real(double_precision) :: NNH2 ! column density [cm-2] (for the self shielding)
real(double_precision) :: NNH3 ! column density [cm-2] (for the self shielding)

logical :: first_step_done = .false. !< do we have currently done the first step of the integrator ?
integer :: NB_OUTPUTS !< Total number of outputs in the simulation
character(len=80) :: OUTPUT_TYPE !< Type of output times sampling. linear, log
!! linear: Output times are linearly spaced\n 
!! log   : Outputs times are log-spaced 
!! Ignored in the case of time evolving physical structure.


! Variables to be used by routines to access 'actual' values of differents physical properties. 
!! This might be not very clear but the thing is : get_temporal_derivatives and get_jacobian have fixed input/output because they are
!! joined as external routines in dlsodes. We can't modify and add inputs in these routines unless we change routines in the ODEPACK
!! module (which is obviously not possible). Thus, we must access these values as global variable, even if this imply changing theses
!! values inside a loop, which is not a proper thing to do in normal, correctly written code. 
real(double_precision) :: actual_gas_temp !< Gas temperature [K]
real(double_precision), dimension(:), allocatable :: actual_dust_temp !< dim(nb_grains) Dust temperature [K]
real(double_precision) :: actual_av !< Visual extinction [mag]
real(double_precision) :: actual_gas_density !< Gas density [part/cm^3]


contains 

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Christophe Cossou
!
!> @date 2014
!
! DESCRIPTION: 
!> @brief Routine to determine array sizes, namely number of reactions, 
!! of species, for gas, grain and in total. 
!! some global size are set (nb_species_for_gas, nb_gas_phase_reactions, nb_species_for_gas, nb_surface_reactions, nb_species, nb_reactions)\n
!! This routine prepare allocation of global dynamical arrays
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!  modified by Wasim Iqbal on 22-02-2017 to add grain size distribution
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine get_array_sizes()

implicit none

integer            :: i
character(len=11)  :: name_species(MAX_COMPOUNDS)
character(len=200) :: line
character(len=80)  :: line_format !< format of one line of gas_reaction.in or grain_reaction.in

integer            :: t_nb_grains_in_gas,t_nb_gas_grain_reactions_in_gas
integer            :: t_nb_species_for_grain !< total number of species involved in grain surface reactions in grain size distribution model
integer            :: t_nb_surface_reactions !< number of reactions on the grain surface in grain size distribution model
!these varaibles are needed as there are few gas phase reactions in grain_reactions.in files
!and also some gas_species in grain_species.in file.

character(len=1), parameter :: comment_character = '!' !< character that will indicate that the rest of the line is a comment
integer :: comment_position !< the index of the comment character on the line. if zero, there is none on the current string
integer :: error !< to store the state of a read instruction

! We get the number of reactions and species
call get_linenumber(filename='gas_species.in', nb_lines=nb_species_for_gas)
call get_linenumber(filename='gas_reactions.in', nb_lines=nb_gas_phase_reactions)

call get_linenumber(filename='grain_species.in', nb_lines=nb_species_for_grain)
call get_linenumber(filename='grain_reactions.in', nb_lines=nb_surface_reactions)

!-----------------------------------------------------------------------------------

open(903, file='gas_species.in', status='old')
open(904, file='gas_reactions.in', status='old')

name_species = ''
t_nb_grains_in_gas=0
t_nb_gas_grain_reactions_in_gas=0

  do
    read(903, '(a)', iostat=error)line
    if (error /= 0) exit
      
    ! We get only what is on the left of an eventual comment parameter
    comment_position = index(line, comment_character)
    
    ! if there are comments on the current line, we get rid of them
    if (comment_position.ne.0) then
      line = line(1:comment_position - 1)
    endif
    
    if (line.ne.'') then
       read(line, '(a)')  name_species(1)    
    endif
    
    
    if ((name_species(1)(1:5).eq.'GRAIN')) then
      t_nb_grains_in_gas = t_nb_grains_in_gas + 1
    endif  
  enddo  

name_species = ''

! Definition of the line format, common to gas_reaction.in and grain_reaction.in
  write(line_format, '(a,i0,a,i0,a)') '(', MAX_REACTANTS, 'A11,1x,', MAX_PRODUCTS, 'A11)' 
  do
    read(904, '(a)', iostat=error) line
    if (error /= 0) exit
      
    ! We get only what is on the left of an eventual comment parameter
    comment_position = index(line, comment_character)
    ! if there are comments on the current line, we get rid of them
    
    if (comment_position.ne.0) then
      line = line(1:comment_position - 1)
    end if
    
    if (line.ne.'') then
      read(line, line_format)  (name_species(I),I=1,MAX_COMPOUNDS)
        
      if ((name_species(1)(1:5).eq.'GRAIN') .or. (name_species(2)(1:5).eq.'GRAIN') .or. &
          (name_species(3)(1:5).eq.'GRAIN') .or. (name_species(4)(1:5).eq.'GRAIN') .or. &
          (name_species(5)(1:5).eq.'GRAIN') .or. (name_species(6)(1:5).eq.'GRAIN') .or. &
          (name_species(7)(1:5).eq.'GRAIN') .or. (name_species(8)(1:5).eq.'GRAIN')) then
        t_nb_gas_grain_reactions_in_gas = t_nb_gas_grain_reactions_in_gas + 1  
      endif
    endif
  enddo  
  
!    write(*,*)"---------------------------------------------------------------"
!    write(*,*)t_nb_species_for_grain, t_nb_surface_reactions
nb_species_for_gas = nb_species_for_gas + t_nb_grains_in_gas * (nb_grains-1)
nb_gas_phase_reactions = nb_gas_phase_reactions + t_nb_gas_grain_reactions_in_gas * (nb_grains-1)

close(903)
close(904)

  
open(901, file='grain_species.in', status='old')
open(902, file='grain_reactions.in', status='old')

name_species = ''
t_nb_species_for_grain=0
t_nb_surface_reactions=0
  do
    read(901, '(a)', iostat=error)line
    if (error /= 0) exit
      
    ! We get only what is on the left of an eventual comment parameter
    comment_position = index(line, comment_character)
    
    ! if there are comments on the current line, we get rid of them
    if (comment_position.ne.0) then
      line = line(1:comment_position - 1)
    endif
    
    if (line.ne.'') then
       read(line, '(a)')  name_species(1)    
    endif
    
    
    if ((name_species(1)(1:1).eq.'J').or.(name_species(1)(1:1).eq.'K')) then
      t_nb_species_for_grain=t_nb_species_for_grain + 1
    endif      
  enddo  
  
name_species = ''

! Definition of the line format, common to gas_reaction.in and grain_reaction.in
  write(line_format, '(a,i0,a,i0,a)') '(', MAX_REACTANTS, 'A11,1x,', MAX_PRODUCTS, 'A11)' 
  do
    read(902, '(a)', iostat=error) line
    if (error /= 0) exit
      
    ! We get only what is on the left of an eventual comment parameter
    comment_position = index(line, comment_character)
    ! if there are comments on the current line, we get rid of them
    
    if (comment_position.ne.0) then
      line = line(1:comment_position - 1)
    end if
    
    if (line.ne.'') then
      read(line, line_format)  (name_species(I),I=1,MAX_COMPOUNDS)
        
      if ((name_species(1)(1:1).eq.'J').or.(name_species(1)(1:1).eq.'K') .or. &
          (name_species(2)(1:1).eq.'J').or.(name_species(2)(1:1).eq.'K') .or. &
          (name_species(3)(1:1).eq.'J').or.(name_species(3)(1:1).eq.'K') .or. &
          (name_species(4)(1:1).eq.'J').or.(name_species(4)(1:1).eq.'K') .or. &
          (name_species(5)(1:1).eq.'J').or.(name_species(5)(1:1).eq.'K') .or. &
          (name_species(6)(1:1).eq.'J').or.(name_species(6)(1:1).eq.'K') .or. &
          (name_species(7)(1:1).eq.'J').or.(name_species(7)(1:1).eq.'K') .or. &
          (name_species(8)(1:1).eq.'J').or.(name_species(8)(1:1).eq.'K') ) then
        t_nb_surface_reactions=t_nb_surface_reactions+1
      endif  
    endif
  enddo  
  
!    write(*,*)"---------------------------------------------------------------"
!    write(*,*)t_nb_species_for_grain, t_nb_surface_reactions
nb_species_for_grain = nb_species_for_grain + t_nb_species_for_grain * (nb_grains-1)
nb_surface_reactions = nb_surface_reactions + t_nb_surface_reactions * (nb_grains-1)
close(901)
close(902)
!-----------------------------------------------------------------------------------

nb_species = nb_species_for_gas + nb_species_for_grain ! The total number of species, sum of species in gas and grain
nb_reactions = nb_gas_phase_reactions + nb_surface_reactions ! The total number of reactions, sum of species in gas and grain

! write(*,*)"nb_species_for_gas",nb_species_for_gas,"nb_gas_phase_reactions",nb_gas_phase_reactions
! write(*,*)"nb_species_for_grain",nb_species_for_grain,"nb_surface_reactions",nb_surface_reactions
! write(*,*)"nb_species,nb_reactions",nb_species,nb_reactions
!   pause
end subroutine get_array_sizes


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
!> @author 
!> Christophe Cossou
!
!> @date 2014
!
! DESCRIPTION: 
!> @brief Routine that allocate global arrays once their sizes are set
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine initialize_global_arrays()

implicit none

! Arrays associated with species
allocate(species_name(nb_species))
species_name(1:nb_species) = ''

allocate(abundances(nb_species, spatial_resolution))
abundances(1:nb_species, 1:spatial_resolution) = 0.d0

allocate(SPECIES_MASS(nb_species))
SPECIES_MASS(1:nb_species) = 0.d0

allocate(STICK_SPEC(nb_species))
STICK_SPEC(1:nb_species) = 0.d0

allocate(THERMAL_HOPING_RATE(nb_species))
THERMAL_HOPING_RATE(1:nb_species) = 0.d0

allocate(CR_HOPING_RATE(nb_species))
CR_HOPING_RATE(1:nb_species) = 0.d0

allocate(ACCRETION_RATES(nb_species))
ACCRETION_RATES(1:nb_species) = 0.d0

allocate(EVAPORATION_RATES(nb_species))
EVAPORATION_RATES(1:nb_species) = 0.d0

allocate(EVAPORATION_RATES_TEMPO(nb_species))
EVAPORATION_RATES_TEMPO(1:nb_species) = 0.d0

allocate(BINDING_ENERGY(nb_species))
BINDING_ENERGY(1:nb_species) = 0.d0

allocate(DIFFUSION_BARRIER(nb_species))
DIFFUSION_BARRIER(1:nb_species) = 0.d0

allocate(GAP_ENERGY_BANDS(nb_species))
GAP_ENERGY_BANDS(1:nb_species) = 0.d0

allocate(FORMATION_ENTHALPY(nb_species))
FORMATION_ENTHALPY(1:nb_species) = 0.d0

allocate(VIBRATION_FREQUENCY(nb_species))
VIBRATION_FREQUENCY(1:nb_species) = 0.d0

allocate(ACC_RATES_PREFACTOR(nb_reactions))
ACC_RATES_PREFACTOR(1:nb_reactions) = 0.d0

allocate(TUNNELING_RATE_TYPE_1(nb_species))
TUNNELING_RATE_TYPE_1(1:nb_species) = 0.d0

allocate(TUNNELING_RATE_TYPE_2(nb_species))
TUNNELING_RATE_TYPE_2(1:nb_species) = 0.d0

allocate(SPECIES_CHARGE(nb_species))
SPECIES_CHARGE(1:nb_species) = 0

allocate(nb_reactions_using_species(nb_species))
nb_reactions_using_species(1:nb_species) = 0

! Variables associated with reactions
allocate(branching_ratio(nb_reactions))
branching_ratio(1:nb_reactions) = 0.d0

allocate(RATE_A(nb_reactions))
RATE_A(1:nb_reactions) = 0.d0

allocate(RATE_B(nb_reactions))
RATE_B(1:nb_reactions) = 0.d0

allocate(RATE_C(nb_reactions))
RATE_C(1:nb_reactions) = 0.d0

allocate(reaction_rates(nb_reactions))
reaction_rates(1:nb_reactions) = 0.d0

allocate(GRAIN_RANK(nb_reactions))
GRAIN_RANK(1:nb_reactions) = 0

allocate(reaction_rates_1D(spatial_resolution,nb_reactions))
reaction_rates_1D(1:spatial_resolution,1:nb_reactions) = 0.d0

allocate(DIFFUSION_RATE_1(nb_reactions))
DIFFUSION_RATE_1(1:nb_reactions) = 0.d0

allocate(DIFFUSION_RATE_2(nb_reactions))
DIFFUSION_RATE_2(1:nb_reactions) = 0.d0

allocate(EVAP_OVER_ACC_RATIO_1(nb_reactions))
EVAP_OVER_ACC_RATIO_1(1:nb_reactions) = 0.d0

allocate(EVAP_OVER_ACC_RATIO_2(nb_reactions))
EVAP_OVER_ACC_RATIO_2(1:nb_reactions) = 0.d0

allocate(ACTIVATION_ENERGY(nb_reactions))
ACTIVATION_ENERGY(1:nb_reactions) = 0.d0

allocate(REACTION_TMIN(nb_reactions))
REACTION_TMIN(1:nb_reactions) = 0.d0

allocate(REACTION_TMAX(nb_reactions))
REACTION_TMAX(1:nb_reactions) = 0.d0

allocate(SURF_REACT_PROBA(nb_reactions))
SURF_REACT_PROBA(1:nb_reactions) = 0.d0

allocate(REACTION_TYPE(nb_reactions))
REACTION_TYPE(1:nb_reactions) = 0

allocate(reactant_1_idx(nb_reactions))
reactant_1_idx(1:nb_reactions) = 0

allocate(reactant_2_idx(nb_reactions))
reactant_2_idx(1:nb_reactions) = 0

allocate(RATE_FORMULA(nb_reactions))
RATE_FORMULA(1:nb_reactions) = 0

allocate(REACTION_ID(nb_reactions))
REACTION_ID(1:nb_reactions) = 0

allocate(REACTION_COMPOUNDS_ID(MAX_COMPOUNDS,nb_reactions))
REACTION_COMPOUNDS_ID(1:MAX_COMPOUNDS,1:nb_reactions) = 0

allocate(REACTION_COMPOUNDS_NAMES(MAX_COMPOUNDS,nb_reactions))
REACTION_COMPOUNDS_NAMES(1:MAX_COMPOUNDS,1:nb_reactions) = ''

! 1D arrays
allocate(grid_sample(spatial_resolution))
grid_sample(1:spatial_resolution) = 0.d0

allocate(gas_temperature(spatial_resolution))
gas_temperature(1:spatial_resolution) = 0.d0

allocate(visual_extinction(spatial_resolution))
visual_extinction(1:spatial_resolution) = 0.d0

allocate(dust_temperature(nb_grains,spatial_resolution))
dust_temperature(1:nb_grains,1:spatial_resolution) = 0.d0

allocate(H_number_density(spatial_resolution))
H_number_density(1:spatial_resolution) = 0.d0

allocate(diffusion_coefficient(spatial_resolution))
diffusion_coefficient(1:spatial_resolution) = 0.d0

allocate(GTODN_1D(nb_grains,spatial_resolution))
GTODN_1D(1:nb_grains,1:spatial_resolution) = 0.d0

allocate(GTODN_0D(nb_grains))
GTODN_0D(1:nb_grains) = 0.d0

allocate(AV_NH_1D(spatial_resolution))
AV_NH_1D(1:spatial_resolution) = 0.d0

allocate(UV_flux_1D(spatial_resolution))
UV_flux_1D(1:spatial_resolution) = 0.d0

allocate(grain_radius_1D(nb_grains,spatial_resolution))
grain_radius_1D(1:nb_grains,1:spatial_resolution) = 0.d0

allocate(NH_z(spatial_resolution))
NH_z(1:spatial_resolution) = 0.d0

allocate(NH2_z(spatial_resolution))
NH2_z(1:spatial_resolution) = 0.d0

allocate(NN2_z(spatial_resolution))
NN2_z(1:spatial_resolution) = 0.d0

allocate(NCO_z(spatial_resolution))
NCO_z(1:spatial_resolution) = 0.d0

allocate(NH2O_z(spatial_resolution))
NH2O_z(1:spatial_resolution) = 0.d0

allocate(NCO2_z(spatial_resolution))
NCO2_z(1:spatial_resolution) = 0.d0

allocate(NN2O_z(spatial_resolution))
NN2O_z(1:spatial_resolution) = 0.d0

allocate(NCH_z(spatial_resolution))
NCH_z(1:spatial_resolution) = 0.d0

allocate(NCH3_z(spatial_resolution))
NCH3_z(1:spatial_resolution) = 0.d0

allocate(NCH4_z(spatial_resolution))
NCH4_z(1:spatial_resolution) = 0.d0

allocate(NOH_z(spatial_resolution))
NOH_z(1:spatial_resolution) = 0.d0

allocate(NHCO_z(spatial_resolution))
NHCO_z(1:spatial_resolution) = 0.d0

allocate(NH2CO_z(spatial_resolution))
NH2CO_z(1:spatial_resolution) = 0.d0

allocate(NCN_z(spatial_resolution))
NCN_z(1:spatial_resolution) = 0.d0

allocate(NHCN_z(spatial_resolution))
NHCN_z(1:spatial_resolution) = 0.d0

allocate(NHNC_z(spatial_resolution))
NHNC_z(1:spatial_resolution) = 0.d0

allocate(NNH_z(spatial_resolution))
NNH_z(1:spatial_resolution) = 0.d0

allocate(NNH2_z(spatial_resolution))
NNH2_z(1:spatial_resolution) = 0.d0

allocate(NNH3_z(spatial_resolution))
NNH3_z(1:spatial_resolution) = 0.d0


! Prime elements
allocate(INITIAL_ELEMENTAL_ABUNDANCE(NB_PRIME_ELEMENTS))
INITIAL_ELEMENTAL_ABUNDANCE(1:NB_PRIME_ELEMENTS) = 0.d0

allocate(PRIME_ELEMENT_IDX(NB_PRIME_ELEMENTS))
PRIME_ELEMENT_IDX(1:NB_PRIME_ELEMENTS) = 0

allocate(species_composition(NB_PRIME_ELEMENTS,nb_species))
species_composition(1:NB_PRIME_ELEMENTS,nb_species) = 0


allocate(rate_tot_des(nb_grains))
rate_tot_des(1:nb_grains) = 0.d0
allocate(rate_tot(nb_grains))  
rate_tot(1:nb_grains) = 0.d0
allocate(rate_tot_acc(nb_grains))
rate_tot_acc(1:nb_grains) = 0.d0

end subroutine initialize_global_arrays



end module global_variables
