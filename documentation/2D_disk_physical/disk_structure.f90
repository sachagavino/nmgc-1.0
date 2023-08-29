PROGRAM disk_structure

implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! definition of variables
!
! temp_z  vertical temperature profile
! H2dens_z=dblarr(nb_points)      ; vertical density profile
! exp_H2dens_z=dblarr(nb_points)  ; exponential of vertical density profile (intermediate value in the calculation)
! Av_z=dblarr(nb_points)          ; vertical visual extinction profile
! zspace=dblarr(nb_points)        ; spatial points

integer, parameter :: nb_points = 64
integer, parameter :: nb_radius = 14

real(kind=8) :: ref_radius, central_mass, tmidplan_ref 
real(kind=8) :: Tatmos_ref, UV_ref, q, sigma, surface_density_ref
real(kind=8) :: surface_density_exponent, small_grains, big_grains
real(kind=8) :: dtogm_up, transition_altitude, Tmidplan, Tatmos
real(kind=8) :: H, int, dens_midplan, Omega2, maximum

real(kind=8), dimension (1:nb_points) :: temp_z, H2dens_z, exp_H2dens_z 
real(kind=8), dimension (1:nb_points) :: Av_z, zspace, GTODN, grain_radius
real(kind=8), dimension (1:nb_points) :: AV_NH

real(kind=8), dimension (1:nb_radius) :: radius_from_center_table
real(kind=8) :: radius_from_center, dtogm_midplan, H2dens_midplan, H_atmos
real(kind=8) :: mass1, mass2, UV_factor, He_ab

character(len=80) :: aa
character(len=11) :: aa2
character(len=9) :: file_name
character(len=9), dimension (1:nb_radius) :: file_table
character(len=20) :: file_name2
character(len=20), dimension (1:nb_radius) :: file_table2
integer :: j, i, iradius
real(kind=8), dimension (1:14) :: inputs
real(kind=8) :: int2

! defining some constants

real(kind=8), parameter :: KB = 1.38054D-16     ! boltzmann constant in erg.K-1
real(kind=8), parameter :: au_conv = 1.49597D13   ! conversion factor of AU in cm
real(kind=8), parameter :: mean_molecular_weight = 2.4D0    
real(kind=8), parameter :: amu = 1.66043D-24      ! atomic mass unit in g
real(kind=8), parameter :: grav_const = 6.668D-8  ! gravitational constant in cm3g-1s-2
real(kind=8), parameter :: NH_to_AV_conversion = 1.6d+21  ! conversion factor of H colmun density to Av (Wagenblast \& Hartquist 1989)
real(kind=8), parameter :: dtogm = 1.D-2                  ! standard dust to gas mass ratio
real(kind=8), parameter :: grain_density = 3.D0           ! volumic density of each grain (in g/cm3)
real(kind=8), parameter :: Pi = 3.1415927D0


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Putting a few things to zero
!

Tmidplan=0.D0                   ! mid plan temperature
Tatmos=0.D0                     ! temperature in the atmosphere where z=4H
H=0.D0                          ! scale height at defined by the mid-plan temperature
int=0.D0
dens_midplan=0.D0
Omega2=0.D0
maximum=0.D0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Reading input file

open(unit = 100, file = 'source_parameters.in', status = 'old', action = 'read')
 
do j = 1, 4
 read(100,*) aa
enddo

do j= 1,14
 read(100,'(29x,e9.3)') int2
 inputs(j)=int2
enddo

close(100)

ref_radius = inputs(1)
central_mass = inputs(2)
Tmidplan_ref = inputs(3)
Tatmos_ref = inputs(4)
UV_ref = inputs(5)
q = inputs(6)
sigma = inputs(7)
surface_density_ref = inputs(8)
surface_density_exponent = inputs(9)
small_grains = inputs(10)
big_grains = inputs(11)
dtogm_up = inputs(12)
transition_altitude = inputs(13)
He_ab = inputs(14)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Defining the radius where we want to compute the structure

radius_from_center_table(1) = 1.D0
radius_from_center_table(2) = 2.D0
radius_from_center_table(3) = 3.D0
radius_from_center_table(4) = 4.D0
radius_from_center_table(5) = 5.D0
radius_from_center_table(6) = 10.D0
radius_from_center_table(7) = 20.D0
radius_from_center_table(8) = 30.D0
radius_from_center_table(9) = 40.D0
radius_from_center_table(10) = 60.D0
radius_from_center_table(11) = 100.D0
radius_from_center_table(12) = 200.D0
radius_from_center_table(13) = 300.D0
radius_from_center_table(14) = 400.D0

file_table(1) = '001AU.txt'
file_table(2) = '002AU.txt'
file_table(3) = '003AU.txt'
file_table(4) = '004AU.txt'
file_table(5) = '005AU.txt'
file_table(6) = '010AU.txt'
file_table(7) = '020AU.txt'
file_table(8) = '030AU.txt'
file_table(9) = '040AU.txt'
file_table(10) = '060AU.txt'
file_table(11) = '100AU.txt'
file_table(12) = '200AU.txt'
file_table(13) = '300AU.txt'
file_table(14) = '400AU.txt'

file_table2(1) = 'parameters_001AU.txt'
file_table2(2) = 'parameters_002AU.txt'
file_table2(3) = 'parameters_003AU.txt'
file_table2(4) = 'parameters_004AU.txt'
file_table2(5) = 'parameters_005AU.txt'
file_table2(6) = 'parameters_010AU.txt'
file_table2(7) = 'parameters_020AU.txt'
file_table2(8) = 'parameters_030AU.txt'
file_table2(9) = 'parameters_040AU.txt'
file_table2(10) = 'parameters_060AU.txt'
file_table2(11) = 'parameters_100AU.txt'
file_table2(12) = 'parameters_200AU.txt'
file_table2(13) = 'parameters_300AU.txt'
file_table2(14) = 'parameters_400AU.txt'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Recompute dust to gas massratio to remove He
! In the following, dtogm is used as a H/dust mass ratio

dtogm_up = dtogm_up*(1.d0+4.D0*He_ab)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Loop over all radius to compute the physical structure

do iradius = 1, nb_radius
!do iradius = 1,1
 radius_from_center = radius_from_center_table(iradius)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! computation of temperature at the given radius

Tmidplan = Tmidplan_ref*(radius_from_center/ref_radius)**(-q)
Tatmos = Tatmos_ref*(radius_from_center/ref_radius)**(-q)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  computation of H in cm (assuming the vertical static equilibrium)
! with Tmidplan

H=sqrt(KB*Tmidplan*(radius_from_center*au_conv)**3.D0/(mean_molecular_weight*amu*grav_const*central_mass))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Computation of H_atmos 
! H_atmos is the scale height corresponding to Tatmos

H_atmos=sqrt(KB*Tatmos*(radius_from_center*au_conv)**3.D0/(mean_molecular_weight*amu*grav_const*central_mass))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! computation of Omega**2 for the vertical gravity

Omega2 = grav_const*central_mass/(radius_from_center*au_conv)**3.0D0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! computation of the box (total height is 4H)
! formula used by FH:

DO i = 0,nb_points-1
	zspace (i+1) = (1.D0-2.D0*float(i)/(2.D0*float(nb_points)-1.D0))*4.D0*H
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Computation of the 1D vertical temperature profile
!
! using the definition by Williams & Best (2014)
! meaning that T is constant above 4 scale height defined
! by the temperature in the mid-plan of the disk. Below this
! height the temperature is computed using a sinus.

do i = 1,nb_points 
	int = Pi*zspace(i)/(2.D0*zspace(1))
	temp_z(i)=Tmidplan+(Tatmos-Tmidplan)*((sin(int))**(2.D0*sigma))
enddo	

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Computation of the 1D vertical density profile
! The density is computing assuming the hydrostatic equilibrium

! computation of the density in the midplan of disk:
H2dens_midplan = surface_density_ref*(radius_from_center/ref_radius)**(-surface_density_exponent)/&
(mean_molecular_weight*amu)/H/sqrt(2.D0*Pi)

!computation of the density (in cm-3) for each z point:

exp_H2dens_z(1)=0.d0      ! the most external density is set to zero as a first step
H2dens_z(1)=exp(exp_H2dens_z(1))
do i=2,nb_points
	exp_H2dens_z(i) = exp_H2dens_z(i-1) - (log(temp_z(i))-log(temp_z(i-1))) - Omega2 * &
	mean_molecular_weight*amu/(KB * temp_z(i)) * zspace(i)*(zspace(i)-zspace(i-1))
	H2dens_z(i)=exp(exp_H2dens_z(i))
enddo

! Rescaling so that the maximum density dens_midplan 
maximum=maxval(H2dens_z)
do i=1,nb_points  
 H2dens_z(i) = H2dens_z(i)/maximum*H2dens_midplan
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Computation of the 1D vertical visual extinction profile
! 
	
! It is assumed that the H column density above 4H is null. 
! This is different from the previous code but with some tests
! I found that the setting an external H column density or not
! would only change for the most external point. This has no impact on the 
! column densities.

AV_z(1)=0.D0
grain_radius(1)=small_grains
gtodn(1)=4.D0*Pi*grain_density*(grain_radius(1)**3)/(3.D0*dtogm_up*1.66053892d-24)
AV_NH(1)=1.D0/NH_to_AV_conversion*(dtogm_up / 1.d-2)
	
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute the mass conservation to compute the correct gtodn
! as a function of altitude

mass1=0.D0
mass2=0.D0

do i=2,nb_points 
 	if (zspace(i)/H.gt.transition_altitude) mass1=mass1+H2dens_z(i)*mean_molecular_weight*amu
	if (zspace(i)/H.le.transition_altitude) mass2=mass2+H2dens_z(i)*mean_molecular_weight*amu
enddo

dtogm_midplan=(dtogm*(mass1+mass2)-dtogm_up*mass1)/mass2
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! define the grain radii as a function of altitude

do i = 2,nb_points 
	if (zspace(i)/H.gt.transition_altitude) then 
		grain_radius(i)=small_grains
		gtodn(i)=4.D0*Pi*grain_density*(grain_radius(i)**3)/(3.D0*dtogm_up*1.66053892d-24)
		AV_NH(i)=1.D0/NH_to_AV_conversion*(dtogm_up / 1.d-2)*(1.e-5/grain_radius(i))
	endif	
    if (zspace(i)/H.le.transition_altitude) then 
    	grain_radius(i)=big_grains
    	gtodn(i)=4.D0*Pi*grain_density*(grain_radius(i)**3)/(3.D0*dtogm_midplan*1.66053892d-24)
    	AV_NH(i)=1.D0/NH_to_AV_conversion*(dtogm_midplan / 1.d-2)*(1.e-5/grain_radius(i))
    endif
    AV_z(i)=AV_z(i-1)+H2dens_z(i)*2.D0*AV_NH(i)*(zspace(i-1)-zspace(i))*(1.e-5/grain_radius(i))
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Computation of the UV irradiation factor at the requested radius
! the UV is divided by two because since we do not make a full 3D 
! transfer we assume that only half of the photons irradiated from the star 
! diffuse into the disk, the other half diffuse in the opposite direction.
! The UV flux coming from the star is assumed to have the same 
! spectrum as the interstellar DRAINE field multiplied by a certain amount and diluted as 1/r^2.

UV_factor=UV_ref/((radius_from_center/ref_radius)**2+(4.D0*H/ref_radius/au_conv)**2)/2.D0

	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! writing the results in a format that can be read by Nautilus
! 

file_name = file_table(iradius)

open (unit=10,file=file_name,action="write",status="new")

write(10,'(a,e9.3)') '! Reference radius AU =',ref_radius
write(10,'(a,e9.3)') '! Midplan temperature (K) at the reference radius AU =',Tmidplan_ref
write(10,'(a,e9.3)') '! Atmosphere (4H) temperature (K) at the reference radius AU =',Tatmos_ref
write(10,'(a,e9.3)') '! Radius where the structure is computed (AU) =',radius_from_center
write(10,'(a,e9.3)') '! Midplan temperature (K) at the requested radius =',Tmidplan
write(10,'(a,e9.3)') '! Atmosphere temperature (K) at the requested radius =',Tatmos
write(10,'(a,e9.3)') '! Mass of the central object (g) =',central_mass
write(10,'(a,e9.3)') '! UV field coming from the star at the reference radius AU (unit: mean ISM field) =',UV_ref
write(10,'(a,e9.3)') '! UV field coming from the star at the requested radius =',UV_factor
write(10,'(a,i3)') '! Vertical resolution = ',nb_points
write(10,'(a,e9.3)') '! Scale Height H in au = ',H/au_conv


write(10,'(a)') '! Distance [AU] ; H Gas density [part/cm^3] ; Gas Temperature [K] ; Visual Extinction [mag] ;& 
& Diffusion coefficient [cm^2/s]; Dust Temperature [K]; 1/abundance of grains ; AV/NH conversion factor ;& 
& radius of grains (cm)'

do i = 1,nb_points 
write(10,'(9e14.5)') zspace(i)/au_conv,H2dens_z(i)*2.D0,Temp_z(i),AV_z(i),0.D0,0.D0,gtodn(i),AV_NH(i),grain_radius(i)
enddo

close(unit=10)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! changing the UV flux in the parameters.in file for Nautilus
! 
file_name2 = file_table2(iradius)

open (unit=100, file='parameters.in', status = 'old', action = 'read')
open (unit=10,file=file_name2,action="write",status="new")

do i = 1,58 
	read(100,'(a80)') aa
	write(10,'(a80)') aa
enddo

read(100,'(a11)') aa2
write(10,'(a11,e9.3)') aa2,	UV_factor

do i = 1,38	
	read(100,'(a80)') aa
	write(10,'(a80)') aa
enddo

close(100)
close(10)


!print,';;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;'
!print,mass1,mass2,dtogm_midplan

ENDDO

end program disk_structure

