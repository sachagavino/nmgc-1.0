!******************************************************************************
! PROGRAM: unitary_tests
!******************************************************************************
!
! DESCRIPTION: 
!> @brief Program that test different functions implemented in nautilus. 
!! Especially routines that were added afterwards, and are not  present in 
!! the original version, thus can't be compared between the two versions.
!!\n\n
!! This program need a test simulation in the sub-folder 'tests'. Output results
!! will be stored in the same folder.
!!\n\n
!! This binary is conceived so that it is run by the python script unitary_tests.py
!
!******************************************************************************
! 

program unitary_tests

  use numerical_types
  use iso_fortran_env, only : error_unit
  use global_variables
  use nautilus_main

    
  implicit none
  
  call initialisation()
  
  call test_read_time_evolution()
  call test_time_interpolation()
  
  contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Christophe Cossou
!
!> @date 7 may 2014
!
! DESCRIPTION: 
!> @brief Test the interpolation of time evolution properties of the structure
!! then compare them to the raw profile, just to be sure the interpolation is
!! correct.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  subroutine test_time_interpolation()
  
    implicit none
    
    integer, parameter :: nb_sample = 100
    real(double_precision), parameter :: t_min = 0.1d0 ! in AU
    real(double_precision), parameter :: t_max = 100.d0! in AU
    real(double_precision), parameter :: step = (t_max / t_min)**(1.d0 / (dfloat(nb_sample)))
    
    real(double_precision) :: time ! Time in Myr
    real(double_precision) :: av, density, gas_temperature, grain_temperature

    logical :: isDefined !< true or false, to test if some files exists or not, for instance
    
    integer :: j ! for loops
    
    write(*,*) 'Test of the time evolution interpolation'
    
    inquire(file='structure_evolution.dat', exist=isDefined)
    
    if (.not.isDefined) then
      write(Error_unit, *) "Warning: The file 'structure_evolution.dat' doesn't exist."
      write(Error_unit, *) "A default one is being generated"
      
      open(10, file='structure_evolution.dat')
      write(10,'(a)') '! time    log(Av)    log(n)    log(T)   '
      write(10,'(a)') '! (Myr)   (mag)      (cm-3)    (K)      '
      write(10,'(a)') '0.000e+00 -1.231e+00 1.813e+00 1.698e+00'
      write(10,'(a)') '5.000e+00 -1.221e+00 1.760e+00 1.715e+00'
      write(10,'(a)') '1.000e+01 -1.511e+00 1.327e+00 1.946e+00'
      write(10,'(a)') '1.500e+01 -1.617e+00 1.166e+00 2.063e+00'
      write(10,'(a)') '2.000e+01 -1.906e+00 7.300e-01 2.463e+00'
      write(10,'(a)') '2.500e+01 -1.781e+00 8.800e-01 2.340e+00'
      write(10,'(a)') '3.000e+01 -1.448e+00 1.367e+00 1.930e+00'
      write(10,'(a)') '3.500e+01 -1.091e+00 1.909e+00 1.656e+00'
      write(10,'(a)') '4.000e+01 -1.202e+00 1.748e+00 1.719e+00'
      write(10,'(a)') '4.500e+01 -1.307e+00 1.583e+00 1.805e+00'
      write(10,'(a)') '5.000e+01 -7.650e-01 2.432e+00 1.480e+00'
      write(10,'(a)') '5.300e+01 -1.024e+00 2.058e+00 1.594e+00'
      close(10)
      
    end if
    
    call init_structure_evolution()
    
    open(10, file='test_time_interpolation.dat')
    write(10,*) '# time, av, density, gas_temperature, grain_temperature'
    do j=1,nb_sample
      time = t_min * step**(j - 1.d0)
      ! We generate cartesian coordinate for the given Semi-major axis
      
      call get_structure_properties_table(time=time*(1e6*YEAR),&!Inputs
        Av=av, density=density, gas_temperature=gas_temperature)
      
      
      write(10,*) time, av, density, gas_temperature
    end do
    close(10)
    
    ! We create associated gnuplot files
    open(10, file="av_interpolation.gnuplot")
    open(11, file="density_interpolation.gnuplot")
    open(12, file="gas_temperature_interpolation.gnuplot")
    
    do j=10,12
      write(j,*) "set terminal pdfcairo enhanced"
    end do
    
    write(10,*) '!rm "av_interpolation.pdf"'
    write(10,*) "set output 'av_interpolation.pdf'"
    write(10,*) 'set ylabel "Visual extinction [mag]"'
    
    write(11,*) '!rm "density_interpolation.pdf"'
    write(11,*) "set output 'density_interpolation.pdf'"
    write(11,*) 'set ylabel "Density [part/cm^3]"'
    write(11,*) 'set logscale y'
    
    write(12,*) '!rm "gas_temperature_interpolation.pdf"'
    write(12,*) "set output 'gas_temperature_interpolation.pdf'"
    write(12,*) 'set ylabel "Temperature [K]"'
    write(12,*) 'set logscale y'
    
    do j=10,12
      write(j,*) 'set xlabel "Time [Myr]"'
      write(j,*) 'set grid'
      write(j,*) 'set logscale x'
    end do

    write(10,*) 'plot "test_time_interpolation.dat" using 1:2 with points linetype 1 pointtype 2 title "Interpolation",\'
    write(10,*) '     "structure_evolution.dat" using 1:(10**$2) with lines title "Profile"'

    write(11,*) 'plot "test_time_interpolation.dat" using 1:3 with points linetype 1 pointtype 2 title "Interpolation",\'
    write(11,*) '     "structure_evolution.dat" using 1:(10**$3) with lines title "Profile"'

    write(12,*) 'plot "test_time_interpolation.dat" using 1:4 with points linetype 1 pointtype 2 title "Interpolation",\'
    write(12,*) '     "structure_evolution.dat" using 1:(10**$4) with lines title "Profile"'

    close(10)
    close(11)
    close(12)
  
  end subroutine test_time_interpolation

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Christophe Cossou
!
!> @date 13 may 2014
!
! DESCRIPTION: 
!> @brief Test the reading of time evolution, compare the actual profile 
!! with the read one (in global variable)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  subroutine test_read_time_evolution()
  
    implicit none
    
    logical :: isDefined !< true or false, to test if some files exists or not, for instance
    
    integer :: j ! for loops
    
    write(*,*) 'Test of read profile from time evolution data file'
    
    inquire(file='structure_evolution.dat', exist=isDefined)
    
    if (.not.isDefined) then
      write(Error_unit, *) "Warning: The file 'structure_evolution.dat' doesn't exist."
      write(Error_unit, *) "A default one is being generated"
      
      open(10, file='structure_evolution.dat')
      write(10,'(a)') '! time    log(Av)    log(n)    log(T)   '
      write(10,'(a)') '! (Myr)   (mag)      (cm-3)    (K)      '
      write(10,'(a)') '0.000e+00 -1.231e+00 1.813e+00 1.698e+00'
      write(10,'(a)') '5.000e+00 -1.221e+00 1.760e+00 1.715e+00'
      write(10,'(a)') '1.000e+01 -1.511e+00 1.327e+00 1.946e+00'
      write(10,'(a)') '1.500e+01 -1.617e+00 1.166e+00 2.063e+00'
      write(10,'(a)') '2.000e+01 -1.906e+00 7.300e-01 2.463e+00'
      write(10,'(a)') '2.500e+01 -1.781e+00 8.800e-01 2.340e+00'
      write(10,'(a)') '3.000e+01 -1.448e+00 1.367e+00 1.930e+00'
      write(10,'(a)') '3.500e+01 -1.091e+00 1.909e+00 1.656e+00'
      write(10,'(a)') '4.000e+01 -1.202e+00 1.748e+00 1.719e+00'
      write(10,'(a)') '4.500e+01 -1.307e+00 1.583e+00 1.805e+00'
      write(10,'(a)') '5.000e+01 -7.650e-01 2.432e+00 1.480e+00'
      write(10,'(a)') '5.500e+01 -1.024e+00 2.058e+00 1.594e+00'
      close(10)
      
    end if
    
    call init_structure_evolution()
    
    open(10, file='test_time_read.dat')
    write(10,*) '# time (Myr), av [mag], density [part/cm^3], gas_temperature (K), grain_temperature (K)'
    do j=1,structure_sample
      
      if (read_dust) then
        write(10,*) structure_time(j)/(1e6*YEAR), 10.d0**(structure_log_Av(j)), 10.d0**(structure_log_density(j)), &
                    10.d0**(structure_log_gas_temperature(j)), 10.d0**(structure_log_dust_temperature(j))
      else
        write(10,*) structure_time(j)/(1e6*YEAR), 10.d0**(structure_log_Av(j)), 10.d0**(structure_log_density(j)), &
                    10.d0**(structure_log_gas_temperature(j))
      endif
      ! if the grain temperature is not defined in the structure_evolution.dat file, grain temperature will appear to be 1K 
      !! since the log will be set to 0
    end do
    close(10)
    
    ! We create associated gnuplot files
    open(10, file="test_av_read.gnuplot")
    open(11, file="test_density_read.gnuplot")
    open(12, file="test_gas_temperature_read.gnuplot")
    
    do j=10,12
      write(j,*) "set terminal pdfcairo enhanced"
    end do
    
    write(10,*) '!rm "test_av_read.pdf"'
    write(10,*) "set output 'test_av_read.pdf'"
    write(10,*) 'set ylabel "Visual extinction [mag]"'
    
    write(11,*) '!rm "test_density_read.pdf"'
    write(11,*) "set output 'test_density_read.pdf'"
    write(11,*) 'set ylabel "Density [part/cm^3]"'
    write(11,*) 'set logscale y'
    
    write(12,*) '!rm "test_gas_temperature_read.pdf"'
    write(12,*) "set output 'test_gas_temperature_read.pdf'"
    write(12,*) 'set ylabel "Temperature [K]"'
    write(12,*) 'set logscale y'
    
    do j=10,12
      write(j,*) 'set xlabel "Time [Myr]"'
      write(j,*) 'set grid'
      write(j,*) 'set logscale x'
    end do

    write(10,*) 'plot "test_time_read.dat" using 1:2 with points linetype 1 pointtype 2 title "Read profile",\'
    write(10,*) '     "structure_evolution.dat" using 1:(10**$2) with lines title "Profile"'

    write(11,*) 'plot "test_time_read.dat" using 1:3 with points linetype 1 pointtype 2 title "Read profile",\'
    write(11,*) '     "structure_evolution.dat" using 1:(10**$3) with lines title "Profile"'

    write(12,*) 'plot "test_time_read.dat" using 1:4 with points linetype 1 pointtype 2 title "Read profile",\'
    write(12,*) '     "structure_evolution.dat" using 1:(10**$4) with lines title "Profile"'
    
    close(10)
    close(11)
    close(12)
    
    if (read_dust) then
      open(13, file="test_grain_temperature_read.gnuplot")
      
      write(13,*) "set terminal pdfcairo enhanced"
      write(13,*) '!rm "test_grain_temperature_read.pdf"'
      write(13,*) "set output 'test_grain_temperature_read.pdf'"
      write(13,*) 'set ylabel "Temperature [K]"'
      write(13,*) 'set logscale y'
      write(13,*) 'set xlabel "Time [Myr]"'
      write(13,*) 'set grid'
      write(13,*) 'set logscale x'
      
      select case(GRAIN_TEMPERATURE_TYPE)
      case('gas')
        write(13,*) 'plot "test_time_read.dat" using 1:5 with points linetype 1 pointtype 2 title "Read profile",\'
        write(13,*) '     "structure_evolution.dat" using 1:(10**$4) with lines title "Profile"'
      
      case('table')
        write(13,*) 'plot "test_time_read.dat" using 1:5 with points linetype 1 pointtype 2 title "Read profile",\'
        write(13,*) '     "structure_evolution.dat" using 1:(10**$5) with lines title "Profile"'

      case default
        write(13,*) 'plot "test_time_read.dat" using 1:5 with points linetype 1 pointtype 2 title "Read profile"'
      end select
      
      close(13)
    endif
    
  
  end subroutine test_read_time_evolution
end program unitary_tests
