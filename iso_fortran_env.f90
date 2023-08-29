!******************************************************************************
! MODULE: iso_fortran_env
!******************************************************************************
!
! DESCRIPTION: 
!> @brief Module that define various global parameters \n\n
!!
!! Nonintrinsic version for Lahey/Fujitsu Fortran for Linux. 
!! See Subclause 13.8.2 of the Fortran 2003 standard. 
!
!******************************************************************************

module iso_fortran_env

  implicit none 
  public 

  integer, parameter :: Character_Storage_Size = 8 
  integer, parameter :: Error_Unit = 0 !< unit to print in the standard error, to be used : write(Error_Unit, *) "Error 404 bla bla bla"
  integer, parameter :: File_Storage_Size = 8 
  integer, parameter :: Input_Unit = 5  !< Input unit integer
  integer, parameter :: IOSTAT_END = -1 
  integer, parameter :: IOSTAT_EOR = -2 
  integer, parameter :: Numeric_Storage_Size = 32 
  integer, parameter :: Output_Unit = 6 !< unit to print in the standard output to be used as : write(Output_Unit, *) "Hello World !"

end module iso_fortran_env