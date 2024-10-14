!******************************************************************************
! MODULE: numerical_types
!******************************************************************************
!
! DESCRIPTION: 
!> @brief Module that define simple and double precision parameters. 
!!\n\n
!!\n real(sp) :: a ! real simple precision
!!\n real(dp) :: b ! real double precision
!
!******************************************************************************

module numerical_types

  implicit none
  
  integer, parameter :: simple_precision = selected_real_kind(6,37) !< simple precision parameter (used with real(sp) :: arg)
  integer, parameter :: double_precision = selected_real_kind(15,307) !< double precision parameter (used with real(dp) :: arg)

end module numerical_types
