!******************************************************************************
! MODULE: utilities
!******************************************************************************
!
! DESCRIPTION: 
!> @brief Module that contains all the routines linked to generic features such
!! as getting the number of lines in a files etc...
!
!******************************************************************************

module utilities

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
!> @brief Routine to retrieve the number of lines of a given file whose
!! filename is passed as an argument. Commented lines are excluded
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine get_linenumber(filename, nb_lines, opt_comment)

  implicit none
  
  ! Input
  character(len=*), intent(in) :: filename !< [in] the filename of the file we want the number of lines
  character(len=1), intent(in), optional :: opt_comment !< [in] optional input parameter to define a 
  !! commenting character different from the default one '!'
  
  ! Output
  integer, intent(out) :: nb_lines !< [out] the number of line of the input file
  
  ! Local
  integer :: error
  logical test
  character(len=80) :: line
  character(len=1) :: comment_character = '!' !< character that will indicate that the rest of the line is a comment
  integer :: comment_position !< the index of the comment character on the line. if zero, there is none on the current string

  ! Treating options
  if (present(opt_comment)) then
    comment_character = opt_comment
  else 
    ! Force definition to '!' because the 'save' option might keep the non default value for other run of the routine
    comment_character = '!'
  endif

  !------------------------------------------------------------------------------
  nb_lines = 0
  
  ! Read in filenames and check for duplicate filenames
  inquire (file=filename, exist=test)
  if (.not.test) then
    write(Error_Unit,'(a,a,a)') 'Error: the file "',trim(filename),'" does not exist.'
    call exit(1)
  end if
  open(15, file=filename, status='old')
  
  error = 0
  do 
    read(15,'(a)',iostat=error) line
    
    if(error.ne.0) then
      exit
    endif
    
    ! We get only what is on the left of an eventual comment parameter
    comment_position = index(line, comment_character)
    
    ! if there are comments on the current line, we get rid of them
    if (comment_position.ne.0) then
      line = line(1:comment_position - 1)
    end if
    
    
    
    if (trim(line).ne.'') then
      nb_lines = nb_lines + 1
    endif
  enddo
  close(15)
  
end subroutine get_linenumber

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Christophe Cossou
!
!> @date 2014
!
! DESCRIPTION: 
!> @brief For a given line passed as an argument, return the number of columns on this line
!! i.e, the number of tokens separated by spaces.
!
!> @return number of columns (separated by spaces) on the input string
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
integer function get_nb_columns(filename)
implicit none

! Inputs
character(len=*), intent(in) :: filename !< [in] datafile name we want to get the number of columns

! Locals
integer :: i, n, toks
character(len=20000) :: line
character(len=1), parameter :: comment_character = '!' !< character that will indicate that the rest of the line is a comment
integer :: comment_position !< the index of the comment character on the line. if zero, there is none on the current string
integer :: error !< to store the state of a read instruction

logical :: isDefined

inquire(file=filename, exist=isDefined)
if (isDefined) then

  open(10, file=filename, status='old')
  
  ! We get the first data line (without taking into account commented or blanck lines)
  line = ''
  do while(line.eq.'')
    read(10, '(a)', iostat=error) line
    
    if (error /= 0) exit
      
    ! We get only what is on the left of an eventual comment parameter
    comment_position = index(line, comment_character)
    
    ! if there are comments on the current line, we get rid of them
    if (comment_position.ne.0) then
      line = line(1:comment_position - 1)
    end if
  enddo
  close(10)
else
  write(Error_unit,*) 'Error: The file ', trim(filename),' does not exist.'
  call exit(1)
endif

i = 1
n = len_trim(line)
toks = 0
get_nb_columns = 0
do while(i <= n)
   do while(line(i:i) == ' ') 
     i = i + 1
     if (n < i) return
   enddo
   toks = toks + 1
   get_nb_columns = toks
   do
     i = i + 1
     if (n < i) return
     if (line(i:i) == ' ') exit
   enddo
enddo
end function get_nb_columns 

end module utilities