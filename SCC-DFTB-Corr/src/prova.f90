MODULE read_tag_dftbp

  USE precision
  USE string_tools
  USE file_tools

! The purpose of this module is reading tag files !
! provided by DFTB+. The foundamental subroutine  !
! is 'getdata'. This subroutine accept as dummy   !
! argument (property,file). 'property' is  the    !
! name  of a property present in the tag file,    !
! while 'file' is the name of tag file to examine.!
! Doing a call to 'getdata(property file)', you   !
! obtain in the main program a variable, namely   !
! 'matrix' of the form: real(kr),dimension(:,:,:) !
! and the variable ifrmt(3) that contains the     !
! dimensions of the variable matrix. The variable !
! matrix contains all the data associated at the  !
! choised property.                               !

  IMPLICIT NONE

  integer(ki),parameter :: maxrank=3

  character(maxchar),dimension(3)       :: form_shrk
  real(kr),allocatable,dimension(:,:,:) :: matrix
  integer(ki)  :: ifrmt(maxrank)

  PRIVATE
  PUBLIC  :: getdata,matrix 
  PUBLIC  :: ifrmt  ! debug 

CONTAINS
  
  SUBROUTINE read_property(name) !property_name
    IMPLICIT NONE
    character(*), intent(in) :: name !name of property (e.g. fillings)
    
    character(maxchar) :: form
    integer(ki) :: nword
    character(maxchar) :: junk
    integer(ki) :: err
    
    ! This loop looks for the property (name) into the tag file
    doread:do
       read(iu,'(A20,A20)',iostat=err) junk, form
       if (err /= 0) stop 'Property not found'
!       print*, '#',name,'#  ', '  #',junk,'#   ', '  #',form,'#' ! debug
       if ( junk .eq. name) then
          form = trim(adjustl(form))
          !  print*, form ! debug

          ! Call to determining the format of the property
          ! form_shrk(1) => variable type
          ! form_shrk(2) => rank
          ! form_shrk(3) => format
          call shrink_string(form,':',form_shrk,nword)
!          print*, nword,'#',form_shrk(:),'#'   ! debug
          ! select the right format and then read all the 
          ! values associated at the property 'name'.
          select case (trim(adjustl(form_shrk(1))))
          case ('real')
             select case (form_shrk(2))
             case ('0')
                print*, 'scalar' ! not implemented
             case ('1')
                print*, 'vector' ! not implemented
             case('2')
!                print*, 'matrix rank 2' ! debug
                ! read the data with the right format and
                ! create the variable matrix
                CALL read_matrix(form_shrk(3))
             case('3')
                print*, 'matrix rank 3' ! not implemented
             case default
                print*, 'Hai rotto il cazzo2' ! not implemented
             end select
          case ('logical')
             print*, 'logical' ! not implemented
          case ('integer')
             print*, 'integer' ! not implemented
          case default
             print*, 'Hai rotto il cazzo1' ! not implemented
          end select
          exit doread
       end if
    end do doread

  END SUBROUTINE read_property
  

  SUBROUTINE read_matrix(form_strg)
    ! read the data with the right format and
    ! create the variable matrix
    IMPLICIT NONE
    character(maxchar),intent(IN) :: form_strg
    character(maxchar) :: frmt(maxrank)

    integer(ki) :: junk
    integer(ki) :: i,j,k ! counter
    
    frmt(:) = '1'
!    print*,  frmt(:)
    ! Determining the matrix dimensions
    call shrink_string(form_strg,',',frmt,junk)

!    print*, '! debug'
    
    ! Converting the matrix dimension from string to number
    ifrmt = char2int(frmt,maxrank) 
!    print*, '! debug'
!    print*, 'frmt(1) ',ifrmt(1),' frmt(2) ',ifrmt(2)   ! debug

    ! Allocating the matrix
    allocate(matrix(ifrmt(1),ifrmt(2),ifrmt(3))) ! Need more options
    matrix = 0.d0
    
!    print*, iu  ! debug
    read(iu,*) (((matrix(i,j,k), i=1,ifrmt(1)), j=1,ifrmt(2)), k=1,ifrmt(3))

!!$    do i=1,ifrmt(1)                          ! debug
!!$       do j=1,ifrmt(2)                       ! debug
!!$       do k=1,ifrmt(3)                       ! debug
!!$          write(*,'(24f8.3)') matrix(i,j,k) ! debug
!!$       end do                                ! debug
!!$       end do                                ! debug
!!$    end do                                   ! debug

  END SUBROUTINE read_matrix

  
  SUBROUTINE getdata(property,file)
    ! This is the only subroutine "callable"
    ! from the main program. This subroutine
    ! (using the above ones) generate a variable
    ! named matrix_tag with the requested data.
    IMPLICIT NONE
    character(*),intent(IN) :: file,property
    
    call open_file(file)
    call read_property(property)
    call clsfile(file)
    
  END SUBROUTINE  getdata

END MODULE read_tag_dftbp
  

program prova
  USE read_tag_dftbp
  USE string_tools
  character*50 :: prop
 character*40 :: vect(50)  ! string_tools
  integer :: nword         ! string_tools

  print*,
  print*, '---------------------------------------------'

  call get_command_argument(1,prop)
!  print*, prop

!!!DEBUGGING: read_property
!  OPEN(5010,file='results.tag',form='FORMATTED',STATUS='OLD')

!  call open_file('results.tag')
!  call read_property(prop)
  
  call getdata(prop,'results.tag')
!  call getdata('eigenvalues','prova.f90')
    do i=1,ifrmt(1)                          ! debug
       do j=1,ifrmt(2)                       ! debug
       do k=1,ifrmt(3)                       ! debug
          write(*,'(24f8.3)') matrix(i,j,k) ! debug
       end do                                ! debug
       end do                                ! debug
    end do                                   ! debug

!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! DEBUGGING: string_tools
!!$  call shrink_string(prop,':',vect,nword)
!!$
!!$  print*, nword
!!$  print*, '#',(vect(i),i=1,nword),'#'
!!!!!!!!!!!!!!!!!!!!!!!!
end program prova
