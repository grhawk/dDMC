MODULE file_tools
  USE precision

  IMPLICIT NONE

  TYPE file2unit
     integer(ki) :: unit
     character*60 :: file
  END type file2unit

  integer(ki),parameter :: initial_value = 5000    ! First opened file unit - 1
  integer(ki),parameter :: maxchar = 60            ! Standard lenght for character variable
  integer(ki),parameter :: maxfile = 15            ! Max num

  type(file2unit),save,dimension(maxfile) :: file_unit 

  
  INTERFACE fiit
     MODULE PROCEDURE getfilename, getunit
  END INTERFACE fiit

  INTERFACE closefile
     MODULE PROCEDURE close_file, close_unit
  END INTERFACE closefile

  PRIVATE
  PUBLIC :: openfile,closefile,fiit



CONTAINS
  character(maxchar) FUNCTION getfilename(input)
    integer(ki),intent(IN) :: input
    
    integer(ki)  :: i    ! counter
    
    doloop:do i = 1,maxfile
       if ( input == file_unit(i)%unit) then
          getfilename = file_unit(i)%file
          exit doloop
       end if
    end do doloop
    
  END FUNCTION getfilename
  
  integer(ki) FUNCTION getunit(input)
    character(*),intent(IN) :: input
    
    integer(ki)  :: i    ! counter
    
    doloop:do i = 1,maxfile
       if ( input == file_unit(i)%file) then
          getunit = file_unit(i)%unit
          exit doloop
       end if
       if( i == maxfile ) stop 'ERROR: File not found in getunit'
    end do doloop
    
  END FUNCTION getunit
  
  SUBROUTINE openfile(file_name,readwrite)
    IMPLICIT NONE
    character(*),optional,intent(IN) :: readwrite
    character(*), intent(IN) :: file_name
    logical :: newold
    integer(ki) :: iu 
    integer(ki)              :: err

    call refresh_file_unit(unit=iu,file=file_name,create=.true.)

!    print*, 'opening ',file_name  ! debug
    
    select case (file_name)
    case('ciccio.brutto')
       ! Personal options for that filename
    case default
       if( .not. present(readwrite)) stop 'ERROR: You have to specify if the file is writable or not'
       if( readwrite == 'write' ) OPEN(iu,file=file_name,action='READWRITE',form='FORMATTED',status='NEW',IOSTAT=err)
       if( readwrite == 'read' ) OPEN(iu,file=file_name,action='READ',form='FORMATTED',STATUS='OLD',IOSTAT=err)
    end select
    
!    print*, iu,file_name,' opened' ! debug

!    do i = 1,10  ! debug
!       read(iu,*) junk  ! debug
!       print*, junk  ! debug
!    end do  ! debug
!    rewind(iu)  ! debug
    
!    print*, err ! debug
    if(err /=  0) call errors('opening',file_name)
    

 !   print*, opened_file ! debug
 !   print*, iu ! debug
  END SUBROUTINE openfile
  
  SUBROUTINE close_file(file_name)
    IMPLICIT NONE
    character(*),intent(IN) :: file_name
    integer(ki) :: unit
    integer(ki) :: i
    
    call refresh_file_unit(unit=unit,file=file_name,create=.false.)

    close(unit)
    

  END SUBROUTINE close_file

  SUBROUTINE close_unit(unit)
    IMPLICIT NONE
    integer (ki) :: unit
    
    call refresh_file_unit(unit,create=.false.)

    close(unit)
    
  END SUBROUTINE close_unit

  SUBROUTINE refresh_file_unit(unit,file,create)
    IMPLICIT NONE
    character(*),intent(IN),optional :: file
    logical,intent(IN) :: create
    integer(ki),optional,intent(INOUT) :: unit
    integer(ki),save         :: opened_files
    integer(ki) :: i


    if( opened_files < 1 ) then
       do i = 1,maxfile
          file_unit(i)%unit = initial_value + i
          file_unit(i)%file = ''
       end do
    end if



    if( .not. present(unit) .and. .not. present(file) ) stop 'ERROR in refreshing file_unit'
    if( .not.  create .and. present(file) )  unit = fiit(file)

    if( .not. create) then
       if(file_unit(unit - initial_value)%file == '' ) stop 'ERROR: File not found'
       file_unit(unit - initial_value)%file = ''
       opened_files = opened_files - 1
    else
       doloop: do i = 1,maxfile
          if( file_unit(i)%file == '' ) then
             file_unit(i)%file = file
             unit = file_unit(i)%unit
             opened_files = opened_files + 1
             exit doloop
          end if
          if( i == maxfile) STOP 'too much opened files'
       end do doloop
    end if
    
!    do i =1,maxfile   ! debug
!       print*, file_unit(i)%unit,' --- ',file_unit(i)%file  ! debug
!    end do  ! debug
!    print*, 'opened files ',opened_files, '   last unit ',unit  ! debug

  END SUBROUTINE refresh_file_unit

  SUBROUTINE errors(type,parameters)
    IMPLICIT NONE
    character(*), INTENT(IN) :: type,parameters
    
    select case (type)
       
    case('opening')
       write(0,10100) parameters
       stop

    end select

10100 format('Error while opening file: ',A30)   
    
  END SUBROUTINE errors
  
END MODULE file_tools

!!$program test
!!$  USE file_tools
!!$
!!$  IMPLICIT NONE
!!$  
!!$  print*, 'open cioacioa in write mode'
!!$  call openfile('cioacioa','write')
!!$  print*, 'open cazcaz in write mode'
!!$  call openfile('cazcaz','write')
!!$  
!!$
!!$  print*, 'write cioacioa'
!!$  write(fiit('cioacioa'),*) 'cioacioa'
!!$  print*, 'write cazcaz'
!!$  write(fiit('cazcaz'),*) 'cazcaz'
!!$
!!$
!!$  print*, 'closing cioacioa'
!!$  call closefile('cioacioa')
!!$  print*, 'opening results.tag'
!!$  call openfile('results.tag','read')
!!$  print*, 'closing cazcaz'
!!$  call closefile('cazcaz')
!!$  print*, 'closing results.tag'
!!$  call closefile('results.tag')
!!$
!!$  print*, 'closeded'
!!$
!!$  
!!$
!!$end program test
