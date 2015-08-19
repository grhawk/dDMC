MODULE read_input
  ! Variable To read
  !  inputtagfile
  !  inputcoofile
  !  atomdatafile
  !  debugflag
  !  dftype
  !  tagtype -> Define the format of the file containing the atomic charges/population

  USE precision
  USE string_tools
  IMPLICIT NONE
  character(kch) :: inputtagfile,inputcoofile,atomdatafile,debugflag,dfparamsfile,tagtype,readgradflag,method
  integer(ki) :: dftype = 0
  logical :: dfprint

  PRIVATE
  PUBLIC :: inputtagfile,inputcoofile,atomdatafile,debugflag,dftype,dfprint,dfparamsfile,tagtype,readgradflag,method
  PUBLIC :: read_stdin, starting_program_announce
  
CONTAINS
  
  SUBROUTINE read_stdin
    IMPLICIT NONE
    character(kch) :: junk,nojunk(2)
    integer(ki) :: err,i,nstring,ninput,nl,dummy(1)

    nl = 0; ninput = 0
    readfile: do
       read(5,'(A64)',iostat=err) junk
!       print*, junk  ! debug
       if( err > 0 ) call die('ERROR: while reading stdout')
       if( err < 0 ) exit
       nl = nl + 1
       
       if( junk(1:1) /= '#' .and. junk(1:1) /= '' ) then
          
          call shrink_string(junk,'=',nojunk,nstring)
!          write(0,* ) junk,nojunk,nstring ! debug
          if( nstring /= 2 ) call die('ERROR: checkpoint 1 in read_input')
          
          select case (trim(adjustl(nojunk(1))))
          case ('tag') 
             inputtagfile = trim(adjustl(nojunk(2)))
             ninput = ninput + 1
          case ('geometry') 
             inputcoofile = trim(adjustl(nojunk(2)))
             ninput = ninput + 1
          case ('atomdata')
             atomdatafile = trim(adjustl(nojunk(2)))
             ninput = ninput + 1
          case ('method')
             method = trim(adjustl(nojunk(2)))
             ninput = ninput + 1
          case ('debugflag')
             debugflag = trim(adjustl(nojunk(2)))
             ninput = ninput + 1
          case ('dfparameters')
             dfparamsfile = trim(adjustl(nojunk(2)))
             ninput = ninput + 1
          case ('dftype')
             dummy = char2int(nojunk(2),1)
             dftype = dummy(1)
             ninput = ninput + 1
          case ('tagtype')
             tagtype = trim(adjustl(nojunk(2)))
          case ('gradient')
             readgradflag = trim(adjustl(nojunk(2)))
          case ('dfprint')
             dfprint = .true.
             dfparamsfile = ''
             exit readfile
          case default
            write(0,*) 'INPUT ERROR IN LINE ',nl
            stop 1
          end select
          
       end if
    end do readfile
    
!    if( ninput /= 6 .and. .not. dfprint ) call die('ERROR: wrong parameters in input')
    
!    call starting_program_announce
    
  END SUBROUTINE read_stdin
  
  SUBROUTINE starting_program_announce
    IMPLICIT NONE
    
    write(0,*) 'Program is starting whit this data:'
    write(0,*) 'file tag: #', inputtagfile,"#"
    write(0,*) 'file geometry: ', inputcoofile
    write(0,*) 'atomic data: ', atomdatafile
    write(0,*) 'dftype: ',dftype
    write(0,*) 'tagtype: ',tagtype
    write(0,*) 'dfparameters: ',dfparamsfile
    write(0,*) 'gradient: ',readgradflag
    write(0,*) 'debug: ',debugflag
    
  END SUBROUTINE starting_program_announce

END MODULE read_input


! # -> Comment
! tag = character
! geometry = character
! atomdata = character
! b0 = real
