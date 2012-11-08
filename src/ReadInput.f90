MODULE ReadInput
  ! Variable To read
  !  inputtagfile
  !  inputcoofile
  !  atomdatafile
  !  b0
  !  A
  !  df  (damping function) (TT,TTfix)

  USE precision
  USE string_tools
  IMPLICIT NONE
  character(kch) :: inputtagfile,inputcoofile,atomdatafile,DampFunc,debugflag
  real(kr)      :: b0,A

  PRIVATE
  PUBLIC :: inputtagfile,inputcoofile,atomdatafile,b0,A,read_stdin,DampFunc,debugflag
!  PUBLIC :: 
  
CONTAINS
  
  SUBROUTINE read_stdin
    IMPLICIT NONE
    character(kch) :: junk,nojunk(2)
    integer(ki) :: err,i,nstring,ninput,nl

    nl = 0; ninput = 0
    do
       read(5,'(A64)',iostat=err) junk
!       print*, junk  ! debug
       if( err > 0 ) stop 'ERROR: while reading stdout'
       if( err < 0 ) exit
       nl = nl + 1
       
       if( junk(1:1) /= '#' .and. junk(1:1) /= '' ) then
          
          call shrink_string(junk,'=',nojunk,nstring)
!          write(0,*) 'asd',junk,nojunk,nstring ! debug
          if( nstring /= 2 ) stop 'ERROR: checkpoint 1 in read_input'
          
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
          case ('b0')
             b0 =  char2real(nojunk(2))
             ninput = ninput + 1
          case ('A')
             A =  char2real(nojunk(2))
             ninput = ninput + 1
          case ('df')
             DampFunc = trim(adjustl(nojunk(2)))
             ninput = ninput + 1
          case ('debugflag')
             debugflag = trim(adjustl(nojunk(2)))
             ninput = ninput + 1
          case default
             write(0,*) 'INPUT ERROR IN LINE ',nl
          end select
          
       end if
    end do
    
    if( ninput /= 7 ) stop 'ERROR: wrong parameters in input'
    
!    call starting_program_announce
    
  END SUBROUTINE read_stdin
  
  SUBROUTINE starting_program_announce
    IMPLICIT NONE
    
    write(0,*) 'Program is starting whit this data:'
    write(0,*) 'file tag: #', inputtagfile,"#"
    write(0,*) 'file geometry: ', inputcoofile
    write(0,*) 'atomic data: ', atomdatafile
    write(0,*) 'b0 value: ', b0
    write(0,*) 'A  value: ', A
    write(0,*) 'Dumping Function: ',DampFunc
    
  END SUBROUTINE starting_program_announce

END MODULE ReadInput


! # -> Comment
! tag = character
! geometry = character
! atomdata = character
! b0 = real
