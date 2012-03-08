MODULE string_tools
  USE precision
  IMPLICIT NONE
  character(30),parameter :: module_name = 'string_tools'  ! For Error subroutine

  integer(ki),parameter :: maxchars=100  ! Max character per string
  integer(ki),parameter :: maxwords=100  ! Max word per string
  integer(ki),parameter :: maxcshrk=5    ! Max character in shrinker

  integer(ki) :: strlen

  PRIVATE
  PUBLIC ::  shrink_count, shrink_string, char2int, char2real!, real2char


CONTAINS
  
!!$  real(kr) FUNCTION real2char(charvar)
!!$    IMPLICIT NONE
!!$    character(*), intent(IN) :: charvar
!!$    
!!$    read(int2char,*) charvar
!!$    
!!$  END FUNCTION int2char

  FUNCTION char2int(charvar,n)
    ! This function converts numbers in string-form 
    ! in integer
    ! The dummy arguments are 
    ! charvar => String
    ! n       => if charvar is a vector, n is the dimension
    !            of that vector, n has to be an integer
    !            if the conversion is not for vector you can 
    !            put n=1 and this subroutine will work well for
    !            you (to test!!)
    IMPLICIT NONE
    integer(ki),intent(IN) :: n
    character(*),intent(IN) :: charvar(n)
    integer(kr) :: char2int(n)
    integer(kr) :: i

    do i = 1,n
       read(charvar(i),*) char2int(i)
    end do
    
  END FUNCTION char2int

  FUNCTION char2real(charvar) !(to test)
    ! Like the char2int but converts the numbers in real
    IMPLICIT NONE
    character(*),intent(IN) :: charvar
    real(kr) :: char2real
!    integer(kr) :: i
    
!    do i = 1,n
       read(charvar,*) char2real
!    end do
  END FUNCTION char2real
  

  SUBROUTINE shrink_count(stringin,cshrk,dims)
    !stringin  -> input string to shrink
    !cshrk     -> shrineker character 
    !dims      -> numer of occuerrences
    !
    ! This subroutine counts the number of the occurrence <cshrk> in
    ! the string <stringin> and return this value <dims>
    IMPLICIT NONE
    character(*),intent(IN)  :: stringin
    character(*),intent(IN)  :: cshrk
    integer(ki), intent(OUT) :: dims

    integer(ki) :: oldpw,pw

    integer(ki) :: i                         

    strlen = len(stringin)                                ! Not only a Control

    i = 1;oldpw = 0; pw = 0                             
    do
       if(i >= maxwords) call errors('exceeds_words')   ! Verify capability of array
       pw = scan(stringin(oldpw+1:strlen),cshrk)        ! Each time one occurrence is found
                                                        ! move the cursor after this occurrence.
                                                        ! pw store all the occurrences position.

       if(oldpw == pw+oldpw) exit                       ! Stop loop when occurrences are finished
       pw = pw + oldpw                                  !+This two lines are necessary to skip the
       oldpw = pw                                       !+just controlled string
       i = i + 1                                        ! counter
    end do
    dims = i - 1                                        ! Since the first value of i is 1 instead 0

  END SUBROUTINE shrink_count


  SUBROUTINE shrink_string(stringin,cshrk,stringout,nword)
    !stringin -> input string to shrink
    !cshrk     -> shrineker character 
    !stringout -> output array with the shrinked string
    !
    ! This routine take a string and a 'separator' as input 
    ! and give a vector with the 'word' and an integer 
    ! indicating the number of the words. The words are 
    ! separated using the separator as delimitng character.
    IMPLICIT NONE
    character(*)             :: stringin    ! Intent avoid to trim the string
    character(*),intent(IN)  :: cshrk
    character(*),intent(OUT) :: stringout(maxwords)
!    character(100) :: asd
    integer(ki)   :: oldpw,occn,start,nword
    integer(ki),dimension(0:maxwords-1)   :: pw   !Position of cshrk characters

    integer(ki) :: i,j  ! counter

!    stringin = ':a1d:a2d:::asdah:'    ! debug

    strlen = len_trim(stringin)
    !print*, 'stringin ',stringin ! debug

    ! Determine all cshrk positions => See subroutine shrink_count
    i = 1;oldpw = 0; pw = 0                             
    do
       if(i >= maxwords) call errors('exceeds_words')   
       pw(i) = scan(stringin(oldpw+1:strlen),cshrk)  
       !print*, i,pw(i),oldpw,pw(i)+oldpw debug
       if(oldpw == pw(i)+oldpw) exit                    
       pw(i) = pw(i) + oldpw                            
       oldpw = pw(i)                                    
       i = i + 1
    end do
    occn = i - 1                                ! REMEMBER: The total number of occurrences is
                                                !           i - 1
    !print*,'occn', occn   ! debug
    if (occn == 0) then                         ! If there is no occurrences, we have only one word, 
       stringout(1) = stringin                  ! so the stringout should contain only this word, the 
       nword = 1                                ! same one of the stringin
       return
    end if

!!$    start = 1;i=1                               ! If the first position is not an occurrence...
!!$    if (pw(1) > 1) start=0                      ! this two lines are needed
!!$    do j = start,occn
!!$       if( pw(j) /= pw(j+1)-1) then             ! Avoid a empty position in stringout if 
!!$                                                ! two occurrence are near
!!$          i = i + 1
!!$          print*, pw(i)
!!$          if( pw(i) == strlen) exit
!!$          if( pw(i) > strlen) call errors('occurrences_out_of_string')
!!$          if( pw(i+1) == 0) pw(i+1) = strlen+1                ! The last word should be entire
!!$          stringout(i+1-start) = stringin(pw(i)+1:pw(i+1)-1)
!!$       end if
!!$    end do

    start=1; nword=0                                                 
    if (pw(1) > 1) start = 0                                ! If the firt position is not an occurrence..
    do i = start,occn                                       ! this two lines are needed                  
       if( pw(i) /= pw(i+1)-1 .and. pw(i) < strlen) then    ! If there is two occurrences in a row, take only
                                                            ! the last one and if the last character of string
                                                            ! is an occurrence, don'mind about it
          if( i == occn ) pw(i+1) = strlen+1                ! If the last occurrence is not in the last position
                                                            ! of string, add a position in pw vector, in this way
                                                            ! we are be able to read the last word also
!          if( nword < 1 .and. pw(i) > 1 .or. pw(i) /= pw(i-1)+1 ) pw(i-1) = 1
          nword = nword + 1                                 ! Counter for the word
          stringout(nword) = stringin(pw(i)+1:pw(i+1)-1)    ! Shrinking between pw(i)+1 and pw(i+1)-1
       end if
    end do
    
!!$    do i = 1,nword                                          
!!$       stringout(i) = trim(adjustl(stringout(i)))
!!$    end do
    !print*, '#',stringout(1),'#' ! debug
    !print*, nword ! debug
!    stop   ! debug
  END SUBROUTINE shrink_string

  
  SUBROUTINE errors(type)   ! This routine manages the errors in the module
    IMPLICIT NONE
    
    character(*),intent(IN) :: type
    
    select case (type)
    case ('exceeds_words')
       write(0,9001) module_name,maxwords
       stop
    case ('exceeds_chars')
       write(0,9002) module_name,maxchars
       stop
    case ('occurrences_out_of_string')
       write(0,9003) module_name,strlen
       stop
    case default
       write(0,9000) module_name
       stop
    end select
    
    ! Error Expression
9000 format('ERROR: UNKNOWN in module ',A15)
9001 format('ERROR: Module: ',A15,' word number exceeds memory ',/,' [maxwords = ',I6,']')
9002 format('ERROR: Module: ',A15,' word number exceeds memory ',/,' [maxchars = ',I6,']')
9003 format('ERROR: Module: ',A15,' found occurrences out of string ',/,' [strlen = ',I6,']')
    
  END SUBROUTINE errors
  
END MODULE string_tools
