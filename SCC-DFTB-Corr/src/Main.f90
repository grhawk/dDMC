MODULE read_result_tag_dftbp
  USE precision
  USE string_tools
  IMPLICIT NONE
  character*30 :: junk!,junk1,junk2
  real :: x
  
CONTAINS
  
  SUBROUTINE read_property(iu,name) !input_unit,property_name
    IMPLICIT NONE
    character*19, intent(in) :: name !name of property (e.g. fillings)
    integer(ki),intent(in)   :: iu   !name of reading unit
    
    character(30) :: form, form_shrk(3)
!    character(15),allocatable :: 
    integer(ki)  :: dims
    
    
    doread:do
       read(iu,*) junk, form
       if ( junk .eq. name) then

          ! Determinare il formato          
          call shrink_count(form,':',dims)
          ! form_shrk(1) => variable type
          ! form_shrk(2) => rank
          ! form_shrk(3) => format
          select case (form_shrk(1))
          case ('real')
             select case (form_shrk(2))
             case ('0')
                print*, 'scalar' ! not implemented
             case ('1')
                print*, 'vector' ! not implemented
             case('2')
                print*, 'matrix rank 2' ! not implemented
             case('3')
                print*, 'matrix rank 2' ! not implemented
             case default
                print*, 'Hai rotto il cazzo' ! not implemented
             end select
          case ('logical')
             print*, 'logical' ! not implemented
          case ('integer')
             print*, 'integer' ! not implemented
          case default
             print*, 'Hai rotto il cazzo' ! not implemented
          end select

             
             
          

          ! Creare un array allocandolo in base al formato
          ! Leggere tutti i numeri e metterli nell'array
          exit doread
       end if
    end do doread

    print*, junk,' -- ',form

  END SUBROUTINE read_property
  
  SUBROUTINE read_write
    
    read(*,*) x
    write(*,*) x
      
  END SUBROUTINE read_write
    
END MODULE read_result_tag_dftbp
  
  
PROGRAM SCC_Disp
  
  USE read_result_tag_dftbp
  
  !  CALL read_write
  
  OPEN(501,file='results.tag',form='FORMATTED',status='OLD')
  
  call read_property(501,junk2)
    
END PROGRAM SCC_Disp
  
