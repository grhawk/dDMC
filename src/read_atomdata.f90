MODULE read_atomdata
  USE precision
  USE file_tools
  
  IMPLICIT NONE
  
  type,public :: atomdata_type
     integer(ki)  :: Z 
     character(5) :: atom_type
     real(kr)     :: polarizability
     real(kr)     :: TSC6
     real(kr)     :: D3C6
  end type atomdata_type
  type(atomdata_type),allocatable,target :: atomdata(:)
  

 INTERFACE findatom
    MODULE PROCEDURE atom_by_name, atom_by_number
  END INTERFACE

  PRIVATE
  PUBLIC :: atomdata, get_atomdata, findatom
  
CONTAINS
  
  SUBROUTINE get_atomdata(file)
    IMPLICIT NONE
    character(*),intent(IN) :: file
    character(1) :: junk
    integer(ki) :: natom,err,nl
    integer(ki) :: i

    call openfile(file,'read')
    
    natom = 0; nl = 0
    do
       read(fiit(file),*,iostat=err) junk
       if( err < 0 ) exit
       if( err > 0 ) stop 'ERROR reading file (get_atomdata 1)'
       if( junk(1:1) /= '#' ) natom = natom + 1
       nl = nl + 1
    end do
    rewind(fiit(file))
    
!    print*, 'natom ',natom,'    nl',nl  ! debug
    
    allocate(atomdata(natom))

    natom = 0
    do i = 1,nl
       read(fiit(file),*) junk
!       print*, junk ! debug
       if( junk /= '#' ) then 
          natom = natom + 1
          backspace(fiit(file))
          read(fiit(file),*,iostat=err) atomdata(natom)%Z, &
               atomdata(natom)%atom_type, atomdata(natom)%polarizability, &
               atomdata(natom)%TSC6, atomdata(natom)%D3C6
!          print*, atomdata(i)%Z  ! debug
          if( err /= 0 ) stop 'ERROR reading file (get_atomdata 2)'
!          print*, 'asd ', i   ! debug
       end if
    end do
    
  END SUBROUTINE get_atomdata

  integer(ki) FUNCTION atom_by_name(atom)
    IMPLICIT NONE
    character(*),intent(IN) :: atom
!    integer(ki),intent(OUT) :: atom_by_name
    integer(ki) :: i

    do i = 1,size(atomdata)
       if( atomdata(i)%atom_type == atom ) then
          atom_by_name = i
          return
       end if
    end do
    
    atom_by_name = 0 ! If atom is not present in the array give 0 as value
    
  END FUNCTION atom_by_name

  integer(ki) FUNCTION atom_by_number(number)
    IMPLICIT NONE
    integer(ki),intent(IN) :: number
    integer(ki) :: i

    do i = 1,size(atomdata)
       if( atomdata(i)%Z == number ) then
          atom_by_number = i
          return
       end if
    end do
    
    atom_by_number = 0 ! If atom is not present in the array give 0 as value
    
  END FUNCTION atom_by_number
 
    


END MODULE read_atomdata


!!$program test
!!$  USE atom_data
!!$
!!$  IMPLICIT NONE
!!$  integer :: i
!!$  integer :: natom = 54
!!$ 
!!$  call get_atomdata('atomdata.data')
!!$  write(*,'(I5,A10,f18.6,2f15.4)') atomdata(:)
!!$  
!!$  print*, findatom(8)
!!$  print*, atomdata(findatom('C'))%TSC6
!!$  
!!$end program test
