MODULE atom_data
  USE precision
  USE file_tools
  
  IMPLICIT NONE
  
  type,public :: atomdata_type
     character(5) :: atom_type
     integer(ki)  :: Z 
     real(kr)     :: polarizability
     real(kr)     :: TSC6
     real(kr)     :: D3C6
  end type atomdata_type
  type(atomdata_type),allocatable,dimension(:) :: atomdata
  

  PRIVATE
  PUBLIC :: atomdata, get_atomdata

CONTAINS
  
  SUBROUTINE get_atomdata(file)
    character(*),intent(IN) :: file
    character(30) :: junk
    integer(ki) :: natom
    integer(ki) :: i

    call openfile(file,'read')
    
    natom = 0
    do
       read(fiit(file),*,iostat=err) junk
       if( err < 0 ) exit
       if( err > 0 ) stop 'ERROR reading file (get_atomdata 1)'
       natom = natom + 1
    end do
    rewind(fiit(file))
    
    allocate(atomdata(natom))

    do i = 1,natom
       read(fiit(file),*,iostat=err) atomdata(i)%atom_type, &
            atomdata(i)%Z,atomdata(i)%polarizability, &
            atomdata(i)%TSC6,atomdata(i)%D3C6
       if( err /= 0 ) stop 'ERROR reading file (get_atomdata 2)'
    end do
    
  END SUBROUTINE get_atomdata
  

END MODULE atom_data
