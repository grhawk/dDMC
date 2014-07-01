MODULE read_gf
  USE precision
  USE file_tools
  
  IMPLICIT NONE
  real(kr),allocatable,dimension(:),target :: pop
  integer(ki),target :: npop
  
  PRIVATE
  PUBLIC pop, npop, read_charge_gf
  
CONTAINS
  
  SUBROUTINE read_charge_gf(file)
    IMPLICIT NONE
    character(*),intent(IN) :: file
    
    integer(ki) :: err,i
    character(kch) :: junk
    
    call openfile(file,'read')
    read(fiit(file),*,iostat = err) npop
    allocate(pop(npop))
    
!    print*, 'asd',npop          ! I don't understand why nothing is printed out

    countline: do i = 1,npop
       read(fiit(file),*,iostat = err) pop(i)
       if( err < 0 ) call die('ERROR: some charge in missing')
       if( err > 0 ) call die('ERROR: while reading charge file')
    end do countline
    call closefile(file)
    
!!$    
!!$    call openfile(file,'read')
!!$    print*, fiit(file)
!!$    i = 1
!!$    readfile: do
!!$       read(fiit(file),iostat = err) pop(i)
!!$       print*, err
!!$       if( err < 0 ) exit
!!$       if( err > 0 ) stop 'ERROR: while reading charge file'
!!$       i = i + 1
!!$    end do readfile
!!$    call closefile(file)

!    print*, npop
!    do i =1,npop
!       print*, pop(i)
!    end do
  END SUBROUTINE read_charge_gf
  
END MODULE read_gf
