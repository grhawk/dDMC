MODULE read_xyz

! Calling get_data(file,dummy) from main (where file !
! is the name of the xyz file) you obtain the number !
! of atom in the file as "dummy" and will be         !
! generated a variable, namely coords, that contains !
! the type of atom and the coordinates of the atom.  !
! This variable is of derived type:                  !
! - coords(i)%atom_type => atom type of i-esim atom  !
! - coords(i)%coord(j) => j-esim coordinate of i-esim!
! atom (of course I have three coordinate for each   !
! atom and, of course, these are disposed as         !
! 1 => x                                             !
! 2 => y                                             !
! 3 => z                                             !
! )                                                  !
!
! 26.04.2012
! Added a function in order to know if an atom is or !
! not an Hydrogen from the atom index in the molecule!
! This function is important to apply the            !
! HH-repulsion correction to dDMC correction.        !

USE precision
USE file_tools

IMPLICIT NONE

TYPE,public :: xyz_coords
   integer(ki)   :: index
   character(5)  :: atom_type
   real(kr),dimension(3) :: coord
END type xyz_coords

character(64) :: commentline

type(xyz_coords),target,allocatable :: coords(:)

PRIVATE
PUBLIC :: get_coords, coords, IsHAtom

CONTAINS

  SUBROUTINE get_coords(file,natom)
    IMPLICIT NONE  
    character(*),intent(IN) :: file
    integer(ki),intent(OUT) :: natom
    character(kch) :: junk

    integer(ki) :: err
    integer(ki) :: i,j,j1,j2
  
    call openfile(file,'read')
    
    read(fiit(file),*) natom
    read(fiit(file),'(A64)') commentline  ! The format prevent reading the new line
!    read(fiit(file),*) junk
!    if( junk(1:1) == 'H' .or. junk(1:1) == 'C' .or.  &
!         junk(1:1) == 'O' .or. junk(1:1) == 'N' .or. &
!         junk(1:1) == 'S' .or. junk(1:1) == 'P' ) backspace(fiit(file))

    allocate(coords(natom)) ! Needs more options
    
!    print*, natom  ! debug

    do i = 1,natom
       coords(i)%index = i
       read(fiit(file),*,iostat=err) coords(i)%atom_type, &
            coords(i)%coord(1),coords(i)%coord(2),coords(i)%coord(3)
!       print*, coords(i), i,natom ! debug
       if( err /=0) exit
    end do

    call closefile(file)

  END SUBROUTINE get_coords

  logical FUNCTION IsHAtom(index)
    IMPLICIT NONE
    integer(ki),intent(IN) :: index
    character(1) :: HAtom(2)
    integer(ki) :: l
    

    HAtom=(/ ' H',' h' /)

    do l = 1,2
       if( trim(adjustl(coords(index)%atom_type)) == trim(adjustl(HAtom(l))) ) IsHAtom = .true.
    end do
    
    
  END FUNCTION IsHAtom
  


END MODULE read_xyz

!!$program test
!!$  USE read_xyz
!!$  USE precision
!!$
!!$  IMPLICIT NONE
!!$  character(30) :: file='91.xyz'
!!$  integer(ki) :: natom,i
!!$
!!$  call get_coords(file,natom)
!!$  
!!$  do i = 1,natom
!!$     write(*,*) coords(i)!coords(i)%atom_type, coords(i)%coord(:)
!!$  end do
!!$
!!$end program test
