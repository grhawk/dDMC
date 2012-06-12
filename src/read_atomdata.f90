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
! Covalent radius
!     real(kr)     :: vdWr
  end type atomdata_type
  type(atomdata_type),allocatable,target :: atomdata(:)

! If the next is .true. the program need a file to read atomic
! properties, if false, the program will use the atomic properties
! reported in this subroutine.
  logical,parameter :: read_from_file=.false. 

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

    if( read_from_file ) then
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
    else
       
       allocate(atomdata(54))

       atomdata%Z = (/ 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20&
            &,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40&
            &,41,42,43,44,45,46,47,48,49,50,51,52,53,54 /)

       atomdata%atom_type = (/ "H", "He", "Li", "Be", "B", "C", "N", "&
            &O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "A&
            &r", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "&
            &Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",&
            & "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd&
            &", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe" /)
       
       atomdata%polarizability = (/ 0.666793, 0.204956, 24.3, 5.60,&
            & 3.03, 1.76, 1.10, 0.802, 0.557, 0.3956, 24.11, 10.6,&
            & 6.8, 5.38, 3.63, 2.90, 2.18, 1.641, 43.4, 22.8, 17.8,&
            & 14.6, 12.4, 11.6, 9.4, 8.4, 7.5, 6.8, 6.2, 5.75, 8.12,&
            & 6.07, 4.31, 3.77, 3.05, 2.484, 47.3, 27.6, 22.7, 17.9,&
            & 15.7, 12.8, 11.4, 9.6, 8.6, 4.8, 7.2, 7.36, 10.2, 7.7,&
            & 6.6, 5.5, 5.35, 4.044 /)
       
       atomdata%TSC6 = (/ 6.498, 1.42, 1392, 227, 99.5, 46.6, 24.2,&
            & 15.6, 9.52, 6.38, 1518, 626, 528, 305, 185, 134, 94.6,&
            & 64.2, 3923, 2163, 1383, 1044, 832, 602, 552, 482, 408,&
            & 373, 253, 284, 498, 354, 246, 210, 162, 130, 4691, 3175&
            &, 38.19, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,&
            & 779, 659, 492, 445, 385, 286 /)

       atomdata%D3C6 = (/ 7.590, 1.560, 1163.450, 257.490, 107.180,&
            & 49.110, 25.270, 15.510, 9.690, 6.290, 1608.030, 683.380&
            &, 540.540, 317.860, 191.690, 134.010, 92.350, 64.650,&
            & 4983.500, 2352.690, 1522.470, 1361.920, 1116.100,&
            & 690.740, 802.750, 491.330, 532.780, 574.740, 337.180,&
            & 340.520, 483.750, 363.550, 262.950, 213.670, 167.130,&
            & 130.400, 6138.780, 3381.370, 2365.890, 1822.720,&
            & 1475.250, 845.900, 1067.020, 598.200, 713.940, 608.500,&
            & 426.750, 468.190, 757.740, 627.570, 492.940, 425.540,&
            & 351.970, 290.220 /)
       

    end if
    
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
