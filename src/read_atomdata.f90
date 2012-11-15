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
     real(kr)     :: vdWr
     real(kr)     :: incharge    !innershell electron number
  end type atomdata_type
  type(atomdata_type),allocatable,target :: atomdata(:)

! If the next is .true. the program need a file to read atomic
! properties, if false, the program will use the atomic properties
! reported in this subroutine.
  logical,parameter :: read_from_file = .false.

 INTERFACE findatom
    MODULE PROCEDURE atom_by_name, atom_by_number
  END INTERFACE

  PRIVATE
  PUBLIC :: atomdata, get_atomdata, findatom
  
CONTAINS
  
  SUBROUTINE get_atomdata(file)
    IMPLICIT NONE
    character(*),intent(IN) :: file
    character(kch) :: junk
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
       
       call closefile(file)

    else

       allocate(atomdata(54))

       atomdata%polarizability = 0.0d+0
       atomdata%TSC6 = 0.0d+0
       atomdata%D3C6 = 0.0d+0

       atomdata%Z = (/ 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20&
            &,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40&
            &,41,42,43,44,45,46,47,48,49,50,51,52,53,54 /)

       do i = 1,54
          if( i <= 2 )                atomdata(i)%incharge = 0
          if( i >  2 .and. i <= 10 )  atomdata(i)%incharge = 2
          if( i > 10 .and. i <= 18 )  atomdata(i)%incharge = 10
          if( i > 18 .and. i <= 36 )  atomdata(i)%incharge = 18
          if( i > 26 .and. i <= 54 )  atomdata(i)%incharge = 36
       end do

       atomdata%atom_type = (/ ' H', 'He', 'Li', 'Be', ' B', ' C', ' N', '&
            & O', ' F', 'Ne', 'Na', 'Mg', 'Al', 'Si', ' P', ' S', 'Cl', 'A&
            &r', ' K', 'Ca', 'Sc', 'Ti', ' V', 'Cr', 'Mn', 'Fe', 'Co', '&
            &Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',&
            & 'Rb', 'Sr', ' Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd&
            &', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', ' I', 'Xe' /)
       
       ! Provide from CRC 10-186
       ! Units = 10^{-24} cm^{3}
       ! Checked on 15.11.2012
       atomdata%polarizability = (/ 0.666793d+0, 0.2050522d+0, 24.33d+0&
            &, 5.60d+0, 3.03d+0, 1.67d+0, 1.10d+0, 0.802d+0, 0.557d+0&
            &, 0.3956d+0, 24.11d+0, 10.6d+0, 6.8d+0, 5.38d+0, 3.63d+0&
            &, 2.90d+0, 2.18d+0, 1.6411d+0, 43.06d+0, 22.8d+0, 17.8d+0,&
            & 14.6d+0, 12.4d+0, 11.6d+0, 9.4d+0, 8.4d+0, 7.5d+0, 6.8d&
            &+0, 6.2d+0, 5.75d+0, 8.12d+0, 6.07d+0, 4.31d+0, 3.77d+0,&
            & 3.05d+0, 2.484d+0, 47.3d+0, 27.6d+0, 22.7d+0, 17.9d+0,&
            & 15.7d+0, 12.8d+0, 11.4d+0, 9.6d+0, 8.6d+0, 4.8d+0, 7.2d&
            &+0, 7.36d+0, 10.2d+0, 7.7d+0, 6.6d+0, 5.5d+0, 5.35d+0,&
            & 4.044d+0 /) 

!!$       ! Provide from CRC 10-186
!!$       ! Units = 10^{-24} cm^{3}
!!$       atomdata%polarizability = (/ 0.666793d+0, 0.204956d+0, 24.3d+0&
!!$            &, 5.60d+0, 3.03d+0, 1.67d+0, 1.10d+0, 0.802d+0, 0.557d+0&
!!$            &, 0.3956d+0, 24.11d+0, 10.6d+0, 6.8d+0, 5.38d+0, 3.63d+0&
!!$            &, 2.90d+0, 2.18d+0, 1.641d+0, 43.4d+0, 22.8d+0, 17.8d+0,&
!!$            & 14.6d+0, 12.4d+0, 11.6d+0, 9.4d+0, 8.4d+0, 7.5d+0, 6.8d&
!!$            &+0, 6.2d+0, 5.75d+0, 8.12d+0, 6.07d+0, 4.31d+0, 3.77d+0,&
!!$            & 3.05d+0, 2.484d+0, 47.3d+0, 27.6d+0, 22.7d+0, 17.9d+0,&
!!$            & 15.7d+0, 12.8d+0, 11.4d+0, 9.6d+0, 8.6d+0, 4.8d+0, 7.2d&
!!$            &+0, 7.36d+0, 10.2d+0, 7.7d+0, 6.6d+0, 5.5d+0, 5.35d+0,&
!!$            & 4.044d+0 /) 

       
       ! All the C6 values are in atomic units
       atomdata%TSC6 = (/ 6.498d+0, 1.42d+0, 1392d+0, 227d+0, 99.5d+0&
            &, 46.6d+0, 24.2d+0, 15.6d+0, 9.52d+0, 6.38d+0, 1518d+0,&
            & 626d+0, 528d+0, 305d+0, 185d+0, 134d+0, 94.6d+0, 64.2d&
            &+0, 3923d+0, 2163d+0, 1383d+0, 1044d+0, 832d+0, 602d+0,&
            & 552d+0, 482d+0, 408d+0, 373d+0, 253d+0, 284d+0, 498d+0,&
            & 354d+0, 246d+0, 210d+0, 162d+0, 130d+0, 4691d+0, 3175d&
            &+0, 38.19d+0, 0.0d+0, 0.0d+0, 0.0d+0, 0.0d+0, 0.0d+0,&
            & 0.0d+0, 0.0d+0, 0.0d+0, 0.0d+0, 779d+0, 659d+0, 492d+0,&
            & 445d+0, 385d+0, 286d+0 /)

       atomdata%D3C6 = (/ 7.590d+0, 1.560d+0, 1163.450d+0, 257.490d+0&
            &, 107.180d+0, 49.110d+0, 25.270d+0, 15.510d+0, 9.690d+0,&
            & 6.290d+0, 1608.030d+0, 683.380d+0, 540.540d+0, 317.860d&
            &+0, 191.690d+0, 134.010d+0, 92.350d+0, 64.650d+0,&
            & 4983.500d+0, 2352.690d+0, 1522.470d+0, 1361.920d+0,&
            & 1116.100d+0, 690.740d+0, 802.750d+0, 491.330d+0,&
            & 532.780d+0, 574.740d+0, 337.180d+0, 340.520d+0,&
            & 483.750d+0, 363.550d+0, 262.950d+0, 213.670d+0,&
            & 167.130d+0, 130.400d+0, 6138.780d+0, 3381.370d+0,&
            & 2365.890d+0, 1822.720d+0, 1475.250d+0, 845.900d+0,&
            & 1067.020d+0, 598.200d+0, 713.940d+0, 608.500d+0,&
            & 426.750d+0, 468.190d+0, 757.740d+0, 627.570d+0,&
            & 492.940d+0, 425.540d+0, 351.970d+0, 290.220d+0 /) 

       !van der Waals radii from CRC
       !Angstrom (10^{-8}cm)
!!$       atomdata%vdWr = (/ 1.10d0, 1.40d0, 10000, 10000, 10000, 1.70d0, 1.55d0, &
!!$            & 1.52d0, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 1.80d0, 10000,&
!!$            &10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, &
!!$            &10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000,&
!!$            & 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000,&
!!$            &10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000 /)

       ! Modified
       atomdata%vdWr = (/ 1.10d0, 1.40d0, 10000d0, 10000d0, 10000d0, 1.70d0, 1.55d0, &
            & 1.52d0, 10000d0, 10000d0, 10000d0, 10000d0, 10000d0, 10000d0, 1.80d0, 1.80d0, 1.8d0,&
            &10000d0, 10000d0, 10000d0, 10000d0, 10000d0, 10000d0, 10000d0, 10000d0, 10000d0, 10000d0, &
            &10000d0, 10000d0, 10000d0, 10000d0, 10000d0, 10000d0, 10000d0, 10000d0, 10000d0,&
            & 10000d0, 10000d0, 10000d0, 10000d0, 10000d0, 10000d0, 10000d0, 10000d0, 10000d0, 10000d0,&
            &10000d0, 10000d0, 10000d0, 10000d0, 10000d0, 10000d0, 10000d0, 10000d0 /)


    end if
    
  END SUBROUTINE get_atomdata

  integer(ki) FUNCTION atom_by_name(atom)
    IMPLICIT NONE
    character(*),intent(IN) :: atom
!    integer(ki),intent(OUT) :: atom_by_name
    integer(ki) :: i

    do i = 1,size(atomdata)
!       print*, '#',trim(atom),'#',trim(adjustl(atomdata(i)%atom_type)),'#'
       if( trim(adjustl(atomdata(i)%atom_type)) == trim(adjustl(atom)) ) then
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
