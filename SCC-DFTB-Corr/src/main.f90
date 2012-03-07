PROGRAM SCC_Disp
  
  USE precision
  USE read_tag_dftbp
  USE read_xyz
  USE read_atomdata
!  USE parser_input  ! Not implemented - not indispensable
!  USE parameter  ! Not implemented - not indispensable
  IMPLICIT NONE
  integer(ki) :: natom

  type(xyz_coords),pointer :: mol1(:),mol2(:)

  integer(ki),pointer  :: Ni_size   ! each time that get_tag_data is called!!
  real(kr),dimension(:),pointer :: Ni

  integer(ki),dimension(:),allocatable:: Zaim
  real(kr),dimension(:,:),allocatable :: C6aim,C6free
  real(kr),dimension(:),allocatable :: polar


  real(kr),dimension(2,2) :: E
  real(kr),dimension(2)  :: C6ab

  real(kr) :: Rab

  integer(ki) :: i,j,k

  character(30) :: inputtagfile,inputcoofile,atomdatafile
  
  inputtagfile = 'results.tag'
  inputcoofile = '91.xyz'
  atomdatafile = 'atomdata.data'
  
  ! Compute weights
  !
  !     W_i = N_i / Z_i
  !
  !
  ! Compute C6AIM
  !
  !     C6AIM = W_i * C6FREE
  !
  !
  ! Compute C6BIAT
  !
  !     C6_AB = 2 * (C6_A * C6_B) / (C6_A + C6_B)
  !
  !
  ! Compute Energy
  ! 
  !    E_disp = sum_A,B (C6_AB / R_AB)

  ! If I would reduce memory consumption, I can compute the C6AIM 
  ! on-the-fly while I compute energy... In this first version I don't 
  ! want to worry about memory consumptio... also because I have only 
  ! very little molecules.

  

  !!!!!!!!!!!!!!! COMPUTE WEIGHTS & C6AIM !!!!!!!!!!!!!!!!

  ! Retrive data from results.tag
  call get_tag_data('atomic_charges',inputtagfile)
  Ni => matrix(:,1,1)    ! WARNING:: If you recall 'get_tag_data' tag matrix change its values
  Ni_size => ifrmt(1)
  
  ! Retrieve data form atomdata.data
  call get_atomdata(atomdatafile) ! Now I can use the atomdata array
  
  ! Retrive name of atom
  !++ To understand what are the two molecules, I separe the red coords in two 
  !++ set. Since I should have the same molecule twice into the file xyz, I assume
  !++ that the first half of coords are associated at one molecule and the second
  !++ are associated at the second.
  call get_coords(inputcoofile,natom)
  if( mod(natom,2) /= 0 ) stop 'ERROR: check 1' ! controls
  mol1 => coords(1:natom/2)      ! WARNING:: If you recall 'get_coords' these values change
  mol2 => coords(natom/2:natom)  !
  

  if( Ni_size /= natom ) stop 'ERROR: check 1' ! Controls
  
  ! Compute C6AIM
  allocate(C6aim(natom,2),C6free(natom,2),Zaim(natom),polar(natom))

  do i = 1,natom                                                  ! Some useful arrays:
     Zaim(i) = atomdata(findatom(coords(i)%atom_type))%Z          ! - Z_i of atoms in xyz file
     C6free(i,1) = atomdata(findatom(coords(i)%atom_type))%TSC6   ! - C6free by TS for "   "
     C6free(i,2) = atomdata(findatom(coords(i)%atom_type))%D3C6   ! - C6free by 3D for "   "
     polar(i) = atomdata(findatom(coords(i)%atom_type))%polarizability ! -polar for    "   "
  end do
 
!  write(*,'(I5 /)') (Zaim(i), i = 1,natom) ! debug
!  write(*,'(f8.3 /)') (C6free(i,1), i = 1,natom) ! debug
!  write(*,'(f8.3 /)') (C6free(i,2), i = 1,natom) ! debug

  
  ! => C6AIM
  do i = 1,2                                                     ! Compute the C6aim form
     C6aim(:,i) = Ni(:) / Zaim(:) * C6free(:,i)                  ! both: TS and 3D
  end do                                                         ! 1 = TS ; 2 = 3D

!  write(*,'(f8.3 /)') (C6eff(i,1), i = 1,natom) ! debug
!  write(*,'(f8.3 /)') (C6eff(i,2), i = 1,natom) ! debug
  



  !!!!!!!!!!!!!!! COMPUTE C6BIAT & ENERGY !!!!!!!!!!!!!!!!

  E = 0.d0
  mol1loop : do i = 1,natom/2          ! loop on first molecule atoms
     mol2loop : do j = natom/2,natom   ! loop on first molecule atoms
        Rab = dist(coords(i)%coord,coords(j)%coord)
!        print*, Rab, i ,j,coords(i)%coord,coords(j)%coord ! debug
!        stop
        TSvsD3loop : do k = 1,2  ! 1 => TS ; 2 => D3
           C6ab(1) = 2 * C6aim(i,k) * C6aim(j,k) / &
                ( C6aim(i,k) + C6aim(j,k) )                  ! REL A
           C6ab(2) = 2 * C6aim(i,k) * C6aim(j,k) / &
                ( ( polar(j) / polar(i) ) * C6aim(i,k) + &
                  ( polar(i) / polar(j) ) * C6aim(j,k) )       ! REL B
           E(:,k) = E(:,k) + C6ab(:) / Rab**6
        end do TSvsD3loop
     end do mol2loop
  end do mol1loop
  
  ! E(1,:) => Energy of relationa A
  ! E(2,:) => Energy of relationa B
  ! E(:,1) => Energy from TS
  ! E(:,2) => Energy from D3

  write(*,*) '#  ',inputcoofile
  write(*,'("A-TS  ",f10.3)') E(1,1)
  write(*,'("A-D3  ",f10.3)') E(1,2)
  write(*,'("B-TS  ",f10.3)') E(2,1)
  write(*,'("B-D3  ",f10.3)') E(2,2)


!  do i = 1,natom/2
!     r12 = dist(mol1(i)%coord,mol2(j)%coord)
!  end do


!  call getdata('atomic_charges','results.tag')
! write(*,*) (((matrix(i,j,k), i=1,ifrmt(1)),j=1,ifrmt(2)),k=1,ifrmt(3))
!  write(*,*) coords
CONTAINS

  real(kr) FUNCTION dist(r1,r2)
!    USE precision
    IMPLICIT NONE
    real(kr),intent(IN),dimension(3) :: r1,r2
    real(kr),dimension(3,2) :: r12
    integer(ki) :: i
    
    r12(:,1) = r1(:)
    r12(:,2) = r2(:)

    do i = 1,3
       dist = dist + ( r12(i,1) - r12(i,2) ) **2
!       print*, dist  ! debug
    end do

    dist = sqrt(dist)
    
  END FUNCTION dist
  
    
!!$  FUNCTION get_C6eff(name,Ni)
!!$    USE atom_data
!!$    IMPLICIT NONE
!!$    real(kr),dimension(2) :: C6eff
!!$    real(kr),intent(IN) :: Ni
!!$    character(*),intent(IN) :: name
!!$    
!!$    
!!$  END FUNCTION Get_C6eff

END PROGRAM SCC_Disp
  
