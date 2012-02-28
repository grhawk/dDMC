PROGRAM SCC_Disp
  
  USE precision
  USE read_tag_dftbp, getdata => get_dtfbp_data
  USE read_xyz
!  USE parser_input  ! Not implemented
  

  integer(ki) :: natom
  type(xyz_coords),pointer :: mol1(:),mol2(:)
  
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
  ! Compute C6MOL
  !
  !     C6_AB = 2 * (C6_A * C6_B) / (C6_A + C6_B)
  !
  !
  ! Compute Energy
  ! 
  !    E_disp = sum_A,B (C6_AB / R_AB)


  

  ! Compute Weights

  call get_dftbp_data('atomic_charges','results.tag')



  ! To understand what are the two molecules, I separe the red coords in two 
  ! set. Since I should have the same molecule twice into the file xyz, I assume
  ! that the first half of coords are associated at one molecule and the second
  ! are associated at the second.
  call get_coords('91.xyz',natom)

! write(*,*) (coords(i), i=1,natom)
!  allocate(mol1(natom/2),mol2(natom/2))
  
  mol1 => coords(1:natom/2)
  mol2 => coords(natom/2:natom)

  do i = 1,natom/2
     r12 = dist(mol1(i)%coord,mol2(j)%coord)
  end do


!  call getdata('atomic_charges','results.tag')
! write(*,*) (((matrix(i,j,k), i=1,ifrmt(1)),j=1,ifrmt(2)),k=1,ifrmt(3))
!  write(*,*) coords
CONTAINS

  real(kr) FUNCTION dist(r1,r2)
    IMPLICIT NONE
    real(kr),intent(IN),dimension(3) :: r1,r2
    real(kr),dimension(3,2) :: r12
    integer(ki) :: i
    
    r12(:,1) = r1(:)
    r12(:,2) = r2(:)

    r12 = r12**2
    
    do i = 1,3
       dist = dist + r12(i,1) - r12(i,2)
    end do

    dist = sqrt(dist)
    
  END FUNCTION dist
  
    
END PROGRAM SCC_Disp
  
