PROGRAM SCC_Disp
  
  USE precision
  USE utils
  USE read_tag_dftbp
  USE read_xyz
  USE read_atomdata
  USE ReadInput  ! TODO better
  USE parameter, only : BohrAngst,HartKcalMol
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
  
  real(kr) :: Rab,ba,bb,bab,damp

  integer(ki) :: i,j,k

!  integer(ki) :: b0                                         --> Red from ReadInput
!  character(30) :: inputtagfile,inputcoofile,atomdatafile   --> Red from ReadInput
  
!  b0 = 1
!  inputtagfile = 'results.tag'
!  inputcoofile = '91.xyz'
!  atomdatafile = 'atomdata.data'

  call read_stdin
!  stop ! debug
  ! Mesure Units:
  ! Ni -> Electron fraction
  ! Coordinates -> Angstrom
  ! Polarizability -> Angstrom^3
  ! C6 coefficients -> Atomic Units
  !
  ! So I think the most simple thing is 
  ! transform all dimension in atomic
  ! units and, maybe, convert energy from 
  ! Hartree/molecule in Kcal/mol

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
  !    E_disp = sum_A sum_B>A F_damp (C6_AB / R_AB)
  !
  ! Which F_damp??
  ! In TS's paper, is used a Fermi-type damping function:
  !
  !  f_damp = (1 + exp[ -d(R_ab/s*R0_ab -1)])^-1
  !
  ! They do a strange discussion about the way to choose R0_ab and R_ab
  ! saying that the VdW by Bondi are radii for AIM and so we should not
  ! use this as Radii for free atoms. 
  ! Speaking with Stephan, "I" decided that it could be simpler start with 
  ! simpler damping function: I will use as first Fdamp the one present in 
  ! J Chem Teory Comp 7 3567 and called f_2n little modified: b(x) is defined as
  ! F(x)*b_ij_asym while I assume F(x) = 1. At the end, I have:
  ! 
  ! f_dump_2n = 1 - exp(-x) * sum_k=0^2n x^k / k!
  ! 
  ! where
  !
  ! x = b_ij,asym * R_ij
  !
  ! where
  ! 
  ! b_ij,asym = 2 * b_ii,asym * b_jj,asym /( b_ii,asym + b_jj,asym )
  !
  ! where
  ! 
  ! b_ii,asym = b_0 * ( 1 / alpha_i,free)^(1/3) * ( V_i,free / V_i,AIM)^(1/3)
  !
  ! the last term, in this case, will be substituited by the "weight" calculated
  ! in this routine:
  !
  ! V_i,free / V_i,AIM = Z_i / N_i (see above)

  ! If I would reduce memory consumption, I can compute the C6AIM 
  ! on-the-fly while I compute energy... In this first version I don't 
  ! want to worry about memory consumption... also because I have only 
  ! very little molecules.

  

  !!!!!!!!!!!!!!! COMPUTE WEIGHTS & C6AIM !!!!!!!!!!!!!!!!

  ! Retrive data from results.tag
  ! + When the spin polarization is used in the calculation with 
  ! + dftb+, the atomic charges are spanned on a matrix of dimension
  ! + natom x 2. Doing some test about the value in the second column, 
  ! + I realized that only the first column is important: I compared the
  ! + value of "Net atomic charges" in the detailed.out file with the
  ! + values in the matrix. I have seen that only the values in the first
  ! + column are coerent with the value in detailed.out. I'm sure that the 
  ! + called "Net atomic charges" are independent from spin polarization because
  ! + I obtained the same value summing the electronic popolutaion of each spin.
  ! + These values of electronic population are present below, still in the 
  ! + detailed.out file.
  call get_tag_data('atomic_charges',inputtagfile) 
  Ni => matrix(:,1,1)    ! WARNING:: If you recall 'get_tag_data' tag matrix change its values
  Ni_size => ifrmt(1)
  
  ! Retrieve data form atomdata.data
  call get_atomdata(atomdatafile) ! Now I can use the atomdata array
  
  ! Retrive name of atom
  call get_coords(inputcoofile,natom)
  
  if( Ni_size /= natom ) stop 'ERROR: check 1' ! Controls
  
  ! Compute C6AIM
  allocate(C6aim(natom,2),C6free(natom,2),Zaim(natom),polar(natom))

  do i = 1,natom                                                  ! Some useful arrays:
     Zaim(i) = atomdata(findatom(coords(i)%atom_type))%Z          ! - Z_i of atoms in xyz file
     C6free(i,1) = atomdata(findatom(coords(i)%atom_type))%TSC6   ! - C6free by TS for "   "
     C6free(i,2) = atomdata(findatom(coords(i)%atom_type))%D3C6   ! - C6free by 3D for "   "
     polar(i) = atomdata(findatom(coords(i)%atom_type))%polarizability / BohrAngst**3 ! -polar for    "   "
  end do
 

!  write(*,'(I5 /)') (Zaim(i), i = 1,natom) ! debug
!  write(*,'(f8.3 /)') (C6free(i,1), i = 1,natom) ! debug
!  write(*,'(f8.3 /)') (C6free(i,2), i = 1,natom) ! debug

  
  ! => C6AIM
  do i = 1,2                                                     ! Compute the C6aim form
     C6aim(:,i) = (Ni(:) / dble(Zaim(:)))**2 * C6free(:,i)                  ! both: TS and 3D
  end do                                                         ! 1 = TS ; 2 = 3D

!  write(*,'(f8.3 /)') (C6eff(i,1), i = 1,natom) ! debug
!  write(*,'(f8.3 /)') (C6eff(i,2), i = 1,natom) ! debug
  



  !!!!!!!!!!!!!!! COMPUTE C6BIAT & ENERGY !!!!!!!!!!!!!!!!

  E = 0.0
  atom1: do i = 1,natom
     atom2: do j = i+1,natom
        Rab = dist(coords(i)%coord,coords(j)%coord)/BohrAngst
        ba = basym(b0,polar(i),Ni(i),Zaim(i))
        bb = basym(b0,polar(j),Ni(j),Zaim(j))
        bab = bmix(ba,bb)
        damp = fdamp(bab,Rab)
!        print*, i,j,Rab,bb,ba,damp ! debug
!        print*, Rab, i ,j,coords(i)%coord,coords(j)%coord ! debug
!        stop ! debug 
        TSvsD3loop : do k = 1,2  ! 1 => TS ; 2 => D3
           C6ab(1) = 2 * C6aim(i,k) * C6aim(j,k) / &
                ( C6aim(i,k) + C6aim(j,k) )                  ! REL A
           C6ab(2) = 2 * C6aim(i,k) * C6aim(j,k) / &
                ( ( polar(j) / polar(i) ) * C6aim(i,k) + &
                ( polar(i) / polar(j) ) * C6aim(j,k) )       ! REL B
           E(:,k) =  E(:,k) - damp * C6ab(:) / Rab**6                ! Not implemented damping function yet. 
        end do TSvsD3loop
     end do atom2
  end do atom1

  
  ! E(1,:) => Energy of relationa A
  ! E(2,:) => Energy of relationa B
  ! E(:,1) => Energy from TS
  ! E(:,2) => Energy from D3

  write(*,*) '#  ',inputcoofile
  write(*,'("A-TS  ",f10.3)') E(1,1) * HartKcalMol
  write(*,'("A-D3  ",f10.3)') E(1,2) * HartKcalMol
  write(*,'("B-TS  ",f10.3)') E(2,1) * HartKcalMol
  write(*,'("B-D3  ",f10.3)') E(2,2) * HartKcalMol

!!$  do i = 1,100     ! debug
!!$     !     write(0,*) bab,bb,ba      ! debug
!!$     kkk =  fdamp(1d0,dble(i)/10)     ! debug
!!$     print*, kkk     ! debug
!!$     !     write(*,*) fdamp(dble(1),dble(i/10))     ! debug
!!$  end do


  
!  do i = 1,natom/2       ! debug
!     r12 = dist(mol1(i)%coord,mol2(j)%coord)      ! debug
!  end do      ! debug


!  call getdata('atomic_charges','results.tag')      ! debug
! write(*,*) (((matrix(i,j,k), i=1,ifrmt(1)),j=1,ifrmt(2)),k=1,ifrmt(3))      ! debug
!  write(*,*) coords      ! debug
CONTAINS

  real(kr) FUNCTION bmix(bi,bj)
    IMPLICIT NONE
    real(kr),intent(IN) :: bi,bj

    bmix = 2 * bi * bj / (bi + bj)
    
  END FUNCTION bmix

  real(kr) FUNCTION basym(b0,alpha,N,Z)
    IMPLICIT NONE
    real(kr),intent(IN) :: b0,alpha,N
    integer(ki),intent(IN) :: Z
    
    basym = b0 * ( 1 / alpha )**(1./3.) * ( real(Z) / N )**(1./3.)
    
!    print*, 'basym', b0,alpha,Z,N  ! debug
    
  END FUNCTION basym
  
  real(kr) FUNCTION fdamp(b,R)
    IMPLICIT NONE
    real(kr),intent(IN) :: b,R
!    real(kr) :: asd
    integer(ki) :: i

    fdamp = 0.0
!    print*, '-----------------> ',R ! debug
    do i = 0,6
       fdamp = fdamp + (b*R)**i / fact(i)
!       print*, i,asd,fact(i),(b*R)**i / fact(i)  ! debug
    end do
    fdamp = 1 - exp( -b * R) * fdamp
!    print*, 'fdamp ',fdamp ! debug
    
  END FUNCTION fdamp
  
END PROGRAM SCC_Disp
  
