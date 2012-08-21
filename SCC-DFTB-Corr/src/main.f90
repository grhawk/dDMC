PROGRAM SCC_Disp
  
  USE precision
  USE utils
  USE read_tag_dftbp
  USE read_xyz
  USE read_atomdata
  USE file_tools
  USE ReadInput  ! TODO better
  USE parameters, only : BohrAngst,HartKcalMol
  IMPLICIT NONE
  integer(ki) :: natom

  type(xyz_coords),pointer :: mol1(:),mol2(:)

  integer(ki),pointer  :: Ni_size   ! each time that get_tag_data is called!!
  real(kr),dimension(:),pointer :: Ni

  integer(ki),dimension(:),allocatable:: Zaim
  real(kr),dimension(:,:),allocatable :: C6aim,C6free
  real(kr),dimension(:),allocatable :: polar,rvdw


  real(kr),dimension(2,2) :: E
  real(kr),dimension(2)  :: C6ab
  
  real(kr) :: Rab0,Rab,ba,bb,bab,damp,hhrep

  integer(ki) :: i,j,k

!  integer(ki) :: Hctrl,l
  logical :: Hi,Hj
  logical :: fixb0
  logical :: TTdf
  logical :: Grd
  logical :: GrTTd

  logical :: debug
  character(20),parameter :: chrgfile = 'chargefile.dat'

!  integer(ki) :: b0                                         --> Red from ReadInput
!  character(30) :: inputtagfile,inputcoofile,atomdatafile   --> Red from ReadInput
  
!  b0 = 1
!  inputtagfile = 'results.tag'
!  inputcoofile = '91.xyz'
!  atomdatafile = 'atomdata.data'

  call read_stdin
!  print*, '#',trim(adjustl(DampFunc)),'#'
  select case ( DampFunc )
  case('TT') 
     TTdf = .true.
  case('TTf')
     TTdf = .true.
     fixb0 = .true.
  case('Gr')
     Grd = .true.
  case('GrTT')
     GrTTd = .true.
  case default
     write(0,*) "ERROR: Selected Dumping function not found"
     stop
  end select

  if( debugflag == 'UP' ) debug=.true.
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
  !
  ! The descibed dumping function don't work well for hydrogen interaction. A known
  ! issue of SCC-DFTB is the overbinding of hydrogen system (see Hobza JCTC 8 (2012), 141).
  ! In order to obtain better result, I use a more complex damping function, 
  ! provided from Tang and Tonnies (J. Chem. Phys. 80 (1984), 8) so that at the end
  ! the energy is:
  ! 
  ! E'_disp = A * exp( x ) * E_disp
  !
  ! where E_disp and x are defined before, while A is a parameter that should be
  ! optimized. Indeed this correction will be used to compute energy only between
  ! H atoms. For all other couples I will use the previous equation for energy: E_disp.

  ! If I would reduce memory consumption, I can compute the C6AIM 
  ! on-the-fly while I compute energy... In this first version I don't 
  ! want to worry about memory consumption... also because I have only 
  ! very little molecules.

  

  !!!!!!!!!!!!!!! COMPUTE WEIGHTS & C6AIM !!!!!!!!!!!!!!!!

  ! Retriving data from results.tag
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
  

  if( debug )then
     call openfile(chrgfile,'write')
     do i = 1,natom
        write(fiit(chrgfile),'(A3,F10.4)') coords(i)%atom_type, Ni(i)
     end do
     call closefile(fiit(chrgfile))
  end if


  if( Ni_size /= natom ) stop 'ERROR: check 1' ! Controls
  
  ! Compute C6AIM
  allocate(C6aim(natom,2),C6free(natom,2),Zaim(natom),polar(natom),rvdw(natom))

  do i = 1,natom                                                  ! Some useful arrays:
     Zaim(i) = atomdata(findatom(coords(i)%atom_type))%Z          ! - Z_i of atoms in xyz file
     C6free(i,1) = atomdata(findatom(coords(i)%atom_type))%TSC6   ! - C6free by TS for "   "
     C6free(i,2) = atomdata(findatom(coords(i)%atom_type))%D3C6   ! - C6free by 3D for "   "
     polar(i) = atomdata(findatom(coords(i)%atom_type))%polarizability / BohrAngst**3 ! -polar for    "   "
     rvdw(i) = atomdata(findatom(coords(i)%atom_type))%vdWr / BohrAngst ! - Bondi radii
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
     Hi = IsHAtom(i)
     atom2: do j = i+1,natom
        Hj = IsHAtom(j)
        Rab = dist(coords(i)%coord,coords(j)%coord)/BohrAngst
        
        if( TTdf ) then
!           write(0,*) 'TTdf TRUE' ! debug
           if( fixb0 ) then

!           write(0,*) 'b0fixed TRUE' ! debug
              
              damp = fdamp(b0,Rab)
              hhrep = hCor(A,b0,Rab)
           else
!           write(0,*) 'b0fixed FAlse' ! debug
              ba = basym(b0,polar(i),Ni(i),Zaim(i))
              bb = basym(b0,polar(j),Ni(j),Zaim(j))
              bab = bmix(ba,bb)
              damp = fdamp(bab,Rab)
              hhrep = hCor(A,bab,Rab)
           end if
        elseif( Grd ) then
!          Rab0 = cubsum(rvdw(i),rvdw(j))
           Rab0 = rvdw(i)+rvdw(j)
           damp = wy2( b0,A,Rab,Rab0)
!           hhrep = hCor(A,bab,Rab)
        elseif( GrTTd ) then
           Rab0 = cubsum(rvdw(i),rvdw(j))
           damp =  GrTTfd(A,b0,Rab,Rab0)

        end if

!        print*, i,j,Rab,bb,ba,damp ! debug
!        print*, Rab, i ,j,coords(i)%coord,coords(j)%coord ! debug
!        stop ! debug 
           
!!$        Hctrl = 0
!!$        do l = 1,2
!!$           if( coords(i)%atom_type == HAtom(l) ) Hctrl = Hctrl + 1
!!$           if( coords(j)%atom_type == HAtom(l) ) Hctrl = Hctrl + 1
!!$           write(0,*) coords(j)%atom_type, coords(i)%atom_type,Hctrl ! debug
!!$           if( Hctrl == 2 ) write(0,*) "Hydrogen Atom Couple" ! debug
!!$           if( Hctrl > 2 )  write(0,*) "Huston, we have a problem!!" ! debug
!!$        end do
!        write(0,*) coords(j)%atom_type, coords(i)%atom_type,Hi,Hj ! debug
!        if( Hi .and. Hj ) write(0,*) "Hydrogen Atom Couple" ! debug


        TSvsD3loop : do k = 1,2  ! 1 => TS ; 2 => D3
           C6ab(1) = 2 * C6aim(i,k) * C6aim(j,k) / &
                ( C6aim(i,k) + C6aim(j,k) )                  ! REL A
           C6ab(2) = 2 * C6aim(i,k) * C6aim(j,k) / &
                ( ( polar(j) / polar(i) ) * C6aim(i,k) + &
                ( polar(i) / polar(j) ) * C6aim(j,k) )       ! REL B

           E(:,k) =  E(:,k) - damp * C6ab(:) / Rab**6   ! Dispersion Energy

           if( Hi .and. Hj ) then
!              print*, 'HH' ! debug
!              print*, "E before: ", E(:,k) ! debug
              E(:,k) = E(:,k) + hhrep  ! HH-repulsion correction
!              print*, "E after: ", E(:,k) ! debug
           end if

        end do TSvsD3loop
     end do atom2
  end do atom1

  
  ! E(1,:) => Energy of relationa A
  ! E(2,:) => Energy of relationa B
  ! E(:,1) => Energy from TS
  ! E(:,2) => Energy from D3

! The results have to be in Hartree to be consistent with DFTB+. 
  write(*,*) '#  ',inputcoofile
  write(*,'("A-TS  ",f20.12)') E(1,1)
  write(*,'("A-D3  ",f20.12)') E(1,2)
  write(*,'("B-TS  ",f20.12)') E(2,1)
  write(*,'("B-D3  ",f20.12)') E(2,2)

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
  
  real(kr) FUNCTION wy2(d,sr,r,r0)
    ! Grimme damping function 
    ! Material Transaction vol. 50 pagg 1664-1670 Corrected with Stephan's supporting information
    IMPLICIT NONE
    real(kr),intent(IN) :: d,sr,r,r0
    
    wy2 = 0.5 * ( 1 + tanh( 0.5*d* (r / ( sr*r0 ) -1) ) )
    
  END FUNCTION wy2
    

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

  real(kr) FUNCTION hCor(A,b,R)
    IMPLICIT NONE
    real(kr),intent(IN) :: A,b,R
    
    hCor = A * exp( -b * R )
!    print*, "Inside Function: ",hcor, A, b ,R ! debug
    
  END FUNCTION hCor

  real(kr) FUNCTION cubsum(a,b)
    IMPLICIT NONE
    real(kr),intent(IN) :: a,b
    
    cubsum = ( a**3 + b**3 ) / ( a**2 + b**3 )

  END FUNCTION cubsum

  !Mixed Fermi+TangToennies DF
  real(kr) FUNCTION  GrTTfd(a,b,R,R0)
    IMPLICIT NONE
    real(kr),intent(IN) :: a,b,R,R0
    print*, a

    GrTTfd = 0.5*( 1 + tanh( 23.0d0 * ( R / ( a * R0 ) - 1 ) ) ) * fdamp(b,R)


  END FUNCTION GrTTfd

END PROGRAM SCC_Disp
  
