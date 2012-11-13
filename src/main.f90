PROGRAM SCC_Disp

  ! Damping function:
  ! from: JCTC: Stephan Steinmann 2009
  ! (no "b" factor")
  ! 
  ! No Hbond correction
  !
  
  USE precision
  USE utils
  USE read_tag_dftbp
  USE read_xyz
  USE read_atomdata
  USE file_tools
  USE read_input  ! TODO better
  USE parameters, only : BohrAngst,HartKcalMol
  IMPLICIT NONE

  integer(ki) :: natom
  real(kr) :: E,C6ab,Rab0,Rab,ba,bb,bab,damp,hhrep

  type(xyz_coords),pointer :: mol1(:),mol2(:)

  integer(ki),pointer  :: Ni_size   ! each time that get_tag_data is called!!
  real(kr),dimension(:),pointer :: Ni

  integer(ki),dimension(:),allocatable:: Zaim
  real(kr),dimension(:),allocatable :: C6aim,C6free
  real(kr),dimension(:),allocatable :: polar,rvdw

  integer(ki) :: i,j

  logical :: Hi,Hj

  logical :: debug
  character(kch),parameter :: chrgfile = 'chargefile.check'
  character(kch),parameter :: coordfile = 'coordinate.check'
  character(kch),parameter :: C6aimfile = 'C6aim.check'
  character(kch),parameter :: C6last = 'C6last.check'
  character(kch),parameter :: energydbg = 'energyloop.check'
  character(kch),parameter :: distances = 'distances.check'

  call read_stdin

  if( debugflag == 'UP' ) debug=.true.

  call get_tag_data('atomic_charges',inputtagfile) 
  Ni => matrix(:,1,1)    ! WARNING:: If you recall 'get_tag_data' tag matrix change its values
  Ni_size => ifrmt(1)
  
  ! Retrieve data form atomdata.data
  call get_atomdata(atomdatafile) ! Now I can use the atomdata array
  
  ! Retrive name of atoms
  call get_coords(inputcoofile,natom)
  

  if( Ni_size /= natom ) stop 'ERROR: check 1' ! Controls
  

  ! => Initialize arrays
  allocate(C6aim(natom),C6free(natom),Zaim(natom),polar(natom),rvdw(natom))
  do i = 1,natom                                                  ! Some useful arrays:
     Zaim(i) = atomdata(findatom(coords(i)%atom_type))%Z          ! - Z_i of atoms in xyz file
     C6free(i) = atomdata(findatom(coords(i)%atom_type))%D3C6   ! - C6free by 3D for "   "
     polar(i) = atomdata(findatom(coords(i)%atom_type))%polarizability / BohrAngst**3 ! -polar for    "   "
     rvdw(i) = atomdata(findatom(coords(i)%atom_type))%vdWr / BohrAngst ! - Bondi radii
     Ni(i) = Ni(i) + atomdata(findatom(coords(i)%atom_type))%incharge
  end do
 
  if( debug )then
     call openfile(chrgfile,'replace')
     call openfile(coordfile,'replace')
     write(fiit(chrgfile),*) '# ATOM_TYPE  MULL_POP   INCHARGE   Z   RvdW  C6free'
     write(fiit(coordfile),*) natom
     write(fiit(coordfile),*)
     do i = 1,natom
        write(fiit(chrgfile),'(A3,2F15.5,I10,2F20.5)') coords(i)%atom_type, Ni(i),atomdata(findatom(coords(i)%atom_type))%incharge,Zaim(i),rvdw(i),C6free(i)
        write(fiit(coordfile),'(A3,3F15.5)') coords(i)%atom_type, coords(i)%coord
     end do
     call closefile(fiit(chrgfile))
     call closefile(fiit(coordfile))
  end if


  

  ! => C6AIM
  C6aim(:) = (Ni(:) / dble(Zaim(:)))**2 * C6free(:)                  ! D3
  
  if( debug )then
     call openfile(C6aimfile,'replace')
     write(fiit(C6aimfile),*) '# ATOM_TYPE   MULL_POP   Z   N/Z  (N/Z)**2  C6aim   C6free'
     do i = 1,natom
        write(fiit(C6aimfile),'(A3,F12.4,I6,4F12.4)') coords(i)%atom_type,Ni(i),Zaim(i),(Ni(i) / dble(Zaim(i))),(Ni(i) / dble(Zaim(i)))**2,C6aim(i),C6free(i)
!        write(fiit(C6aimfile),*) coords(i)%atom_type,Ni(i),Zaim(i),(Ni(:) / dble(Zaim(:))),(Ni(:) / dble(Zaim(:)))**2,C6aim(i),C6free(i)
     end do
     call closefile(fiit(C6aimfile))
  end if


  ! ==> STARTING ATOMIC LOOP <== !
  if( debug ) call openfile(energydbg,'replace')
  if( debug ) write(fiit(energydbg),*) '# AT_TP(i)    AT_TP(j)    Rab    Rab0    damp   C6(i)   C6(j)  C6ab    E   Etot'
  if( debug ) call openfile(distances,'replace')
  if( debug ) write(fiit(distances),*) '# AT_TP(i)    AT_TP(j)    xi   yi   zi   xj   yj    zj Rab'

  E = 0.0
  atom1: do i = 1,natom/2
     Hi = IsHAtom(i) ! To use hhrep
     atom2: do j = natom/2+1,natom
        Hj = IsHAtom(j) ! To use hhrep
        
        Rab = dist(coords(i)%coord,coords(j)%coord)/BohrAngst
        
        Rab0 = cubsum(rvdw(i),rvdw(j))
        !Rab0 = rvdw(i)+rvdw(j)
        damp =  FdTTdf(Rab,Rab0)
!        write(77,*) damp
!        damp = 1.d0
        hhrep = 0.0d0

        
        ! => Mixing C6aim (REL. A)
        C6ab = 2 * C6aim(i) * C6aim(j) / ( C6aim(i) + C6aim(j) )
        
        
        E =  E - damp * C6ab / Rab**6   ! Dispersion Energy
           
        if( debug ) write(fiit(energydbg),'(2A3,2X,10(X,F15.5))') coords(i)%atom_type, coords(j)%atom_type, Rab, Rab0, damp, C6aim(i), C6aim(j), C6ab, damp * C6ab / Rab**6, E
        if( debug ) write(fiit(distances),'(2A3,2X,7(X,F15.5))') coords(i)%atom_type, coords(j)%atom_type, coords(i)%coord, coords(j)%coord, Rab
        
        if( Hi .and. Hj ) then
           E = E + hhrep  ! HH-repulsion correction
        end if

     end do atom2
  end do atom1

  if( debug ) call closefile(energydbg)
  if( debug ) call closefile(distances)

  if (  debug )then
     call openfile(c6last,'replace')
     write(fiit(c6last),*) '#ATOM_TYPE     C6free    C6aim'
     do i = 1,natom
        write(fiit(c6last),'(A3,2F10.4)') coords(i)%atom_type, C6free(i), C6aim(i)
     end do
     call closefile(fiit(c6last))
  end if

  
  write(*,'(f20.12)') E*HartKcalMol

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
  
  real(kr) FUNCTION hCor(A,b,R)
    IMPLICIT NONE
    real(kr),intent(IN) :: A,b,R
    
    hCor = A * exp( -b * R )
!    print*, "Inside Function: ",hcor, A, b ,R ! debug
    
  END FUNCTION hCor

  real(kr) FUNCTION cubsum(a,b)
    IMPLICIT NONE
    real(kr),intent(IN) :: a,b
    
    cubsum = 2*( a**3 + b**3 ) / ( a**2 + b**2 )

  END FUNCTION cubsum

  real(kr) FUNCTION  FdTTdf(R,R0)
    IMPLICIT NONE
    real(kr) :: a,b,s,R,R0,TT,Fd

    CALL read_parameters(b,a,s)
    
!    b0 = 10
!    a=1
!    s=23

    Fd = 0.5*( 1.d0 + tanh( s * ( R / ( a * R0 ) - 1.d0 ) ) )
    TT = 1.d0 - ( exp( -b * R ) * (1.d0 + b*R + (b*R)**2.d0/2.d0 + (b*R)**3.d0/6.d0 + (b*R)**4.d0/24.d0 + (b*R)**5.d0/120.d0 + (b*R)**6.d0/720.d0) )
    FdTTdf = Fd * TT
!    write(88,*) Gr,TT,GrTTfd
    write(88,*) b,a,s

  END FUNCTION FdTTdf


  SUBROUTINE read_parameters(a,b,c)
    IMPLICIT NONE
    real(kr),intent(OUT) :: a,b,c
    integer(ki) :: err
    character(kch) :: filename='parameters.dat'

    err = 0
    call openfile(filename,'read')
    read(fiit(filename),*,iostat=err) a
    read(fiit(filename),*,iostat=err) b
    read(fiit(filename),*,iostat=err) c
    call closefile(filename)
    if( err /= 0 ) call die('ERROR: reading parameters.dat')

  END SUBROUTINE read_parameters
    

  SUBROUTINE die(msg)
    IMPLICIT NONE
    character(*),intent(IN) :: msg
    
    write(0,*) msg
    stop
    
  END SUBROUTINE die

END PROGRAM SCC_Disp
  
