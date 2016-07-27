PROGRAM dDMC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!                         _ ___  __  __  ___                      !!!
!!!                      __| |   \|  \/  |/ __|                     !!!
!!!                     / _` | |) | |\/| | (__                      !!!
!!!                     \__,_|___/|_|  |_|\___|                     !!!
!!!                                                                 !!!
!!! Provide a dispersion correction based on Mulliken charges.      !!!
!!! Read DOI: 10.1002/qua.24887 for details.                        !!!
!!!                                                                 !!!
!!! Contact Riccardo for any question about the code or for         !!!
!!! signaling bugs at riccardo.petraglia@gmail.com                  !!!
!!!                                                                 !!!
!!! This code is also available at:                                 !!!
!!!  * https://github.com/grhawk/dDMC/                              !!!
!!!                                                                 !!!
!!! Look the examples for instructions.                             !!!
!!!                                                                 !!!
!!! As regression test compare what you obtain in the ddmc.tag      !!!
!!! file with the ddmc.tag.chk file. A script should be available   !!!
!!! for that but it is not! ;)                                      !!!
!!!                                                                 !!!
!!! Version 1.0                                                     !!!
!!!                                                                 !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  ! The program works with atomic units but print 
  ! the final energy and forces are in kcal/mol and kcal/mol/ang
  
  USE precision
  USE utils
  USE read_tag_dftbp
  USE read_xyz
  USE read_atomdata
  USE file_tools
  USE read_input
  USE dampingfunctions
  USE read_gf
  USE iosockets
  USE parameters, only : BohrAngst,HartKcalMol
  IMPLICIT NONE

  integer(ki) :: natom
  real(kr) :: E,C6ab,Rab0,Rab,damp,hhrep,dfpr

  type(xyz_coords),pointer :: mol1(:),mol2(:)

  integer(ki),pointer  :: Ni_size   ! each time that get_tag_data is called!!
  real(kr),dimension(:),pointer :: Ni

  integer(ki),dimension(:),allocatable:: Zaim
  real(kr),dimension(:),allocatable :: C6aim,C6free
  real(kr),dimension(:),allocatable :: polar,rvdw,basym_ii
  real(kr),dimension(:,:),allocatable :: grad
  real(kr8),dimension(:,:),allocatable :: Coord0
  real(kr8),dimension(3,3)             :: latvec = 0.0  ! Lattice vector

  integer(ki) :: i,j

  logical :: Hi,Hj

  logical :: debug = .false.
  character(kch),parameter :: chrgfile = 'chargefile.check'
  character(kch),parameter :: coordfile = 'coordinate.check'
  character(kch),parameter :: C6aimfile = 'C6aim.check'
  character(kch),parameter :: C6last = 'C6last.check'
  character(kch),parameter :: energydbg = 'energyloop.check'
  character(kch),parameter :: distances = 'distances.check'
  character(kch),parameter :: dampingfunc = 'damping.check'
  character(kch),parameter :: bvalues = 'basym.check'
  character(kch),parameter :: excelfile = 'excelfile.check'
  character(kch),parameter :: graddbg = 'gradloop.check'

  character(kch),parameter :: tagfile = 'dDMC.tag'
!  integer(ki) :: dftype

  ! Socket stuff
  LOGICAL                       :: tSocket = .true.
  INTEGER                       :: socketPort = 0
  type(socketData),pointer      :: pSocket
  character(lc)                 :: SocketHost
  integer                       :: SocketVerb
  integer                       :: protocol
                                         
  real(kr) :: totalStress(3,3) = 0.0
  real(kr) :: cellVol = 0.0
  
  integer :: iRid, iORid                  !* Socket ID
  
  ! rid ... The index of the request, i.e. the replica that the force calculation is for (integer corresponding to the bead index)


  ! Read the input
  call read_stdin
  
  ! Initialize some stuff
  if( debugflag == 'UP' ) debug=.true.
  CALL initialize_dfmodule(dftype, dfparamsfile, method)
  if( dfprint ) CALL printDf()

  ! Print announce if debug is up
  if( debug ) call starting_program_announce

  ! Retrieve data from atomdata.data
  call get_atomdata(atomdatafile)
  
  ! Retrive name of atoms
  call get_coords(inputcoofile,natom)

  select case( tagtype )
  case ('dftbp') 
     call dftbp_tag()
  case ('column')
     call read_charge_gf(inputtagfile) ! Read the charge from a column.
     Ni => pop
     Ni_size => npop
  case default
     call dftbp_tag()
  end select
  
  if( debug )then 
     write(0,*) 'Ni_size: ',Ni_size 
     write(0,*) 'natom: ',natom 
  end if

  if( Ni_size /= natom ) call die('ERROR: Number of charge is different from the number of atoms!') ! Controls
  
  ! Socket initialization
  tSocket = .true.
  socketPort = 0  
  SocketHost = "/tmp/ipi_F2_AP-dDMC"
  SocketVerb = 55
  protocol =  IPI_PROTOCOL1
  iRid = -1
  iORid = -2  
  if( tSocket )then
    write(*,*) "Initialising for socket communication to host ", trim(SocketHost)    
    if ( socketPort == 0 ) then
      call create( pSocket, nAtom, SocketHost, SocketVerb, protocol )
    else
      call create( pSocket, nAtom, SocketHost, SocketVerb, protocol, SocketPort)
    end if
  end if

  
  ! => Initialize arrays
  allocate(C6aim(natom),C6free(natom),Zaim(natom),polar(natom),rvdw(natom),basym_ii(natom),grad(natom,3), Coord0(3,natom))
  do i = 1,natom                                                                         ! Some useful arrays:
     Zaim(i) = atomdata(findatom(coords(i)%atom_type))%Z                                 ! - Z_i of atoms in xyz file
     C6free(i) = atomdata(findatom(coords(i)%atom_type))%D3C6                            ! - C6free by 3D for "   "
     polar(i) = atomdata(findatom(coords(i)%atom_type))%polarizability / BohrAngst**3    ! - polar for        "   "
     rvdw(i) = atomdata(findatom(coords(i)%atom_type))%vdWr / BohrAngst                  ! - Bondi radii
     basym_ii(i) = basym(polar(i),Ni(i),Zaim(i))                                         ! - bii of the damping function for atom i (it has to be multiplied by b0).
  end do

  ! => C6AIM
  C6aim(:) = (Ni(:) / dble(Zaim(:)))**2 * C6free(:)                  ! D3
  

  if( debug )then
     call openfile(bvalues,'replace')                                     
     call openfile(chrgfile,'replace')
     call openfile(coordfile,'replace')
     call openfile(C6aimfile,'replace')
     write(fiit(bvalues),*) '#ATOM_TYPE_i  b_asym  alpha  N  Z' 
     write(fiit(chrgfile),*) '# ATOM_TYPE  MULL_POP   INCHARGE   Z   RvdW  C6free'
     write(fiit(C6aimfile),*) '# ATOM_TYPE   MULL_POP   Z   N/Z  (N/Z)**2  C6aim   C6free'
     write(fiit(coordfile),*) natom
     write(fiit(coordfile),*)
     do i = 1,natom
        write(fiit(chrgfile),'(A3,2F15.5,I10,2F20.5)') coords(i)&
             &%atom_type, Ni(i),atomdata(findatom(coords(i)&
             &%atom_type))%incharge,Zaim(i),rvdw(i),C6free(i) 
        write(fiit(coordfile),'(A3,3F15.5)') coords(i)%atom_type, coords(i)%coord
        write(fiit(bvalues),'(1A3,3F8.3,I5)') coords(i)%atom_type, basym_ii(i), polar(i), Ni(i), Zaim(i)
        write(fiit(C6aimfile),'(A3,F12.4,I6,4F12.4)') &
             & coords(i)%atom_type,Ni(i),Zaim(i), (Ni(i) / dble(Zaim(i))),&
             & (Ni(i) / dble(Zaim(i)))**2,C6aim(i),C6free(i)
     end do
     call closefile(fiit(chrgfile))
     call closefile(fiit(coordfile))
     call closefile(bvalues)
     call closefile(fiit(C6aimfile))
  end if


  ! ==> STARTING ATOMIC LOOP <== !
  if( debug ) call openfile(energydbg,'replace')
  if( debug ) write(fiit(energydbg),*) '#AT_TP(i)    AT_TP(j)    Rab    Rab0    damp   C6(i)   C6(j)  C6ab    E_ij   Etot'
  if( debug ) call openfile(graddbg,'replace')
  if( debug ) write(fiit(graddbg),*) '#AT_TP(i)    AT_TP(j)   i  j    Rab    Rab0  Rvdw(i)  Rvdw(j)  dfp    C6ab    Fij'
  if( debug ) call openfile(distances,'replace')
  if( debug ) write(fiit(distances),*) '# AT_TP(i)    AT_TP(j)    xi   yi   zi   xj   yj    zj Rab'
  if( debug ) call openfile(dampingfunc,'replace')
  if( debug ) write(fiit(dampingfunc),*) '#ATOM_TYPE_i ATOM_TYPE_j   R0  R  b0  a  s  basymi  basymj   bij  Fd   TT   df'
  if( debug ) call openfile(excelfile,'replace')
  if( debug ) write(fiit(excelfile),*) '#AT_TP(i)  AT_TP(j) Rab  Rvdw(&
       &i) Rvdw(j)  Rab0  alpha(i)  alpha(j)  C6free(i)  C6free(j)  Za&
       &im(i)  Zaim(j)  Ni(i)  Ni(j)  C6aim(i)  C6aim(j)  C6ab  basym(&
       &i)  basym(j)  damp  E  Etot' 

mainloop: do while (.true.)

  call retrieve(pSocket, Coord0, latvec, iRid) ! iRid is the index of the request, i.e. the replica that the force calculation is for 
  
  do i = 1,natom
    coords(i)%coord(1) = Coord0(1,i)
    coords(i)%coord(2) = Coord0(2,i)
    coords(i)%coord(3) = Coord0(3,i)
  end do
  
  ! Can be deleted
  if ( iRid >= 0 ) then
    if ( iRid /= iORid ) then ! iORid = -2 in this version of the program (Revision 4913)
!       qInput = q0
      iORid = iRid
    end if
  end if

  E = 0.0d0
  atom1: do i = 1,natom
     Hi = IsHAtom(i) ! To use hhrep
     atom2: do j = i+1,natom
        Hj = IsHAtom(j) ! To use hhrep
        
        Rab = dist(coords(i)%coord,coords(j)%coord)/BohrAngst
        
!        Rab0 = cubsum(rvdw(i),rvdw(j))
        Rab0 = rvdw(i)+rvdw(j)
        
        if(debug) then
          damp =  df(coords(i)%atom_type,coords(j)%atom_type,basym_ii(i),basym_ii(j),Rab,Rab0,dampingfunc)
        else
          damp =  df(coords(i)%atom_type,coords(j)%atom_type,basym_ii(i),basym_ii(j),Rab,Rab0)
        end if
        
        ! => Mixing C6aim (REL. A)
        C6ab = 2 * C6aim(i) * C6aim(j) / ( C6aim(i) + C6aim(j) )
        
        
        E =  E - damp * C6ab / Rab**6   ! Dispersion Energy
        
        if( debug ) write(fiit(energydbg),'(2A3,2X,10(X,F15.5))') &
             & coords(i)%atom_type, coords(j)%atom_type, Rab, Rab0, &
             & damp, C6aim(i), C6aim(j), C6ab, -damp * C6ab / Rab**6, E
        if( debug ) write(fiit(excelfile),'(2A3,2X,20(X,F15.5))') &
             & coords(i)%atom_type, coords(j)%atom_type, Rab, rvdw(i)&
             &, rvdw(j), Rab0, polar(i), polar(j), C6free(i),&
             & C6free(j), real(Zaim(i)),real(Zaim(j)), Ni(i), Ni(j),&
             & C6aim(i), C6aim(j), C6ab, basym_ii(i), basym_ii(j),&
             & damp, damp * C6ab / Rab**6, E 
        if( debug ) write(fiit(distances),'(2A3,2X,7(X,F15.5))') &
             & coords(i)%atom_type, coords(j)%atom_type, &
             & coords(i)%coord, coords(j)%coord, Rab
        
        if( Hi .and. Hj ) then
           hhrep = 0d0
           E = E + hhrep  ! HH-repulsion correction
        end if

     end do atom2
     
  end do atom1


!    write(*,*) 'ATOM i     ATOM j     Grad'
  if (readgradflag == 'UP') then 
     grad(:,:) = 0.d0
     do i = 1,natom
        do j = 1,natom
           if (j /= i) then
              Rab = dist(coords(i)%coord,coords(j)%coord)/BohrAngst
              Rab0 = rvdw(i)+rvdw(j)
              dfpr = dfp(basym_ii(i),basym_ii(j),Rab,Rab0)
              C6ab = 2 * C6aim(i) * C6aim(j) / ( C6aim(i) + C6aim(j) )
              grad(i,:) = grad(i,:) - (dfpr * C6ab * (coords(i)%coord - coords(j)%coord)/BohrAngst/Rab)
              if( debug ) write(fiit(graddbg),'(2A3,2X,2I4,9(X,E20.12))') &
                   & coords(i)%atom_type, coords(j)%atom_type, i, j, Rab, Rab0, rvdw(i), rvdw(j), &
                   & dfpr, C6ab, dfpr * C6ab * (coords(i)%coord - coords(j)%coord)/BohrAngst/Rab
!                write(*,'(2I4,3E30.10)') i,j,grad(i,:)
           end if
        end do
     end do
  end if
  
  if( debug ) call closefile(energydbg)
  if( debug ) call closefile(distances)
  if( debug ) call closefile(dampingfunc)
  if( debug ) call closefile(excelfile)
  if( debug ) call closefile(graddbg)

  if ( debug )then
     call openfile(c6last,'replace')
     write(fiit(c6last),*) '#ATOM_TYPE     C6free    C6aim'
     do i = 1,natom
        write(fiit(c6last),'(A3,2F10.4)') coords(i)%atom_type, C6free(i), C6aim(i)
     end do
     call closefile(fiit(c6last))
  end if

  write(*,'("Energy:  " f20.12 "    kcal/mol")') E*HartKcalMol ! E should be Hartree (atomic) units and HartKcalMol is conversion factor from Hart to KcalMol
  if (readgradflag == 'UP') then 
     write(*,'("Forces:")')
     write(*,'("Atom    x     y    z")')
     write(*,'(A3, 3g20.12)') (coords(i)%atom_type,grad(i,:)*27.211399000000004/0.529177249, i = 1,natom) ! 27.211399000000004/0.529177249 conversion to eV per Angstrom, 0.529177249: Bohr ~ Angstrom conversion
  end if
  
 ! E already in atomic units
!  print *,'SHAPE1', SHAPE(grad(i,:),(/3,natom/)))
!  print *,'SHAPE2', SHAPE(RESHAPE(-grad(i,:),(/3,natom/)))
 
 call send(pSocket, E, RESHAPE(-grad(:,:),(/3,natom/)), totalStress*cellVol)

  CALL openfile(tagfile,'replace')
  write(fiit(tagfile),*) 'correction_energy     1'
  write(fiit(tagfile),'(E30.20)') E
  if (readgradflag == 'UP') then
    write(fiit(tagfile),*) 'forces    3,',natom
    do i = 1,natom
      write(fiit(tagfile),'(3E30.20)') (grad(i,j), j=1,3)
    end do
  end if
  CALL closefile(tagfile)

!   write(0,*) NEW_LINE('C'),'The Programs Ends Correctly!'
  
end do mainloop

CONTAINS
  
  SUBROUTINE dftbp_tag()
    call get_tag_data('atomic_charges',inputtagfile) 
    Ni => matrix(:,1,1)    ! WARNING:: If you recall 'get_tag_data' tag matrix change its values
    Ni_size => ifrmt(1)
    do i = 1,natom
       Ni(i) = Ni(i) + atomdata(findatom(coords(i)%atom_type))%incharge  ! - Population for atom i (calculated as population on the outest shell + population in the inner shells)
    end do
  END SUBROUTINE dftbp_tag


  real(kr) FUNCTION basym(alpha,N,Z)
    IMPLICIT NONE
    real(kr),intent(IN) :: alpha,N
    integer(ki),intent(IN) :: Z
    
    basym = ( 1 / alpha )**(1./3.) * ( real(Z) / N )**(1./3.)

  END FUNCTION basym
  

  real(kr) FUNCTION cubsum(a,b)
    IMPLICIT NONE
    real(kr),intent(IN) :: a,b
    
    cubsum = 2*( a**3 + b**3 ) / ( a**2 + b**2 )

  END FUNCTION cubsum

    

END PROGRAM dDMC

   
