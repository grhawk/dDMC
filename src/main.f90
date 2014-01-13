PROGRAM SCC_Disp

  ! Damping function:
  ! from: JCTC: Stephan Steinmann 2009
  ! (with "b" factor depending on atom couple)
  ! 
  ! No Hbond correction
  !

  ! The program works with atomic units but print 
  ! the final energy in kcal/mol
  
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
  real(kr) :: E,C6ab,Rab0,Rab,damp,hhrep

  type(xyz_coords),pointer :: mol1(:),mol2(:)

  integer(ki),pointer  :: Ni_size   ! each time that get_tag_data is called!!
  real(kr),dimension(:),pointer :: Ni

  integer(ki),dimension(:),allocatable:: Zaim
  real(kr),dimension(:),allocatable :: C6aim,C6free
  real(kr),dimension(:),allocatable :: polar,rvdw,basym_ii

  integer(ki) :: i,j

  logical :: Hi,Hj

  logical :: debug
  character(kch),parameter :: chrgfile = 'chargefile.check'
  character(kch),parameter :: coordfile = 'coordinate.check'
  character(kch),parameter :: C6aimfile = 'C6aim.check'
  character(kch),parameter :: C6last = 'C6last.check'
  character(kch),parameter :: energydbg = 'energyloop.check'
  character(kch),parameter :: distances = 'distances.check'
  character(kch),parameter :: dampingfunc = 'damping.check'
  character(kch),parameter :: bvalues = 'basym.check'
  character(kch),parameter :: excelfile = 'excelfile.check'

  character(kch),parameter :: tagfile = 'dDMC.tag'


  call read_stdin
  if( dfprint ) CALL printDf()

  if( debugflag == 'UP' ) debug=.true.

  call get_tag_data('atomic_charges',inputtagfile) 
  Ni => matrix(:,1,1)    ! WARNING:: If you recall 'get_tag_data' tag matrix change its values
  Ni_size => ifrmt(1)
  
  ! Retrieve data from atomdata.data
  call get_atomdata(atomdatafile) ! Now I can use the atomdata array
  
  ! Retrive name of atoms
  call get_coords(inputcoofile,natom)
  

  if( Ni_size /= natom ) stop 'ERROR: check 1' ! Controls
  

  ! => Initialize arrays
  allocate(C6aim(natom),C6free(natom),Zaim(natom),polar(natom),rvdw(natom),basym_ii(natom))
  do i = 1,natom                                                                         ! Some useful arrays:
     Zaim(i) = atomdata(findatom(coords(i)%atom_type))%Z                                 ! - Z_i of atoms in xyz file
     C6free(i) = atomdata(findatom(coords(i)%atom_type))%D3C6                            ! - C6free by 3D for "   "
     polar(i) = atomdata(findatom(coords(i)%atom_type))%polarizability / BohrAngst**3    ! - polar for        "   "
     rvdw(i) = atomdata(findatom(coords(i)%atom_type))%vdWr / BohrAngst                  ! - Bondi radii
     Ni(i) = Ni(i) + atomdata(findatom(coords(i)%atom_type))%incharge                    ! - Population for atom i (calculated as population on the outest shell + population in the inner shells)
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
  if( debug ) call openfile(distances,'replace')
  if( debug ) write(fiit(distances),*) '# AT_TP(i)    AT_TP(j)    xi   yi   zi   xj   yj    zj Rab'
  if( debug ) call openfile(dampingfunc,'replace')
  if( debug ) write(fiit(dampingfunc),*) '#ATOM_TYPE_i ATOM_TYPE_j   R0  R  b0  a  s  basymi  basymj   bij  Fd   TT   dfR0 '
  if( debug ) call openfile(excelfile,'replace')
  if( debug ) write(fiit(excelfile),*) '#AT_TP(i)  AT_TP(j) Rab  Rvdw(&
       &i) Rvdw(j)  Rab0  alpha(i)  alpha(j)  C6free(i)  C6free(j)  Za&
       &im(i)  Zaim(j)  Ni(i)  Ni(j)  C6aim(i)  C6aim(j)  C6ab  basym(&
       &i)  basym(j)  damp  E  Etot' 


  E = 0.0d0
  atom1: do i = 1,natom
     Hi = IsHAtom(i) ! To use hhrep
     atom2: do j = i+1,natom
        Hj = IsHAtom(j) ! To use hhrep
        
        Rab = dist(coords(i)%coord,coords(j)%coord)/BohrAngst
        
!        Rab0 = cubsum(rvdw(i),rvdw(j))
        Rab0 = rvdw(i)+rvdw(j)
        damp =  df(dftype,basym_ii(i),basym_ii(j),Rab,Rab0)
!        write(77,*) damp
!        damp = 1.d0

        
        ! => Mixing C6aim (REL. A)
        C6ab = 2 * C6aim(i) * C6aim(j) / ( C6aim(i) + C6aim(j) )
        
        
        E =  E - damp * C6ab / Rab**6   ! Dispersion Energy
           
        if( debug ) write(fiit(energydbg),'(2A3,2X,10(X,F15.5))') &
             & coords(i)%atom_type, coords(j)%atom_type, Rab, Rab0, &
             & damp, C6aim(i), C6aim(j), C6ab, damp * C6ab / Rab**6, E
        if( debug ) write(fiit(excelfile),'(2A3,2X,20(X,F15.5))')&
             & coords(i)%atom_type, coords(j)%atom_type, Rab, rvdw(i)&
             &, rvdw(j), Rab0, polar(i), polar(j), C6free(i),&
             & C6free(j), real(Zaim(i)),real(Zaim(j)), Ni(i), Ni(j),&
             & C6aim(i), C6aim(j), C6ab, basym_ii(i), basym_ii(j),&
             & damp, damp * C6ab / Rab**6, E 
        if( debug ) write(fiit(distances),'(2A3,2X,7(X,F15.5))') &
             & coords(i)%atom_type, coords(j)%atom_type, &
             & coords(i)%coord, coords(j)%coord, Rab
        
        if( Hi .and. Hj ) then
!           print*, 'Two H   ',coords(i)%atom_type,coords(j)%atom_type ! debug
           hhrep = hcor(Rab)
!           hhrep = 0d0
           E = E + hhrep  ! HH-repulsion correction
        end if

     end do atom2
  end do atom1

  if( debug ) call closefile(energydbg)
  if( debug ) call closefile(distances)
  if( debug ) call closefile(dampingfunc)
  if( debug ) call closefile(excelfile)

  if (  debug )then
     call openfile(c6last,'replace')
     write(fiit(c6last),*) '#ATOM_TYPE     C6free    C6aim'
     do i = 1,natom
        write(fiit(c6last),'(A3,2F10.4)') coords(i)%atom_type, C6free(i), C6aim(i)
     end do
     call closefile(fiit(c6last))
  end if

  
  write(*,'(f20.12)') E*HartKcalMol
  CALL openfile(tagfile,'replace')
  write(fiit(tagfile),*) 'correction_energy     1'
  write(fiit(tagfile),'(E30.20)') E*HartKcalMol
  CALL closefile(tagfile)
  

CONTAINS
  

  real(kr) FUNCTION bmix(bi,bj)
    IMPLICIT NONE
    real(kr),intent(IN) :: bi,bj

    bmix = 2 * bi * bj / (bi + bj)
    
  END FUNCTION bmix


  real(kr) FUNCTION basym(alpha,N,Z)
    IMPLICIT NONE
    real(kr),intent(IN) :: alpha,N
    integer(ki),intent(IN) :: Z
    
    basym = ( 1 / alpha )**(1./3.) * ( real(Z) / N )**(1./3.)
!    if (debug) write(fiit(bvalues),'(1A3,3F8.3,I5)') coords(i)%atom_type, basym, alpha, N, Z
!    print*, 'basym', b0,alpha,Z,N  ! debug
    
  END FUNCTION basym
  
  real(kr) FUNCTION hCor(R)
    IMPLICIT NONE
    real(kr),intent(IN) :: R
    real(kr) :: A,b!,dummy

!    CALL read_parameters(A,b,dummy)

    A = 1.49689540654504
    b = 1.67162251877518
    
    hCor = A * exp( -b * R )
    
!    print*, "Inside Function: ",hcor, A, b ,R ! debug
    
  END FUNCTION hCor

  real(kr) FUNCTION cubsum(a,b)
    IMPLICIT NONE
    real(kr),intent(IN) :: a,b
    
    cubsum = 2*( a**3 + b**3 ) / ( a**2 + b**2 )

  END FUNCTION cubsum

  real(kr) FUNCTION  df(type,basymi,basymj,R,R0)
    IMPLICIT NONE
    integer(ki),intent(IN) :: type
    real(kr),intent(IN) :: basymi,basymj,R,R0
    real(kr) :: a,b0,s,TT,Fd,bij,x,bx

!    b0 = 10
!    a=1
!    s=23
    select case (type)
    case (1)
       ! From the first report

       b0 = 1.5d0
       a = 0d0
       s=0d0
       
       bij = bmix(b0*basymi,b0*basymj)
       bx = bij*R
       
       TT = 1.d0 - &
            & ( exp( -bx ) * (1.d0 + bx + (bx)**2.d0/2.d0 + (bx)**3.d0/6.d0 + &
            & (bx)**4.d0/24.d0 + (bx)**5.d0/120.d0 + (bx)**6.d0/720.d0) )

       df = TT

       if(debug)write(fiit(dampingfunc),'(2A3,11F15.8)')coords(i)&
            &%atom_type,coords(j)%atom_type,R0,R,b0,a,s,basymi,basymj,bij&
            &,Fd,TT,df 
       
       
    case (2)
       !FdTTdf

       !CALL read_parameters(b0,a,s)

       b0 = 2.18206081886510
       a =  1.12451132211179
       s =  34.9266956797606
    
       bij = bmix(b0*basymi,b0*basymj)
       x = bij*R
       
       Fd = 0.5*( 1.d0 + tanh( s * ( x / ( a * R0 ) - 1.d0 ) ) )
       bx = Fd * bij * R
       
       TT = 1.d0 - &
            & ( exp( -bx ) * (1.d0 + bx + (bx)**2.d0/2.d0 + (bx)**3.d0/6.d0 + &
            & (bx)**4.d0/24.d0 + (bx)**5.d0/120.d0 + (bx)**6.d0/720.d0) )
       
       df = TT
       
!       if(debug)write(fiit(dampingfunc),*) 'DF Type: ',type

       if(debug)write(fiit(dampingfunc),'(2A3,11F15.8)')coords(i)&
            &%atom_type,coords(j)%atom_type,R0,R,b0,a,s,basymi,basymj,bij&
            &,Fd,TT,df 

    case (3)
       !Fd*TTdf

       !CALL read_parameters(b0,a,s)

       b0 = 2.18206081886510
       a =  1.12451132211179
       s =  34.9266956797606
    
       bij = bmix(b0*basymi,b0*basymj)
       x = bij*R
       
       Fd = 0.5*( 1.d0 + tanh( s * ( x / ( a * R0 ) - 1.d0 ) ) )
       bx = bij * R
       
       TT = 1.d0 - &
            & ( exp( -bx ) * (1.d0 + bx + (bx)**2.d0/2.d0 + (bx)**3.d0/6.d0 + &
            & (bx)**4.d0/24.d0 + (bx)**5.d0/120.d0 + (bx)**6.d0/720.d0) )
       
       df = TT*Fd
       
!       if(debug)write(fiit(dampingfunc),*) 'DF Type: ',type

       if(debug)write(fiit(dampingfunc),'(2A3,11F15.8)')coords(i)&
            &%atom_type,coords(j)%atom_type,R0,R,b0,a,s,basymi,basymj,bij&
            &,Fd,TT,df 
       
    case DEFAULT
       df = 1E10

    end select
    
!    TT = 1.d0 - &
!         & ( exp( -bij * R ) * (1.d0 + bij*R + (bij*R)**2.d0/2.d0 + (bij*R)**3.d0/6.d0 + &
!         & (bij*R)**4.d0/24.d0 + (bij*R)**5.d0/120.d0 + (bij*R)**6.d0/720.d0) )
          

!    FdTTdf = Fd * TT
!    FdTTdf = Fd
!    FdTTdf = 1.0d0


  END FUNCTION Df
  

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


  SUBROUTINE printDf()
    ! To print all the available damping function at once.
    IMPLICIT NONE
    integer(ki) :: i,j,k
    integer(kr) :: np
    real(kr) :: start_d,increment,end_d,d
    real(kr),allocatable :: df_value(:)
    integer(kr),parameter :: max_dftype = 10
    character(kch) :: format_string
    
    allocate(df_value(max_dftype))
    
    np = 10000
    start_d = 0.1d0
    end_d   = 8d0
    increment = (end_d - start_d)/np
    
    d = start_d
    do j = 1,np
       dftype: do i = 1,max_dftype
          if( df(i,1d0,1d0,start_d,1d9) > 1E9 ) exit dftype
          df_value(i) = df(i,0.5d0,0.5d0,d,1.9d0)
       end do dftype
       write(format_string,*) i
       format_string = '('//trim(format_string)//'F15.8)'
       write(*,format_string) d, (-df_value(k)/d**6*1000, k=1,i-1)
       d = d + increment
    enddo
    
    
    
    stop
    
  END SUBROUTINE printDf
  
END PROGRAM SCC_Disp

   
