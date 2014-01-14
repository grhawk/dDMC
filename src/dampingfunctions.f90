MODULE dampingfunctions

  USE precision
  USE file_tools
  USE read_input

  IMPLICIT NONE
  TYPE,private :: df_function
     integer(ki)    :: index
     character(kch) :: name
     integer(ki)    :: n_par
     real(kr),allocatable :: parameters(:)
  END type df_function
  
! PUBLIC
  logical :: read_param_from_file = .false.
!  integer(ki) :: dftype

! PRIVATE
  integer(ki),parameter :: df_num = 3 ! number of available damping function
  logical :: hh_read_param = .false., df_read_param = .false.
  type(df_function) :: df_functions(0:df_num)
  

  PRIVATE
  PUBLIC read_param_from_file, dftype
  PUBLIC initialize_dfmodule,df,printDf,hcor
  
CONTAINS
  
  SUBROUTINE initialize_dfmodule()
    CALL df_creation()
    
    if ( read_param_from_file ) then
!       write(*,*) 'dftype -> ',dftype
       if( dftype > 100 ) then
          dftype = dftype - 100
          CALL read_params(df_functions(0)%n_par,df_functions(0)%parameters)
       else
          CALL read_params(df_functions(dftype)%n_par,df_functions(dftype)%parameters)
       end if

    end if
  END SUBROUTINE initialize_dfmodule

  SUBROUTINE df_creation()
    IMPLICIT NONE
    integer(ki) :: type
    
    type = 0
    df_functions(type)%index = 0
    df_functions(type)%name  = 'HH'
    df_functions(type)%n_par = 2
    allocate(df_functions(type)%parameters(df_functions(type)%n_par))
    df_functions(type)%parameters = (/ 1.49689540654504, 1.67162251877518 /)
    
    type = 1
    df_functions(type)%index = 1
    df_functions(type)%name  = 'TT'
    df_functions(type)%n_par = 1
    allocate(df_functions(type)%parameters(df_functions(type)%n_par))
    df_functions(type)%parameters = (/ 1.5d0 /)
!                                        b0

    type = 2
    df_functions(type)%index = 2
    df_functions(type)%name  = 'TT(Fd(bijR)bijR)'
    df_functions(type)%n_par = 3
    allocate(df_functions(type)%parameters(df_functions(type)%n_par))
    df_functions(type)%parameters = (/ 2.18206081886510, 1.12451132211179, 34.9266956797606 /)
!                                          b0                   a               s

    type = 3
    df_functions(type)%index = 3
    df_functions(type)%name  = 'TT(bijR)*Fd(bijR)'
    df_functions(type)%n_par = 3
    allocate(df_functions(type)%parameters(df_functions(type)%n_par))
    df_functions(type)%parameters = (/ 2.18206081886510, 1.12451132211179, 34.9266956797606 /)
!                                          b0                   a               s

  END SUBROUTINE df_creation


  real(kr) FUNCTION hCor(R)
    IMPLICIT NONE
    real(kr),intent(IN) :: R
    real(kr) :: A,b!,dummy

    A = df_functions(0)%parameters(1)
    b = df_functions(0)%parameters(2)

    hCor = A * exp( -b * R )
    
!    print*, "Inside Function: ",hcor, A, b ,R ! debug
    
  END FUNCTION hCor



  real(kr) FUNCTION  df(basymi,basymj,R,R0)
    IMPLICIT NONE
    real(kr),intent(IN) :: basymi,basymj,R,R0
    real(kr) :: TT,Fd,bij,x,bx

    select case (dftype)
    case (1)
       ! From the first report
       
       bij = bmix(df_functions(dftype)%parameters(1)*basymi, &
            &     df_functions(dftype)%parameters(1)*basymj)
       bx = bij*R
       
       TT = 1.d0 - &
            & ( exp( -bx ) * (1.d0 + bx + (bx)**2.d0/2.d0 + (bx)**3.d0/6.d0 + &
            & (bx)**4.d0/24.d0 + (bx)**5.d0/120.d0 + (bx)**6.d0/720.d0) )

       df = TT

       ! if(debug)write(fiit(dampingfunc),'(2A3,11F15.8)')coords(i)&
       !      &%atom_type,coords(j)%atom_type,R0,R,b0,a,s,basymi,basymj,bij&
       !      &,Fd,TT,df 
       
       
    case (2)
       !FdTTdf
!                            b0
       bij = bmix(df_functions(dftype)%parameters(1)*basymi, &
            & df_functions(dftype)%parameters(1)*basymj)
       x = bij*R
       
!       Fd = 0.5*( 1.d0 + tanh( s * ( x / ( a * R0 ) - 1.d0 ) ) )
       Fd = 0.5*( 1.d0 + tanh( df_functions(dftype)%parameters(3) * &
            & ( x / ( df_functions(dftype)%parameters(2) * R0 ) - 1.d0 ) ) )
       bx = Fd * bij * R
       
       TT = 1.d0 - &
            & ( exp( -bx ) * (1.d0 + bx + (bx)**2.d0/2.d0 + (bx)**3.d0/6.d0 + &
            & (bx)**4.d0/24.d0 + (bx)**5.d0/120.d0 + (bx)**6.d0/720.d0) )
       
       df = TT
       
!       if(debug)write(fiit(dampingfunc),*) 'DF Type: ',type

    case (3)
       !Fd*TTdf
    
       bij = bmix(df_functions(dftype)%parameters(1)*basymi, &
            & df_functions(dftype)%parameters(1)*basymj)
       x = bij*R
       
       Fd = 0.5*( 1.d0 + tanh( df_functions(dftype)%parameters(3) * &
            & ( x / ( df_functions(dftype)%parameters(2) * R0 ) - 1.d0 ) ) )

       bx = bij * R
       
       TT = 1.d0 - &
            & ( exp( -bx ) * (1.d0 + bx + (bx)**2.d0/2.d0 + (bx)**3.d0/6.d0 + &
            & (bx)**4.d0/24.d0 + (bx)**5.d0/120.d0 + (bx)**6.d0/720.d0) )
       
       df = TT*Fd
       
! !       if(debug)write(fiit(dampingfunc),*) 'DF Type: ',type

!        ! if(debug)write(fiit(dampingfunc),'(2A3,11F15.8)')coords(i)&
!        !      &%atom_type,coords(j)%atom_type,R0,R,b0,a,s,basymi,basymj,bij&
!        !      &,Fd,TT,df, type
       
    case DEFAULT
       df = 1E10
       
    end select

    ! if(debug)write(fiit(dampingfunc),'(2A3,11F15.8)')coords(i)&
    !      &%atom_type,coords(j)%atom_type,R0,R,b0,a,s,basymi,basymj,bij&
    !      &,Fd,TT,df,type
    
    
!    TT = 1.d0 - &
!         & ( exp( -bij * R ) * (1.d0 + bij*R + (bij*R)**2.d0/2.d0 + (bij*R)**3.d0/6.d0 + &
!         & (bij*R)**4.d0/24.d0 + (bij*R)**5.d0/120.d0 + (bij*R)**6.d0/720.d0) )
          

!    FdTTdf = Fd * TT
!    FdTTdf = Fd
!    FdTTdf = 1.0d0


  END FUNCTION Df
  


  SUBROUTINE read_params(n,parameters)
    IMPLICIT NONE
    real(kr),intent(OUT) :: parameters(:)
    integer(ki),intent(IN) :: n
    integer(ki) :: err,i
    character(kch) :: filename='parameters.dat'
    real(kr) :: tmp
    
!    write(*,*) n

    err = 0
    call openfile(filename,'read')
    do i = 1,n
       read(fiit(filename),*,iostat=err) parameters(i)
!       write(*,*) 'p -> ',parameters(i)
    end do
    call closefile(filename)
    if( err /= 0 ) call die('ERROR: reading parameters.dat')
    
  END SUBROUTINE read_params



  real(kr) FUNCTION bmix(bi,bj)
    IMPLICIT NONE
    real(kr),intent(IN) :: bi,bj

    bmix = 2 * bi * bj / (bi + bj)
    
  END FUNCTION bmix




  SUBROUTINE die(msg)
    IMPLICIT NONE
    character(*),intent(IN) :: msg
    
    write(0,*) msg
    stop
    
  END SUBROUTINE die


  SUBROUTINE printDf()
    ! To print all the available damping function at once.
    IMPLICIT NONE
    integer(ki) :: j,k
    integer(kr) :: np
    real(kr) :: start_d,increment,end_d,d
    real(kr),allocatable :: df_value(:)
    character(kch) :: format_string
    
    allocate(df_value(df_num))
    
    np = 10000
    start_d = 0.1d0
    end_d   = 8d0
    increment = (end_d - start_d)/np
    
    d = start_d
    do j = 1,np
       do dftype = 1,df_num
          df_value(dftype) = df(0.5d0,0.5d0,d,1.9d0)
!          write(*,*) df(0.5d0,0.5d0,d,1.9d0)
       end do
       write(format_string,*) df_num+1
       format_string = '('//trim(format_string)//'F15.8)'
       write(*,format_string) d, (-df_value(k)/d**6*1000, k=1,df_num)
!       write(*,*) d, (-df_value(k)/d**6*1000, k=1,df_num)
       d = d + increment
    enddo
    stop
        
  END SUBROUTINE printDf

END MODULE dampingfunctions

! PROGRAM test
!   ! To compile this module has testdf program
!   ! you need to decomment the line of the program part
!   ! and run in a shell the following comamnd:
!   ! $ gfortran precision.f90 file_tools.f90 dampingfunctions.f90  -o testdf
!   USE dampingfunctions
!   IMPLICIT NONE
  
!   read_param_from_file = .false.
!   read_param_from_file = .true.
!   dftype = 1
  
!   CALL initialize_dfmodule()
!   write(*,*) df(1d0,1d0,2d0,3d0) 
!   write(*,*) df(1d0,1d0,3d0,3d0) 
!   write(*,*) df(1d0,1d0,4d0,3d0) 

! END PROGRAM test