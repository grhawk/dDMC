!> Routines to make socket contact with an external code and
!! communicate data back and forward from DFTB+ to the external code.
module iosockets
#include "allocate.h"
#include "assert.h"
  use accuracy
  use message
  use f90sockets
  implicit none
  
  private
  
  !> type for the socket set-up and comms.
  type socketData
    private
    integer                :: verbosity = 0       ! level of verbosity
    integer                :: Socket              ! socket number
    integer                :: nAtom     = 0       ! expected number of atoms
    logical                :: tInit     = .false. ! Initialisation of variables
  end type socketData
  
  integer, parameter :: IPI_PROTOCOL1 = 1 ! Definition of the message
                                          ! format and data to
                                          ! exchange
  integer, parameter :: IPI_PROTOCOL2 = 2 ! Definition of the message
                                          ! format and data to
                                          ! exchange (replica indexes)
  integer, parameter :: IPI_MSGLEN  = 12 ! length of strings expected
                                         ! for i-pi messages

  integer, save :: prtcl ! Keep trace of the protocol used among all the module  
  character(len=120) :: error_string !* Used to return runtime diagnostics
  
  !> Creates socket instances
  interface create
    module procedure socketCreate
  end interface create
  
  !> Send data
  interface send
    module procedure socketSend
  end interface send
  
  !> Retrieve data
  interface retrieve
    module procedure socketRetrieve
  end interface retrieve
  
  public :: socketData, IPI_PROTOCOL1, IPI_PROTOCOL2
  public :: create, send, retrieve
  
contains
  
  !> Connect to socket comms.
  !! \param self Instance.
  !! \param natom number of atoms initially, used as a sanity check later.
  !! \param host host name.
  !! \param verbosity level of detail on communications.
  !! \param protocol type of message headers and data to use,
  !!        currently only IPI_PROTOCOL1 understood.
  !! \param port port to connect to if using an internet protocol, if
  !!        absent, its a file system connection
  subroutine socketCreate(self, nAtom, host, verbosity, protocol, port)
    type(socketData), pointer       :: self
    integer, intent(in)             :: nAtom
    integer, intent(in)             :: verbosity
    character(lc), intent(in)       :: host
    integer, intent(in)             :: protocol
    integer, intent(in), optional   :: port
    
    logical :: tUnix

    ASSERT(.not.associated(self))

    INITALLOCATE_P(self)

    ASSERT(.not.self%tInit)
    
    ASSERT(natom > 0)
    ASSERT(verbosity >= 0)
    
    self%nAtom = nAtom    
    self%verbosity = verbosity
    
    if ((protocol /= IPI_PROTOCOL1) .and. (protocol /= IPI_PROTOCOL2)) then
      write(*,*)protocol, IPI_PROTOCOL1
      call error("Unknown message protocol")
    end if

    prtcl = protocol
    
    if( verbosity > 0 )then
      write(*,*) "socketCreate: Opening socket for two-way communication&
          & with a server."
    end if
    
    tUnix = .true.
    if (present(port)) then ! internet socket to be opened
      ASSERT(port > 0)
      tUnix = .false.
      if( verbosity > 0 )then
        write(*,*) 'Establishing an internet connection to'
        write(*,*) 'Host: ', trim(host)
        write(*,*) 'Port: ', port
      end if
    else
      if( verbosity > 0 )then
        write(*,*)'Establishing UNIX socket connection to', trim(host)
      end if
    end if
    
    if( tUnix ) then ! local socket comms can be used
      call connect_socket(self%socket, host)
    else
      call connect_socket(self%socket, host, port)
    end if
    
    if( self%Verbosity > 0 )then
      write(*,*) "socketCreate: ...Done"
    end if
    
    self%tInit = .true.
    
  end subroutine socketCreate
  
  !> Retrieve data from an external program via a socket
  !! \param self Instance.
  !! \param coord Atomic coordinates.
  !! \param cell lattice vectors.
  !! \note all data in atomic units, and currently assumes the number
  !! of atoms is fixed according to the socketCreate call.
  subroutine socketRetrieve(self, coord, cell, rid)
    type(socketData), pointer    :: self
    integer , intent(out)        :: rid ! The index of the request, i.e. the replica that the force calculation is for (integer corresponding to the bead index)
    real(dp), intent(out)        :: coord(:,:) ! Atomic coordinates
    real(dp), intent(out)        :: cell(3,3)  ! Lattice vectors
    
    character(len=IPI_MSGLEN)    :: header
    integer :: nAtom
    
    ! single precision in the communication interface for some reason
    real(rdp), allocatable :: commsBuffer1(:)
    real(rdp)              :: commsBuffer2(9)
    integer :: commsBufferI
    character(len=:), allocatable :: commsBufferC

    logical :: isinit = .false.
    ASSERT(associated(self))
    
    ASSERT(self%tInit)
    ASSERT(size(coord,dim=1) == 3)
    
    nAtom = size(coord,dim=2)
    
    if (nAtom /= self%nAtom) then
      write (error_string, "(1X,A,2I4)") "Mismatch in number of atoms in&
          & socketRetrieve", nAtom,self%nAtom
      call error(error_string)
    end if
    
    ALLOCATE_(commsBuffer1, (self%nAtom*3))
    
    if( self%Verbosity > 0 )then
      write(*,*) "socketRetrieve: Retrieving data from socket... "
    end if
    
    ! wait for anything other than 'STATUS' state from the interface,
    ! returning state 'READY' in the meanwhile
    do while(.true.)
      
      call readbuffer(self%socket, header, IPI_MSGLEN)

      select case (prtcl)
        
        case ( IPI_PROTOCOL1 )
          select case (trim(header))
          case ('STATUS')
            if( .not. isinit ) then
              call writebuffer(self%socket,'NEEDINIT    ',IPI_MSGLEN)
              if( self%Verbosity > 2 ) write(*,*) "SocketRetrieve: write to socket: NEEDINIT"
            else
              call writebuffer(self%socket,'READY       ',IPI_MSGLEN)
              if( self%Verbosity > 2 ) write(*,*) "SocketRetrieve: write to socket: READY"
            end if
            
          case ('INIT')
            call readbuffer(self%socket, rid)
            print*, 'HERE'
            if( self%Verbosity > 2 ) write(*,*) "SocketRetrieve: read from socket: rid"
            if( self%Verbosity > 3 ) write(*,*) rid
            call readbuffer(self%socket, commsBufferI)
            if( self%Verbosity > 3 ) write(*,*) 'Init msg length:'
            if( self%Verbosity > 3 ) write(*,*) commsBufferI
            if( .not. allocated(commsBufferC) ) allocate(character(len=commsBufferI) :: commsBufferC)
            call readbuffer(self%socket, commsBufferC, commsBufferI)
            if( self%Verbosity > 3 ) write(*,*) 'Init msg:'
            if( self%Verbosity > 3 ) write(*,*) commsBufferC
            isinit = .true.
            
          case default
            exit
          end select

        case( IPI_PROTOCOL2 )

          rid = -1
          
          if( self%Verbosity > 2 ) then
            write(*,*) "SocketRetrieve: read from socket: ",trim(header)
          end if
          
          if (trim(header)/='STATUS') then
            exit
          end if
          
          call writebuffer(self%socket,"NEEDINIT       ",IPI_MSGLEN)
          
          if( self%Verbosity > 2 ) then
            write(*,*) "SocketRetrieve: write to socket: READY"
          end if

        case default
          exit
        end select
        
    enddo
    
    if( self%Verbosity > 2 ) then
      write(*,*) "SocketRetrieve: read from socket: ", trim(header)
    end if
    
    ! expecting positions data
    if (trim(header) /= "POSDATA") then
      call error('coordsFromSocket: Unexpected message from server')
    end if
    
    ! lattice vector data
    call readbuffer(self%socket, commsBuffer2, 9)
    cell = RESHAPE(commsBuffer2,(/3,3/))
    
    if( self%Verbosity > 2 ) then
      write(*,*) "SocketRetrieve: read from socket: cell"
      if( self%Verbosity > 3 ) then        
        write(*,'((3f12.6))') cell
      end if
    end if
    
    ! inverse lattice vectors
    call readbuffer(self%socket, commsBuffer2, 9)
    
    if( self%Verbosity > 2 )  then
      write(*,*) "SocketRetrieve: read from socket: inverse of cell"
    end if
    
    ! number of atomic coordinates
    call readbuffer(self%socket, natom)    
    
    if( self%Verbosity > 2 ) then
      write(*,*) "SocketRetrieve: read from socket: number of atoms"
      if( self%Verbosity > 3 ) then
        write(*,*) natom
      end if
    end if
    
    if (nAtom /= self%nAtom) then
      write (error_string, "(1X,A,2I4)")"Mismatch in number of atoms received",&
          & nAtom,self%nAtom
      call error(error_string)
    end if
    
    ! read actual coordinates
    call readbuffer(self%socket, commsBuffer1, natom*3)
    coord = reshape(commsBuffer1, (/ 3 , natom /) )
    
    if( self%Verbosity > 2 ) then
      write(*,*) "SocketRetrieve: read from socket: atomic positions"    
      if( self%Verbosity > 4 ) then
        write(*,'((3f12.6))') coord
      end if
    end if
    
    if( self%Verbosity > 0 )then
      write(*,*) "SocketRetrieve: ...Done"
    end if
    
    DEALLOCATE_(commsBuffer1)
    
  end subroutine socketRetrieve
  
  !> Send data to an external program via a socket
  !! \param self Instance.
  !! \param energy total energy.
  !! \param forces total forces.
  !! \param stress cell stresses.
  !! \note all data in atomic units, and currently assumes the number
  !! of atoms is fixed according to the socketCreate call.
  subroutine socketSend(self, energy, forces, stress)
    type(socketData), pointer :: self
    real(dp), intent(in) :: energy            ! Total Energy
    real(dp), intent(in) :: forces(:,:)   ! Atomic Forces
    real(dp), intent(in) :: stress(3,3)       ! Stress Tensor
    
    character(len=IPI_MSGLEN) :: header
    integer :: nAtom
    
    ASSERT(associated(self))
    
    ASSERT(self%tInit)
    ASSERT(size(forces,dim=1) == 3)
    
    nAtom = size(forces,dim=2)
    
    if (nAtom /= self%nAtom) then
      write (error_string, "(1X,A,2I4)") "Mismatch in number of atoms in&
          & socketSend", nAtom,self%nAtom
      call error(error_string)
    end if
    
    if( self%Verbosity > 0 )then
      write(*,*) "SocketSend: Sending data to socket... "
    end if
    
    ! wait for anything other than a 'STATUS' state from the interface,
    ! returning state 'HAVEDATA' in the meanwhile
    listen: do while (.true.)
      
      call readbuffer(self%socket, header, IPI_MSGLEN)
      
      if( self%Verbosity > 2 ) then
        write(*,*) "SocketSend: read from socket: ", trim(header)
      end if
      
      if (trim(header)/='STATUS') then
        exit listen
      end if
      
      ! announce that we have available data
      call writebuffer(self%socket,"HAVEDATA    ",IPI_MSGLEN)
      
      if( self%Verbosity > 2 ) then
        write(*,*) "SocketSend: write to socket: HAVEDATA"
      end if
      
    enddo listen
    
    if( self%Verbosity > 2 ) then
      write(*,*) "SocketSend: read from socket: ", trim(header)
    end if
    
    ! expecting to send force data
    if (trim(header)/='GETFORCE') then
      call error('forcesToSocket: Error in socket communication!') 
    end if
    
    call writebuffer(self%socket, "FORCEREADY  ",IPI_MSGLEN)
    if( self%Verbosity > 2 ) then
      write(*,*) "SocketSend: write to socket: FORCEREADY"
    end if
    
    ! transmit total energy
    call writebuffer(self%socket, energy)
    if( self%Verbosity > 2 ) then
      write(*,*) "SocketSend: write to socket: energy"
      if( self%Verbosity > 3 ) then
        write(*,*) energy
      end if
    end if
    
    ! transmit number of atoms we have
    call writebuffer(self%socket, natom)
    if( self%Verbosity > 2 ) then
      write(*,*) "SocketSend: write to socket: natom"
    end if
    
    ! transmit forces
    call writebuffer(self%socket, RESHAPE(forces, (/ 3*natom /)), 3*natom)
    if( self%Verbosity > 2 ) then
      write(*,*) "SocketSend: write to socket: forces"
      if( self%Verbosity > 4 ) then 
        write(*,'((3f12.6))') forces
      end if
    end if
    
    ! transmit stress
    call writebuffer(self%socket, RESHAPE(stress, (/9/)) ,9)
    
    if( self%Verbosity > 2 ) then
      write(*,*) "SocketSend: write to socket: stress"
      if( self%Verbosity > 3 ) then
        write(*,'((3f12.6))') stress
      end if
    end if
    
    ! i-pi can also receive an arbitrary string, that will be printed
    ! out to the "extra" trajectory file. this is useful if you want
    ! to return additional information, e.g.  atomic charges, wannier
    ! centres, etc. one must return the number of characters, then the
    ! string. here we just send back zero characters.
    call writebuffer(self%socket,0)
    
    if( self%Verbosity > 2 ) then
      write(*,*) "SocketSend: 0: nothing else to send"
    end if
    
    if( self%Verbosity > 0 )then
      write(*,*) "SocketSend: ...Done"
    end if
    
  end subroutine socketSend
  
end module iosockets
