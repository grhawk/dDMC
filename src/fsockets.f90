!F90 ISO_C_BINGING wrapper for socket communication.

!Copyright (C) 2013, Michele Ceriotti

!Permission is hereby granted, free of charge, to any person obtaining
!a copy of this software and associated documentation files (the
!"Software"), to deal in the Software without restriction, including
!without limitation the rights to use, copy, modify, merge, publish,
!distribute, sublicense, and/or sell copies of the Software, and to
!permit persons to whom the Software is furnished to do so, subject to
!the following conditions:

!The above copyright notice and this permission notice shall be included
!in all copies or substantial portions of the Software.

!THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
!EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
!MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
!IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
!CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
!TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
!SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


!Contains both the functions that transmit data to the socket and read the data
!back out again once finished, and the function which opens the socket initially.

!Functions:
!   open_socket: Opens a socket with the required host server, socket type and
!      port number.
!   write_buffer: Writes a string to the socket.
!   read_buffer: Reads data from the socket.

MODULE F90SOCKETS
  USE ISO_C_BINDING
  
  IMPLICIT NONE
  
  private
  public :: connect_socket, readbuffer, writebuffer
  
  integer, parameter :: dp = kind(1.0d0) !* precision of the real data type
  
  interface connect_socket
    module procedure connect_isocket
    module procedure connect_nsocket
  end interface connect_socket
  
  INTERFACE writebuffer
    MODULE PROCEDURE writebuffer_s, &
        writebuffer_d, writebuffer_dv, &
        writebuffer_i
    
  END INTERFACE writebuffer
  
  INTERFACE readbuffer
    MODULE PROCEDURE readbuffer_s
    MODULE PROCEDURE readbuffer_dv
    MODULE PROCEDURE readbuffer_d
    MODULE PROCEDURE readbuffer_i
  END INTERFACE readbuffer
  

  ! interface bindings for C routines
  INTERFACE
    
    SUBROUTINE connect_cinetsocket(psockfd, host, port) BIND(C, &
        & name="connect_inet_socket")
      USE ISO_C_BINDING
      INTEGER(KIND=C_INT)                      :: psockfd, port
      CHARACTER(KIND=C_CHAR), DIMENSION(*)     :: host
      
    END SUBROUTINE connect_cinetsocket
    
    SUBROUTINE connect_cnixsocket(psockfd, host) BIND(C, &
        & name="connect_unix_socket")
      USE ISO_C_BINDING
      INTEGER(KIND=C_INT)                      :: psockfd
      CHARACTER(KIND=C_CHAR), DIMENSION(*)     :: host
      
    END SUBROUTINE connect_cnixsocket
    
    SUBROUTINE writebuffer_csocket(psockfd, pdata, plen) BIND(C, &
        & name="writebuffer")
      USE ISO_C_BINDING
      INTEGER(KIND=C_INT)                      :: psockfd
      TYPE(C_PTR), VALUE                       :: pdata
      INTEGER(KIND=C_INT)                      :: plen
      
    END SUBROUTINE writebuffer_csocket
    
    SUBROUTINE readbuffer_csocket(psockfd, pdata, plen) BIND(C, &
        & name="readbuffer")
      USE ISO_C_BINDING
      INTEGER(KIND=C_INT)                      :: psockfd ! The id of the socket that will be read from.
      TYPE(C_PTR), VALUE                       :: pdata   ! The storage array for data read from the socket.
      INTEGER(KIND=C_INT)                      :: plen    ! The length of the data in bytes.
      
    END SUBROUTINE readbuffer_csocket
    
  END INTERFACE
  
CONTAINS
  
  ! internet socket
  SUBROUTINE connect_isocket(psockfd, host, port)
    INTEGER, INTENT(IN) :: port
    INTEGER, INTENT(OUT) :: psockfd
    CHARACTER(LEN=*), INTENT(IN) :: host
    
    INTEGER(KIND=C_INT) :: Cport
    
    Cport = port
    CALL connect_cinetsocket(psockfd, trim(host)//C_NULL_CHAR, Cport)
    
  END SUBROUTINE connect_isocket
  
  ! unix socket
  SUBROUTINE connect_nsocket(psockfd, host)      
    INTEGER, INTENT(OUT) :: psockfd
    CHARACTER(LEN=*), INTENT(IN) :: host
    
    CALL connect_cnixsocket(psockfd, trim(host)//C_NULL_CHAR)
  END SUBROUTINE connect_nsocket
  
  ! write a float (8-byte real)
  SUBROUTINE writebuffer_d (psockfd, fdata)
    USE ISO_C_BINDING
    INTEGER, INTENT(IN)                      :: psockfd
    REAL(dp), INTENT(IN)                     :: fdata
    
    REAL(KIND=C_DOUBLE), TARGET              :: cdata
    
    cdata = fdata
    CALL writebuffer_csocket(psockfd, c_loc(cdata), 8)
  END SUBROUTINE writebuffer_d
  
  ! write an integer (4 byte integer)
  SUBROUTINE writebuffer_i (psockfd, fdata)
    USE ISO_C_BINDING
    INTEGER, INTENT(IN)                      :: psockfd, fdata
    
    INTEGER(KIND=C_INT), TARGET              :: cdata
    
    cdata = fdata
    CALL writebuffer_csocket(psockfd, c_loc(cdata), 4)
  END SUBROUTINE writebuffer_i
  
  ! write a string
  SUBROUTINE writebuffer_s (psockfd, fstring, plen)
    USE ISO_C_BINDING
    INTEGER, INTENT(IN)                      :: psockfd
    CHARACTER(LEN=*), INTENT(IN)             :: fstring
    INTEGER, INTENT(IN)                      :: plen
    
    INTEGER                                  :: ii
    CHARACTER(LEN=1, KIND=C_CHAR), TARGET    :: cstring(plen+1)
    
    DO ii = 1,plen
      cstring(ii) = fstring(ii:ii)
    ENDDO
!    cstring(plen+1) = C_NULL_CHAR
    
    CALL writebuffer_csocket(psockfd, c_loc(cstring(1)), plen)
  END SUBROUTINE writebuffer_s
  
  SUBROUTINE writebuffer_dv(psockfd, fdata, plen)
    USE ISO_C_BINDING  
    INTEGER, INTENT(IN)                     :: psockfd, plen
    REAL(dp), INTENT(IN), TARGET            :: fdata(plen)
    
    REAL(KIND=C_DOUBLE), TARGET              :: cdata(plen)
    
    cdata = fdata
    CALL writebuffer_csocket(psockfd, c_loc(cdata(1)), 8*plen)
  END SUBROUTINE writebuffer_dv
  
  SUBROUTINE readbuffer_d (psockfd, fdata) ! Read (8-bytes) real number from socket
    USE ISO_C_BINDING
    INTEGER, INTENT(IN)                      :: psockfd
    REAL(dp), INTENT(OUT)                    :: fdata
    
    REAL(KIND=C_DOUBLE), TARGET              :: cdata
    
    CALL readbuffer_csocket(psockfd, c_loc(cdata), 8)
    fdata=cdata
  END SUBROUTINE readbuffer_d
  
  SUBROUTINE readbuffer_i (psockfd, fdata) ! Read (4-bytes) integer from socket
    USE ISO_C_BINDING
    INTEGER, INTENT(IN)                      :: psockfd
    INTEGER, INTENT(OUT)                     :: fdata
    
    INTEGER(KIND=C_INT), TARGET              :: cdata
    
    CALL readbuffer_csocket(psockfd, c_loc(cdata), 4)
    fdata = cdata
  END SUBROUTINE readbuffer_i
  
  SUBROUTINE readbuffer_s (psockfd, fstring, plen)
    USE ISO_C_BINDING
    INTEGER, INTENT(IN)                      :: psockfd
    CHARACTER(LEN=*), INTENT(OUT)            :: fstring
    INTEGER, INTENT(IN)                      :: plen
    
    INTEGER                                  :: ii
    CHARACTER(LEN=1, KIND=C_CHAR), TARGET    :: cstring(plen)
    
    CALL readbuffer_csocket(psockfd, c_loc(cstring(1)), plen)
    fstring=""
    DO ii = 1,plen
      if (cstring(ii) == C_NULL_CHAR) then
        exit
      end if
      fstring(ii:ii) = cstring(ii)
    ENDDO
  END SUBROUTINE readbuffer_s
  
  SUBROUTINE readbuffer_dv(psockfd, fdata, plen)
    USE ISO_C_BINDING  
    INTEGER, INTENT(IN)                :: psockfd, plen
    REAL(dp), INTENT(OUT)              :: fdata(plen)
    
    REAL(KIND=C_DOUBLE), TARGET        :: ctmp(plen)
    
    CALL readbuffer_csocket(psockfd, c_loc(ctmp(1)), 8*plen)
    fdata = ctmp
    
  END SUBROUTINE readbuffer_dv
  
END MODULE F90SOCKETS
