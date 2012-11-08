MODULE utils
  USE precision
  IMPLICIT NONE
  
CONTAINS
  real(kr) FUNCTION dist(r1,r2)
    IMPLICIT NONE
    real(kr),intent(IN),dimension(3) :: r1,r2
    integer(ki) :: i

    dist=0.0d0
    do i = 1,3
       dist = dist + ( r1(i) - r2(i) )**2
    end do

    dist = sqrt(dist)
    
  END FUNCTION dist

  ! If needed it could be created an interface
  ! in order to compute with, nominally, the 
  ! same function factorial for real number.
  integer(ki) FUNCTION fact(n)
    IMPLICIT NONE
    integer(ki),intent(IN) :: n
    integer(ki) :: i

    fact = 1
    do i = 1,n
       fact = fact * i
    end do

  END FUNCTION fact
  
END MODULE utils

!!$PROGRAM test
!!$  USE utils
!!$  IMPLICIT NONE
!!$  integer :: n
!!$  real(8),dimension(3) :: r1,r2
!!$  
!!$  n=6
!!$  print*, fact(n)
!!$  ! result: 720
!!$
!!$  r1 = (/ 0,0,0 /)
!!$  r2 = (/ 1.,2.5,1.5 /)
!!$  print*, dist(r1,r2)
!!$  ! result: 3.0822
!!$
!!$END PROGRAM test
