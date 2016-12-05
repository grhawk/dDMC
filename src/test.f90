program test
implicit none
integer :: sum, a

integer, parameter :: dp = kind(0.0d0) !* precision of the real data type (8)
integer, parameter :: rdp = kind(0.0d0) !* Real double precision - do not edit (8)
  integer, parameter :: kr4 = selected_real_kind(6,37)     ! (4)
  integer, parameter :: kr8 = selected_real_kind(15,307)   ! (8)
  
integer, dimension(2,3) :: b

print *, "dp = ", dp
print *, "rdp = ", rdp
print *, "kr4 = ", kr4
print *, "kr8 = ", kr8

b = 0.0
print *,b
print *,SIZE(b,dim=1)
print *,SIZE(b,dim=2)
 

end program