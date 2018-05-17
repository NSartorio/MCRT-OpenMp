! reading photon sources:

implicit none

character(len=len('1s_x')) :: x, y, z

integer :: i

integer, parameter :: nsources = 2

real, parameter :: s_x(1) = 0
real, parameter :: s_y(1) = 0
real, parameter :: s_z(1) = 0

real, parameter :: s_x(2) = 0.5
real, parameter :: s_y(2) = 0
real, parameter :: s_z(2) = 0

real, dimension (:), allocatable :: xsource

real, dimension (:), allocatable :: ysource

real, dimension (:), allocatable :: zsource

print*, s_x1, s_x2

do i= 1, nsources
   write(x, '( "s_x", I1.1)') i
   write(y, '( "s_y", I1.1)') i
   write(z, '( "s_z", I1.1)') i
   print*, x
!   xsource(i) = x
!   ysource(i) = y
!   zsource(i) = z
end do 

end 
