! reading photon sources:

implicit none

character(len=len('1s_x')) :: x, y, z

integer :: i

integer, parameter :: nsources = 2

real, parameter :: 1s_x = 0
real, parameter :: 1s_y = 0
real, parameter :: 1s_z = 0

real, parameter :: 2s_x = 0.5
real, parameter :: 2s_y = 0
real, parameter :: 2s_z = 0

integer, parameter :: nsources = 3

real, dimension (:), allocatable :: xsource

real, dimension (:), allocatable :: ysource

real, dimension (:), allocatable :: zsource



do i= 1, nsources
   write(x, '(I4.4, "s_x")') i
   write(y, '(I4.4, "s_y")') i
   write(z, '(I4.4, "s_z")') i
   xsource(i) = x
   ysource(i) = y
   zsource(i) = z
end do 

end 
