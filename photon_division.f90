!Testing photon

use omp_lib

!$OMP DO
do i = 1, 1000
...
enddo
!$OMP END D
