!according to Kowal's code
!Add:
!    use parameters, only : get_parameter_integer, get_parameter_real  
!    call get_parameter_integer("nn"  , nn   )
!    in = nn
!    jn = nn
!    kn = nn
!    call get_parameter_integer("ng"  , ng   )
!    call get_parameter_real   ("xmin", xmin )
!    call get_parameter_real   ("xmax", xmax )
!    call get_parameter_real   ("ymin", ymin )
!   call get_parameter_real   ("ymax", ymax )

use omp_lib

!######### from blast file #########
!density (1 code unit is how much?)
integer, parameter ::      unitdens=3.1d3

!unit length (1 code unit is how much?) 1pc (in cm)
real*8, parameter ::      length=3.086e18

!cross section for diffuse photons
real*8, parameter ::      sigd=6.3d-18 ! in cm^2 

real*8 ::  sigr0

integer :: in = 3, jn =3, kn = 4

!photon source location:

real, parameter :: source_x = 0
real, parameter :: source_y = 0
real, parameter :: source_z = 0

!parameters in the module (not to alter)

integer, dimension(:), allocatable :: seed_rank 

integer :: remain

real*8, save :: nxp,nyp,nzp
!$OMP THREADPRIVATE(nxp,nyp,nzp)

integer :: xmax = 0.5, ymax = 0.5, zmax = 0.5

real, dimension (:), allocatable :: xface

real, dimension (:), allocatable :: yface

real, dimension (:), allocatable :: zface

real, dimension (:, :, :), allocatable :: rhokap

real, dimension (:, :, :), allocatable :: nfrac

real, dimension (:, :, :), allocatable :: dens

integer, save :: iseed
!$OMP THREADPRIVATE(iseed)

integer, save :: irank
!$OMP THREADPRIVATE(irank)

integer, save :: nphotons = 1011
!$OMP THREADPRIVATE(nphotons)

integer :: tseed = 8888
tseed = -abs(tseed)


!allocate arrays*********************
allocate(rhokap(in,jn,kn))
rhokap = 0
allocate(dens(in,jn,kn))
dens = 1
allocate(nfrac(in,jn,kn))
allocate(xface(in + 1))
allocate(yface(in + 1))
allocate(zface(in + 1))


!***************start openmp
!$omp parallel 
isize = omp_get_num_procs ( )   !finds n of cores
irank = omp_get_thread_num ( )  !finds which processor I am
!$omp end parallel

print*, -ran2(tseed), tseed

print*, "number of cells", in*jn*kn
print*, rhokap(:,:,:)

!define different random numbers for each core
!create a seed array
allocate(seed_rank(isize))
  do i = 1, isize
   seed_rank(i) = int(-ran2(tseed)*2.**31)
  end do

print*, "the seeds are:", seed_rank
!$omp parallel
!send array as a different seed for each array!
    iseed = seed_rank(irank + 1)
   
print*, "i'm processor", irank, "with iseed", iseed
!$omp end parallel
deallocate(seed_rank)

!Do we need thi bit? won't openmp automatically deal with this?
!$omp parallel 
!photons divided between threads:

!send photons to different cores  
    remain = mod(nphotons,isize)         !finding remaining photons
    nphotons = int(nphotons/isize)     !dividing photons in the # of cores
!distributing remaining photons among cores

    if (irank.lt.remain) then    !if rank is less than the remainder add a further
        nphotons = nphotons + 1        !photon to that core 
    end if

    print*, "i'm processor", irank, "with nphotons", nphotons
!$omp end parallel 

!Set up grid faces based on cell numeber

      do i=1,in+1
         xface(i)=(i-1)*2.*xmax/in
      end do
      do i=1,jn+1
         yface(i)=(i-1)*2.*ymax/jn
      end do
      do i=1,kn+1
         zface(i)=(i-1)*2.*zmax/kn
      end do

print*, "faces", xface(:)


print*, "my source location is ", source_x, source_y, source_z

     
!emitting photons from origin

!$omp parallel
print*, "i'm processor", irank, "with iseed", iseed
call mcemit(nxp,nyp,nzp,iseed)
!$omp end parallel


  stop
end

!###################################################################################
      FUNCTION ran2(idum)
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      real ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,&
     & IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,&
     & NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=min(AM*iy,RNMX)
      return
      END
!############################################################################################


subroutine mcemit(nxp,nyp,nzp,iseed)

implicit none
!uses the monte carlo method to select a direction
!of travel for the photon and sets the photon
!at the position of the source.

!communicated variables
 real*8 :: nxp,nyp,nzp
!!$OMP THREADPRIVATE(nxp,nyp,nzp)
!local variables
!random number
 real  ran2
!!$OMP THREADPRIVATE(ran2)
!the value of two pi
 real*8 , parameter :: twopi = 8.*atan(1.)
!the seed of each processor
 integer iseed
!!$OMP THREADPRIVATE(iseed)
 real*8 sint,cost,sinp,cosp,phi


!***** emit photon isotropically from point source location
print*, "i'm processor with iseed", iseed

!angle theta (measured from z axis) between 0 and pi
 cost=2.*ran2(iseed)-1        !value btw -1 and 1
 sint=sqrt(1.-cost*cost)      !value btw 0 and 1

!angle phi between 0 and 2pi
 phi=twopi*ran2(iseed)
 cosp=cos(phi)
 sinp=sin(phi)

!***** Set photon direction cosines for direction of travel 
 nxp=sint*cosp  
 nyp=sint*sinp
 nzp=cost

print*, "my photon direction is", nxp, nyp, nzp

 return
end


!############################################################################################



!###############################################################################








