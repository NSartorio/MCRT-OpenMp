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

!parameters in the module (not to alter)

integer, dimension(:), allocatable :: seed_rank 

integer :: remain

integer :: xmax = 1, ymax = 1, zmax = 1

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

!set up uniform density grid

print*, "density grid", dens(:,:,:)

   do i = 1,in
    do j = 1, jn
     do k = 1, kn
        dens(i,j,k) = dens(i,j,k)*unitdens
     end do
    end do
  end do

print*, "density grid", dens

!fing the cross section to be used for length unit in question:
      sigr0=sigd*length !cross section*length unit (cm^3)

!set up neutral fractions and opacities on the grid
print*, "sigd", sigd
print*, "length", length
print*, "sigr", sigr0

      do i=1,in
       do j=1,jn
        do k=1,kn
!set neutral fraction very low (quicker convergence)
           nfrac(i,j,k)=1.e-6
!opacity = density*neutral fraction*cross section
           rhokap(i,j,k)=dens(i,j,k)*nfrac(i,j,k)*sigr0
        end do
       end do
      end do

print*, "opacity grid", rhokap

  stop
end
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
