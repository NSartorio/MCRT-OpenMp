!! Mcion module
!! 
!! This module computes the effect of radiation  
!! from one or multiple stars. The radiation is
!! transfer uses a Monte Carlo method and currently
!! uses a two temperature scheme (one temperature for 
!! ionized material - Ti - and one for the neutral 
!! material - Tn.
!!

module mcion

implicit none

use omp_lib

!include this file with parameters that will be set at ".in" file

     include 'blast.txt'

!*****Constants**************************************************************

!seed for ran2 random number generator
  integer, parameter :: iseed=7686
  iseed=-abs(iseed)  ! must be negative for ran2

!probability of emission of an ionizing photon after absorption by an H atom
  integer, parameter :: pionize=0.38

!

!***** Parameter declarations ***********************************************

real*8 xp,yp,zp,nxp,nyp,nzp,sint,cost,sinp,cosp,phi
real*8 fi,fq,fu,fv
real*8 xface(nxg+1),yface(nyg+1),zface(nzg+1)
real*8 rhokap(nxg,nyg,nzg)
real*8 jmean(nxg,nyg,nzg)
real*8 jmean_tot(nxg,nyg,nzg)
real*8 nfrac(nxg,nyg,nzg)
real*8 ntot(nxg,nyg,nzg)
real*8 tem(nxg,nyg,nzg)
real*8 pres(nxg,nyg,nzg)
      integer nphotons, niter,nph,iseed,is,jcount,nscatt,totscatt,xl,yl
      integer xcell,ycell,zcell,tflag,aflag,i,j,k,mubin,iter
      integer pesc,nphot
      integer nsource   !&&&&
      integer nbub,d_absorb
      integer Qphot, remain, tnc
      integer ierr, irank, isize,l
      real*8 start,finish
      integer, dimension(:), allocatable, save :: seed_rank 
      !$OMP THREADPRIVATE(seed_rank)
      real*8 kappa,albedo,hgg,pl,pc,sc,xmax,ymax,zmax,rimage,tema,temb
      real*8 pi,twopi,fourpi,g2,ximage,yimage,ph,wgt1,wgt,tau,tau1,imtot
      real*8 sigma,jfac,Q49,sigr0,sigd,alpha_a,length,vcell
      real*8 zstars,sigdustR,pionize,pdust, Q490, time
      integer nph_tot
      real*8 unitdens, unitmassdens,unitElow,UnitEhigh
      real*8 conversion, factor, rhosum
      real*8 xsource(1),ysource(1),zsource(1)   !&&&&
      real*8 lsource(1),rbub(1),lumtot   !&&&&
      real  ran2

!*****************OPEN MP*********************************

integer, dimension(:), allocatable :: seed_rank 
integer, save :: iseed
!$OMP THREADPRIVATE(iseed)
integer, save :: irank
!$OMP THREADPRIVATE(irank)
integer :: tseed = 1234
tseed = -abs(tseed)


!**************** BASIC PARAMETERS ***********************

!size of the box in x, y and z axis
      xmax=5.
      ymax=5.
      zmax=5.
!unit mass density
      unitmassdens=unitdens*1.67d-24 !mass of H in g
!cross section for diffuse photons
      sigr0=sigd*length !cross section*length unit (cm^3)
!2*pi  
      twopi=8.*atan(1.)
!cell volume
      vcell = (2.*xmax/nxg)*(2.*ymax/nyg)*(2.*zmax/nzg)   
!factor that multiplies the intensity for ionization fraction calc.
      factor = (1.d49*Q49*sigd)/(vcell*alpha_a*length*length)
!total number of cells
      tnc = nxg*nyg*nzg

!***************start openmp

!$omp parallel 
isize = omp_get_num_procs ( )   !finds n of cores
irank = omp_get_thread_num ( )  !finds which processor I am
!$omp end parallel

!******************	HEADER ******************************

!define different random numbers for each core
!create a seed array
allocate(seed_rank(isize))
  do i = 1, isize
   seed_rank(i) = int(-ran2(tseed)*2.**31)
  end do

!print*, "the seeds are:", seed_rank
!$omp parallel
    iseed = seed_rank(irank + 1)
!$omp end parallel
deallocate(seed_rank)

!$omp parallel 
!photons divided between threads:

!send photons to different cores  
    remain = mod(nphotons,isize)         !finding remaining photons
    nphotons = int(nphotons/isize)     !dividing photons in the # of cores
!distributing remaining photons among cores

    if (irank.lt.remain) then    !if rank is less than the remainder add a further
        nphotons = nphotons + 1        !photon to that core 
    end if

!    print*, "i'm processor", irank, "with nphotons", nphotons
!$omp end parallel 

!**************  Loop through x, y, and z to set up grid density
    do i=1,nxg
      do j=1,nyg
        do k=1,nzg
           open(unit=7, file="density1.dat")
           read(7,*), ntot(i,j,k)
!           ntot(i,j,k) = ntot(i,j,k)*unitdens
!           rhokap(i,j,k)=0.
!           tem(i,j,k) = 0.
        end do
      end do
    end do
    close(7)

!_--_--_Set up density grid_--_--_--_--_--_--_--_--_--_--_--_--_--
      call gridset(xface,yface,zface,rhokap,rhosum,sigr0,length,&
                                &xmax,ymax,zmax,nfrac,ntot)

!_--_--_Set up point source locations, luminosities, total luminosity
      call sources(xsource,ysource,zsource,lsource,lumtot)

      !print*,"lumtot", lumtot

!_--_--_Divide the photons into cores_--_


!--------------------START ITERATIONS--------------------!



do iter = 1, niter     !**********LOOP over iterations
!print*, "iteration", iter
!  if (iter.ge.7) then
!     nphotons = 10*nphotons
!  end if

!set counter to zero
  jcount = 0
    do is=1,nsource !**********LOOP over sources

!calculate number of photons per source
          nph=int(nphotons*lsource(is)/lumtot)
        
          do j=1,nph !**********LOOP over photons from each source

!print @ which point we are
           jcount=jcount+1
!           if((mod(jcount,2).eq.0).and.(irank.eq.0)) then
!             print *, 'iteration ',iter,',',jcount,' photons completed'
!           end if
           
!***** Release photon from point source *******************************
!*************** Cartesian Grid 
          xp = xsource(is)*xmax
          yp = ysource(is)*xmax
          zp = zsource(is)*xmax
!*************** Linear Grid
          xcell=int(nxg*(xp+xmax)/(2.*xmax))+1
          ycell=int(nyg*(yp+ymax)/(2.*ymax))+1
          zcell=int(nzg*(zp+zmax)/(2.*zmax))+1
!emit photon
          call emit(nxp,nyp,nzp,twopi,iseed)

          sigma=0.44 ! cross section for source photons: flux averaged
          nscatt=0   ! it hasn't been scattered yet
          tflag=0    !it is within the grid and it is ionizing

!******* Find first scattering location for random tau.
          call tauint(xp,yp,zp,nxp,nyp,nzp,xmax,ymax,zmax,&
     &xface,yface,zface,rhokap,ntot,jmean,&
     &xcell,ycell,zcell,sigma,tflag,iseed, irank)

!--------------------------start do while---------------------------------!
          do while(tflag.eq.0)

!_______________________________________________
!tflag signalizes when a photon no longer      |
!need to be tracked, that is, if a photon has  |
!exited the grid or it was re-emitted as       |
!a non-ionizing photon                         |
!----------------------------------------------!

!check if it is going to be ionizing or not
              if(ran2(iseed).le.pionize) then

!photon absorbed and re-emitted as ionizing photon: change
!cross section and re-emit isotropically
                 call emit(nxp,nyp,nzp,twopi,iseed)

                 sigma=1.  !cross section for diffuse photons
              else
!********* photon absorbed and re-emitted as a non-ionizing photon:terminate 
                 !print*, "terminate: photon re-emitted as non-ionizing"
                 tflag = 1
                 goto 100 ! terminate
              end if

!************ Find next scattering location
              call tauint(xp,yp,zp,nxp,nyp,nzp,xmax,ymax,zmax,&
     &xface,yface,zface,rhokap,ntot,jmean,&
     &xcell,ycell,zcell,sigma,tflag,iseed, irank)

           end do	
!--------------------------end do while---------------------------------!
100       continue
          end do !**********END LOOP over photons
        end do   !**********END LOOP over sources

        call MPI_Barrier( MPI_COMM_WORLD, ierr)
!reduce - only processor 0 has reduced values
        call MPI_REDUCE(jmean,jmean_tot,tnc,MPI_REAL8,MPI_SUM,&
          0,MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(nph,nph_tot,1,MPI_INT,MPI_SUM,&
          0,MPI_COMM_WORLD,ierr)
!processor 0 calculates nfrac!
        IF (irank .eq. 0) THEN
              !print*,"jmean tottal", jmean_tot(65,65,64)
              !print*,"nph tottal", nph_tot
              jfac = factor/nph_tot
              !print*,"nph tottal", jfac
              !print*,"nfrac", nfrac(65,65,64)
              !print*,"ntot", ntot(65,65,64)
              call ionize(jmean_tot, jmean, rhokap,nfrac,ntot,jfac,sigr0)
              !print*,"jmean tottal", jmean_tot(65,65,64)
             !print*,"jmean tottal", jmean(65,65,64)
              !print*,"nfrac", nfrac(65,65,64)
              !print*, ntot(65,65,64)
        ENDIF
!let other processors know the change in jmean (set to 0 again) and
!the updated rhokap
        call mpi_bcast(jmean, tnc, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(rhokap, tnc, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

!let other processors have access to the nfrac so they
!can calculate the temperature and pressure!
        call mpi_bcast(nfrac, tnc, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

end do !**********END LOOP over iterations


!**************update the opacity (rhokap)************************

!_________________________________________________________________
!use the sum of path lengths to compute ionizing fraction        |
!and adjust opacity accordingly (ie. the fraction of H ionized   |
!does not contribute to the opacity).                            |
!this will be the opacity grid the photons of the next           |
!iteration will see.                                             |
!----------------------------------------------------------------!

!IF (irank .eq. 0) THEN 

!******calculate the temperature and pressure
      tem(:,:,:) = tema
      ntot(:,:,:) = ntot(:,:,:)/unitdens 
      pres(:,:,:) = ntot(:,:,:)*(tema/tema)         
       do i = 1,nxg	
         do j = 1,nyg
            do k = 1,nzg
!            if (nfrac(i,j,k) .lt. 0.5) then 
      tem(i,j,k) = tema+temb*(1.-nfrac(i,j,k))
      pres(i,j,k) = (2.-nfrac(i,j,k))*pres(i,j,k)*(tem(i,j,k)/tema)*0.01
      if (pres(i,j,k).lt.0) then
            endif
             end do
          end do
       end do
!IF (irank .eq. 0) THEN 
!print*, "temperature", tem(65,65,64)
!print*, "nfrac", nfrac(65,65,64)
!print*, "pressure", pres(65,65,64)
!endif

!output ntot, temperature and ionized fraction -------------------
 !     open(unit=10,file='ntot@50.dat')
 !             do i=1,nxg 
!	       write(10,*)(ntot(i,k,50),k=1,nyg)
 !             end do
 !     close(10)
!      open(unit=10,file='tem@50.dat')
!           do i=1,nxg 
!	       write(10,*)(tem(i,k,50),k=1,nyg)
!           end do
!      close(10)


!      open(unit=10,file='nfrac@50.dat')
!              do i=1,nxg 
!	       write(10,*)(nfrac(i,k,50),k=1,nyg)
!              end do
!      close(10)
!ENDIF
return
  

end subroutine


