!!******************************************************************************
!!
!!  This file is part of the GODUNOV source code, a program to perform
!!  Newtonian or relativistic magnetohydrodynamical simulations.
!!
!!  Copyright (C) 2007-2013 Grzegorz Kowal <grzegorz@amuncode.org>
!!
!!  This program is free software: you can redistribute it and/or modify
!!  it under the terms of the GNU General Public License as published by
!!  the Free Software Foundation, either version 3 of the License, or
!!  (at your option) any later version.
!!
!!  This program is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!  GNU General Public License for more details.
!!
!!  You should have received a copy of the GNU General Public License
!!  along with this program.  If not, see <http://www.gnu.org/licenses/>.
!!
!!******************************************************************************
!!
!! module: MCRT
!!
!!  This module provides subroutines and variables related to radiation tranfer
!!  of photons from one or more sources and computes the resulting ionization 
!!  of the cells and currently uses a two temperature scheme (one temperature for 
!! ionized material - Ti - and one for the neutral material - Tn).
!! 
!!
!!******************************************************************************

module mcrt

implicit none

!=== MODULE VARIABLES ===
!
! nphotons = total number of photons (all sources)
! nsource = number of sources
! nphot = photons from each source
! length = the length unit (given in cm) used for the simulation, currently set to 1 pc.
! xmax, ymax and zmax = half the size of the x, y and z axis. For instance if xmax = 5 then my simulation box is 10 parsecs long in the x direction.
! tema = temperature of the neutral gas in Kelvin
! temb = temperature of the ionized gas in Kelvin
! sigd = x-sect for diffuse photons in $cm^2$ set to 6.3d-18
! alpha\_a= in $cm^3$/s recombination coefficient (for T=8000K using 5.25e-13 and for T=10000K using  4.18e-13)
! vcell = volume of the cell
! nfrac = fraction of neutral hydrogen
! rhokap = opacity 
! ntot = total (ionized + neutral) number density of hydrogen 
! jmean = sum of the pathlengths travelled by all photons in cells they've crossed during one interaction
! temp = temperature
! pres = pressure
! xface, yface, zface = faces of the cells in each axis
!

!variables that we will no longer need when properly merged with GC

real*8 xp,yp,zp,nxp,nyp,nzp,sint,cost,sinp,cosp,phi
real*8 xface(nxg+1),yface(nyg+1),zface(nzg+1) !y
real*8 rhokap(nxg,nyg,nzg) !y
real*8 jmean(nxg,nyg,nzg) !y
real*8 jmean_tot(nxg,nyg,nzg)
real*8 nfrac(nxg,nyg,nzg) !y
real*8 ntot(nxg,nyg,nzg)  !y
real*8 tem(nxg,nyg,nzg) 
real*8 pres(nxg,nyg,nzg)
integer nphotons, niter,nph,iseed,is,jcount,nscatt,totscatt,xl,yl
integer xcell,ycell,zcell,tflag,aflag,i,j,k,mubin,iter
integer pesc,nphot
integer nsource   !y
integer nbub,d_absorb
integer Qphot, remain, tnc
integer ierr, irank, isize,l
real*8 start,finish
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



!********************************************************************************
!----This section will not need to be here as grid already created in Godunov Code


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

!**********  Linear Cartesian grid. Set up grid faces ****************

      do i=1,nxg+1
         xface(i)=(i-1)*2.*xmax/nxg
      end do
      do i=1,nyg+1
         yface(i)=(i-1)*2.*ymax/nyg
      end do
      do i=1,nzg+1
         zface(i)=(i-1)*2.*zmax/nzg
      end do

      rhosum=0.
!**************  Loop through x, y, and z to set up grid density.  
      do i=1,nxg
       do j=1,nyg
        do k=1,nzg
           x=xface(i)-xmax+xmax/nxg
           y=yface(j)-ymax+ymax/nyg
           z=zface(k)-zmax+zmax/nzg
           nfrac(i,j,k)=1.e-6
           rhokap(i,j,k)=ntot(i,j,k)*nfrac(i,j,k)*sigr0
           rhosum=rhosum+ntot(i,j,k)
        end do
       end do
      end do
































