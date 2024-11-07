! This file is a part of DUDI-heliocentric, the Fortran-90 implementation 
! of the two-body model for the dynamics of dust ejected from an atmosphereless
! body moving around the Sun
! Version 1.0.0
! This is free software. You can use and redistribute it 
! under the terms of the GNU General Public License (http://www.gnu.org/licenses/)
! If you do, please cite the following paper
! Anastasiia Ershova and Jürgen Schmidt, 
! Two-body model for the spatial distribution of dust ejected from
! an atmosphereless body, 2021, A&A, 650, A186 

! Author: Anastasiia Ershova
! E-mail: vveyzaa@gmail.com

! File: phaethon_input.f90
! Description: This module contains the routines used for inputing
! the dust ejection parameters and the desirable points of interest
! to compute the number density of dust ejected from NEA 3200 Phaethon

module phaethon_input

contains

	! input the paths to the impact-ejecta maps 
	! and the corresponding heliocentric distances
   subroutine get_maps_data(rhels, fnames)
      use const
      implicit none
      integer, parameter :: N = 4
      real(8), intent(out) :: rhels(N)
      character(len = *), intent(out) :: fnames(N)
      integer i
      character(len = 1) ic

      rhels = (/0.96d0, 0.94d0, 0.93d0, 0.91d0/)
      do i = 1, N
         write(ic,'(I1)') i+4
         fnames(i) = './input_data_files/impact_ejecta_maps/COS3_Yield_ALL_map_' &
            // ic // '.sav_yield_orb_pha.txt'
      enddo

   end subroutine get_maps_data


	! set the parameters of the sources modeling dust ejection
	! from the asteroid moving along its orbit
	! fname is the name of the file
	! containing the asteroid ephemeridae (positions and velocities)
   subroutine get_moving_sources(fname, Np, Neph, Nlin, &
                                                  sources, comet)
      use const
      use define_types
      use help
      use distributions_fun
      implicit none
      integer, parameter :: Nmaps = 9
      integer, parameter :: N = 500
      integer, intent(in) :: Np, Neph, Nlin
      ! Phaethon's radius is used as a normalization factor when 
      ! calculating the dust production rate
      real(8), parameter :: Rast = 2.9e3		
      type(source_properties), intent(out) :: sources(Np)
      type(ephemeris), intent(out) :: comet(Np)
      integer i, ii, iii
      real(8) moment(Np), dNlin
      real(8) rates(nlons, nlats), totrate
      character(*), intent(in) :: fname
      character(len = 93), dimension(Nmaps) :: fnames
      real(8) rhels(Nmaps), rhel1, rhel2
      integer mapind1, mapind2

      call get_maps_data(rhels, fnames)
      mapind1 = 1 ; mapind2 = mapind1 + 1
      call read_first_ratemap(fnames(mapind1), rhel1)
      call read_ratemap(fnames(mapind2), rhel2)
       
      dNlin = dble(Nlin)
      open(200, file = fname, status = 'old')
      ! reading the 1st line of the ephemeridae
      read(200,*) moment(1), comet(1)%coords, comet(1)%Vastvec
      do i = Nlin+1, Np, Nlin
		! reading subsequent lines of the ephemeridae
         read(200,*) moment(i), comet(i)%coords, comet(i)%Vastvec
         ! linearly interpolating in between
         if(Nlin > 1) then
            forall(ii = (i-Nlin+1):(i-1))
               moment(ii) = moment(i-Nlin) &
                  + (moment(i) - moment(i-Nlin)) * dble(ii-i+Nlin) / dNlin

               comet(ii)%Vastvec = comet(i-Nlin)%Vastvec &
                  + (comet(i)%Vastvec - comet(i-Nlin)%Vastvec) * dble(ii-i+Nlin) / dNlin

               comet(ii)%coords = comet(i-Nlin)%coords &
                  + (comet(i)%coords - comet(i-Nlin)%coords) * dble(ii-i+Nlin) / dNlin
            endforall
         endif
      enddo
      close(200)
      ! shifting the timeline
      forall(i = 1:Np) moment(i) = moment(i) - moment(1)


      ! inputing the coordinates of the point sources representing the asteroid
      do i = 1, Np
         sources(i)%rrM = comet(i)%coords 
         sources(i)%r = norma3d(comet(i)%coords)
         if(sources(i)%r < rhel2) then
             rhel1 = rhel2
             mapind1 = mapind2
             mapind2 = mapind2 + 1
             rmap1 = rmap2
             call read_ratemap(fnames(mapind2), rhel2)
          endif
          call ratematr_interpolate(sources(i)%r, rhel1, rhel2)
         comet(i)%Vast = norma3d(comet(i)%Vastvec)    ! asteroid speed at position i
         ! generating Ns points uniformly distributed over a unit sphere
         
		sources(i)%alphaM = acos(sources(i)%rrM(3) / sources(i)%r)
		sources(i)%betaM = atan(sources(i)%rrM(2), &
									sources(i)%rrM(1))
		sources(i)%symmetry_axis = sources(i)%rrM(1) / sources(i)%r
		sources(i)%zeta = 0d0
		sources(i)%eta = 0d0
		sources(i)%ud%ud_shape = 1
		sources(i)%ud%umin = 2d0 / AUdays2SI
		sources(i)%ud%umax = 2399d0 / AUdays2SI
		sources(i)%ejection_angle_distr = 3
		sources(i)%Tj = moment(i)
		sources(i)%dtau = 0d0
		
		! integrate the number density of impact ejecta over the matrix
		call integrate_over_matrix(totrate)
		! converting the number density to flux
		! see Eq. 3 from the Szalay et al, 2016 (asteroid on a spherical orbit)
		totrate = totrate / 0.31d0 / 7.2e-3 / 4d0 / pi * Rast**2
		! converting flux to the number of ejected particles
		sources(i)%Nparticles = totrate &
		                    * (moment(2) - moment(1)) * s_in_day
      enddo

   end subroutine get_moving_sources




   subroutine integrate_over_matrix(integral)
      use const
      use distributions_fun
      implicit none
      real(8), intent(out) :: integral
      integer i, ii
      real(8) rint, rint1, polangle(nlats)

      integral = 0d0

      polangle = halfpi - lats
      do ii = 2, nlats
         rint = 0d0
         ! integrating over a ring at iith latitude
         do i = 2, nlons
            rint = rint + (ratemap(i,ii) + ratemap(i-1,ii)) &
               * (lons(i) - lons(i-1))
         enddo
         ! closing the ring
         rint = rint + (ratemap(1,ii) + ratemap(nlons,ii)) * (lons(1) - lons(nlons) + twopi)
         ! trapezoid formula is used
         rint = rint / 2d0

         rint1 = 0d0
         !integrating over a ring at ii-1th latitude
         do i = 2, nlons
            rint1 = rint1 + (ratemap(i,ii-1) + ratemap(i-1,ii-1)) &
               * (lons(i) - lons(i-1))
         enddo
         ! closing the ring
         rint1 = rint1 + (ratemap(1,ii-1) + ratemap(nlons,ii-1)) &
            * (lons(1) - lons(nlons) + twopi)
         rint1 = rint1 / 2d0   ! the factor from the trapezoid formula

         ! trapezoid formula
         integral = integral + (rint + rint1) * sin(polangle(ii)) &
            * (polangle(ii-1) - polangle(ii)) / 2d0
      enddo
      rint = 0d0
      ! integrating over the uppermost ring
      do i = 2, nlons
         rint = rint + (ratemap(i,1) + ratemap(i-1,1)) &
            * (lons(i) - lons(i-1))
      enddo
      rint = rint + (ratemap(1,1) + ratemap(nlons,1)) &
         * (lons(1) - lons(nlons) + twopi)
      rint = rint / 2d0
      integral = integral &
         + rint * sin((-halfpi - lats(1)) / 2d0) &
         * (-halfpi - lats(1))

      ! integrating over the nlatsth ring
      rint = 0d0
      do i = 2, nlons
         rint = rint + (ratemap(i,nlats) + ratemap(i-1,nlats)) &
            * (lons(i) - lons(i-1))
      enddo
      rint = rint + (ratemap(1,nlats) + ratemap(nlons,nlats)) &
         * (lons(1) - lons(nlons) + twopi)
      integral = integral &
         + rint * sin((halfpi - lats(nlats)) / 2d0) &
         * (halfpi - lats(nlats))

   end subroutine integrate_over_matrix


	

	! Inputing a 2-d array of the points where the number density
	! will be computed
   subroutine get_points(points, nt1, nt2, resolution, lastrM, &
                                          VVast, Vast, cntrpx, cntrpy)
      use const
      use define_types
      use help
      implicit none
      integer, intent(in) :: nt1, nt2
      real(8), intent(in) :: resolution(2), cntrpx, cntrpy
      type(position_in_space), intent(out) :: points(nt1,nt2)
      real(8), intent(in) :: lastrM(3), VVast(3), Vast
      integer i, ii
      real(8) heldist
      real(8) tmpvec(3)
      real(8) zvec(3), xvec(3), yvec(3)

      ! the CS to compare with Szalay et al, 2019
       zvec = (/0d0, 0d0, 1d0/)
       xvec = lastrM  / norma3d(lastrM)
       zvec = zvec - zvec * dot_product(xvec, zvec)
       zvec = zvec / norma3d(zvec)
       yvec = vector_product(zvec, xvec)

      xvec = xvec * resolution(1) / AU
      yvec = yvec * resolution(2) / AU
      tmpvec = lastrM - nt1 * xvec * cntrpx - nt2 * yvec * cntrpy

      do ii = 1, nt2
         do i = 1, nt1
            points(i,ii)%rvector = tmpvec + (i * xvec) + (ii * yvec)
            points(i,ii)%rvector = points(i,ii)%rvector
            points(i,ii)%r = norma3d(points(i,ii)%rvector)
            points(i,ii)%alpha = acos(points(i,ii)%rvector(3) / points(i,ii)%r)
            points(i,ii)%beta = atan(points(i,ii)%rvector(2), points(i,ii)%rvector(1))
         enddo
      enddo

   end subroutine get_points

	
	! From the given table of particle radii and corresponding beta
	! values, interpolate the value of beta corresponding to Rg
   subroutine beta_from_Rg(beta, Rg)
      use const
      use help
      implicit none
      real(8), intent(out) :: beta
      real(8), intent(in) :: Rg 
      integer, parameter :: N = 500
      real(8) Rgs(N), betas(N), mass1
      integer i

      open(300, file = "input_data_files/beta_vs_forsterite_Rg.dat", status = "old")
         do i = 1, N
            read(300,*) Rgs(i), betas(i)
         enddo
      close(300)
      beta = LiNTERPOL(N, betas, Rgs, Rg)
   
   end subroutine beta_from_Rg



end module phaethon_input
