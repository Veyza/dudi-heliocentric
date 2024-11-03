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

! File: data_in.f90
! Description: Contains subroutines for inputting data defining the dynamical
! problem and essential for calculations.

! Author: Anastasiia Ershova
! E-mail: vveyzaa@gmail.com

module data_in

contains

    ! reads comet's ephemeridae (Np rows) from the file fname 
    ! formatted as 7 columns:
    ! moment [days], 3d Cartesian position [AU], 3d Cartesian velocity [AU/days]
    ! Makes Ns sources at each of the Np positions. The sources are located at
    ! the comet radius Rast_AU [AU] from the comet center and uniformly
    ! distributed over a sphere
   subroutine get_sources(fname, Np, Ns, sources, comet, Rast_AU)
      use const
      use define_types
      use help
      use distributions_fun
      implicit none
      integer, intent(in) :: Np, Ns
      real(8), intent(in) :: Rast_AU
      type(source_properties), intent(out) :: sources(Np,Ns)
      type(ephemeris), intent(out) :: comet(Np)
      integer i, ii, iii
      real(8) moment(Np), xyz(Ns,3), tmp, tmpvec(3)
      character(*), intent(in) :: fname

      open(200, file = fname, status = 'old')
      do i = 1, Np, 1
         read(200,*) moment(i), comet(i)%coords, comet(i)%Vastvec
      enddo
      close(200)

      ! inputing the coordinates and orientations of sources on the surface of the asteroid
      do i = 1, Np, 1
         ! asteroid speed at position i
         comet(i)%Vast = norma3d(comet(i)%Vastvec)
         ! generating Ns points uniformly distributed over a unit sphere
         call uniform_random_sphere(Ns, xyz)
         do ii = 1, Ns
            sources(i,ii)%symmetry_axis = xyz(Ns+1-ii,:)
            sources(i,ii)%rrM = comet(i)%coords &
                + sources(i,ii)%symmetry_axis * Rast_AU
            sources(i,ii)%r = norma3d(sources(i,ii)%rrM)
            sources(i,ii)%alphaM = acos(sources(i,ii)%rrM(3) / sources(i,ii)%r)
            sources(i,ii)%betaM = atan(sources(i,ii)%rrM(2), &
               sources(i,ii)%rrM(1))
            tmp = dot_product(sources(i,ii)%symmetry_axis, &
               sources(i,ii)%rrM/sources(i,ii)%r)
            if(tmp > 1d0) then
               tmp = 1d0
            endif
            if(tmp < -1d0) then
               tmp = -1d0
            endif
            sources(i,ii)%zeta = acos(tmp)
            if(sources(i,ii)%zeta > 1d-8) then
               sources(i,ii)%eta = &
                         azimuth(sources(i,ii)%symmetry_axis, sources(i,ii)%rrM)
            else
               sources(i,ii)%eta = 0d0
            endif
            sources(i,ii)%ud%ud_shape = 0
            sources(i,ii)%ud%umin = 5d0 / AUdays2SI
            sources(i,ii)%ud%umax = 100d0 / AUdays2SI
            sources(i,ii)%ejection_angle_distr = 2
            sources(i,ii)%Tj = moment(i)
            sources(i,ii)%dtau = 1d-2 / s_in_day
            sources(i,ii)%Nparticles = 1d5 + 1d5 &
             * dot_product(comet(i)%Vastvec, sources(i,ii)%symmetry_axis) &
               / comet(i)%Vast
         enddo
      enddo

   end subroutine get_sources



   ! returns N points
   ! uniformly distributed over a sphere
   subroutine uniform_random_sphere(N, xyz)
      use const
      use help
      implicit none
      integer, intent(in) :: N
      real(8), intent(out) :: xyz(N,3)
      integer i
      real(8) polar_angle(N), longitude(N)

      do i = 1, N
         polar_angle(i) = rand(0)
         polar_angle(i) = acos(-2d0 * polar_angle(i) + 1d0)
         longitude(i) = rand(0)
         longitude(i) = longitude(i) * twopi
      enddo

      do i = 1, N
         xyz(i,3) = cos(polar_angle(i))
         xyz(i,1) = sin(polar_angle(i)) * cos(longitude(i))
         xyz(i,2) = sin(polar_angle(i)) * sin(longitude(i))
      enddo

   end subroutine uniform_random_sphere


   ! obtains Cartesian coordinates of a unit vector
   ! in the moon-centered coordinate system
   ! aligned with the ejection symmetry axis
   ! the source's position (rrM - Cartesian coordinates, betaM - eastern longitude),
   ! the jet's zenith angle (zeta) and azimuth (eta) are known
   subroutine jet_direction(betaM, zeta, eta, rrM, jetdir)
      use const
      use help
      real(8), intent(in) :: betaM, zeta, eta, rrM(3)
      real(8), intent(out) :: jetdir(3)
      real(8) xj(3), yj(3), rtmp(3)
      real(8) xout, yout, zout, tmpang

      rtmp = rrM / norma3d(rrM)
      tmpang = 3d0 * halfpi - betaM
      if(zeta /= 0d0) then
         ! rotate the CS around z-axis so that rtmp would have zero x-axis component
         call eulrot(0d0, 0d0, tmpang, rtmp(1), rtmp(2), rtmp(3), &
            xout, yout, zout, .FALSE.)
         rtmp(1) = xout ; rtmp(2) = yout ; rtmp(3) = zout
         ! xj is the vector normal to rtmp
         xj(1) = 0d0
         xj(2) = rtmp(3)
         xj(3) = -rtmp(2)
         xj = xj / norma3d(xj)
         ! yj is the 3rd vector to complete the basis of the horizontal CS at the point rtmp
         yj = vector_product(rtmp,xj)

         ! we know Cartesian coordinates of the desired vector jetdir in the horizontal CS
         ! we know coordinates of the basis vectores of the horizontal CS
         ! in the desired CS
         ! then the coordinates of jetdir in the dsired CS are found
         ! as a sum of the basis vectors with the appropriate coefficients
         jetdir = sin(zeta) * cos(eta) * xj &
            - sin(zeta) * sin(eta) * yj &
            + cos(zeta) * rtmp
         jetdir = jetdir / norma3d(jetdir)

         ! transforming back to the initial CS
         call eulrot(0d0, 0d0, tmpang, rtmp(1), rtmp(2), rtmp(3), &
            xout, yout, zout, .TRUE.)
         rtmp(1) = xout ; rtmp(2) = yout ; rtmp(3) = zout

         call eulrot(0d0, 0d0, tmpang, jetdir(1), jetdir(2), jetdir(3), &
            xout, yout, zout, .TRUE.)
         jetdir(1) = xout ; jetdir(2) = yout ; jetdir(3) = zout
      else
         jetdir = rtmp
      endif

   end subroutine jet_direction



    ! returns a 2d grid of points in the comet's orbital plane
    ! the point with coordinates ´center´ is in the center of this grid
    ! the grid has dimensions of nt1 and nt2
    ! the distance between the grid nodes along the two dimensions
    ! is resolution(1) and resolution(2) respectively
   subroutine orbital_plane_grid(points, nt1, nt2, resolution, comet, center)
      use const
      use define_types
      use help
      implicit none
      integer, intent(in) :: nt1, nt2
      real(8), intent(in) :: resolution(2), center(3)
      type(ephemeris), intent(in) :: comet
      type(position_in_space), intent(out) :: points(nt1,nt2)
      integer i, ii
      real(8) tmpvec(3)
      real(8) zvec(3), xvec(3), yvec(3)

      ! Coordinate system orts
      ! z-axis is normal to the asteroid orbital plane
      zvec = vector_product(comet%Vastvec, comet%coords)
      zvec = zvec / norma3d(zvec)
      ! x-axis points along the asteroid heliocentric radius at the moment tnow
      xvec = comet%coords / norma3d(comet%coords)
      yvec = vector_product(zvec, xvec)

      xvec = xvec * resolution(1) / AU
      yvec = yvec * resolution(2) / AU
      tmpvec = center - nt1 * xvec * 0.5d0 - nt2 * yvec * 0.5d0

      do ii = 1, nt2
         do i = 1, nt1
            points(i,ii)%rvector = tmpvec + (i * xvec) + (ii * yvec)
            points(i,ii)%rvector = points(i,ii)%rvector
            points(i,ii)%r = norma3d(points(i,ii)%rvector)
            points(i,ii)%alpha = acos(points(i,ii)%rvector(3) / points(i,ii)%r)
            points(i,ii)%beta = atan(points(i,ii)%rvector(2), points(i,ii)%rvector(1))    
         enddo
      enddo

   end subroutine orbital_plane_grid


    ! returns the reduced gravitational parameter \mu_R
    ! calculated from the grain radius Rg [m] and the radiation pressure
    ! efficiecy Qpr
   subroutine reduced_gravitational_parameter(Rg, Qpr, muR)
      use const
      implicit none
      real(8), intent(in) :: Rg, Qpr
      real(8), intent(out) :: muR
      real(8) radiation

      radiation = 3d0 / 16d0 / pi * Lsun * Qpr / Rg / rho/ speed_of_light
      muR = gravity_constant * Msun - radiation

      muR = muR / AU**3 * s_in_day**2

   end subroutine reduced_gravitational_parameter



end module data_in
