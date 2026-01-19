!===============================================================
! Program: phaethon_3d
! Description: Compute dust number density on a 3-D grid and
!              save each Rg result as a raw binary file.
!===============================================================
program phaethon_3d
   use const
   use help
   use distributions_fun
   use define_types
   use DUDIhc
   use phaethon_input
   use data_out           ! for read_* ratemap helpers and interpolation
   use data_in            ! for get_moving_sources, runge_kutta_point_position
   use OMP_LIB
   implicit none

   integer, parameter :: Neph = 2000
   integer, parameter :: Nlin = 9
   integer, parameter :: Np = (Neph-1)*Nlin + 1
   integer, parameter :: Nmaps = 4
   integer, parameter :: Nrgs = 4

   ! Grid size (adjust as needed; memory ~ nx*ny*nz*8 bytes per array)
   integer :: nx, ny, nz
   real(8) :: resolution(3)
   real(8), parameter :: centerpositionx = 0.5d0
   real(8), parameter :: centerpositiony = 0.5d0
   real(8), parameter :: centerpositionz = 0.5d0

   ! Time and dynamics
   real(8) :: tnow, muR, dt, beta
   real(8) :: dtlim2, dtlim3
   integer :: idt

   ! Indices
   integer :: i_R, i, j, k, i_p
   integer :: mapind1, mapind2

   ! Data
   type(source_properties), allocatable :: sources(:)
   type(position_in_space), allocatable :: points(:,:,:)
   type(ephemeris), allocatable :: comet(:)
   real, allocatable :: density(:,:,:), tmp_res(:,:,:)

   real(8) :: cloudcentr(3)
   character(len=66) :: fname
   character(len=64) :: fnameout
   character(len=93) :: fnames(Nmaps)
   real(8) :: rhels(Nmaps), rhel1, rhel2
   real(8) :: Rgs(0:Nrgs)

   !--- Setup ---
   fname = 'input_data_files/Phaethon_2025-02-22_last_int=10min_ECLIPJ2000.dat'

   Rgs = (/0.2d0, 0.3d0, 1d0, 1.2d0, 2.5d0/)

   ! Grid definition (example spacing equal on all axes)
   nx = 100 ; ny = 100 ; nz = 100
   resolution(1) = 16d3        ! meters
   resolution(2) = resolution(1)
   resolution(3) = resolution(1)

   allocate(points(nx,ny,nz), tmp_res(nx,ny,nz), density(nx,ny,nz), &
            sources(Np), comet(Np))

   ! Sources along trajectory and ephemerides
   call get_moving_sources(fname, Np, Nlin, sources, comet)
   tnow = sources(Np)%Tj

   ! Build 3-D grid around the last comet position
   call get_points_3d(points, nx, ny, nz, resolution, comet(Np)%coords, &
                      centerpositionx, centerpositiony, centerpositionz)

   ! Impact–ejecta maps and their heliocentric distances
   call get_maps_data(rhels, fnames)
   mapind1 = 1 ; mapind2 = 2
   call read_first_ratemap(fnames(mapind1), rhel1)
   call read_ratemap(fnames(mapind2), rhel2)

   !--- Loop over grain radii (β/μ_R) ---
   do i_R = 0, Nrgs
      density = 0d0

      if (i_R > 0) then
         call beta_from_Rg(beta, Rgs(i_R))
      else
         beta = 0d0   ! large grains approximation for Rg=99 μm bin
      end if
      muR = GMsun * (1d0 - beta)

      ! Time-window limits for contributing sources (same logic as 2-D)
      dtlim2 = resolution(1) * nx * (1d0 - centerpositionx) / AU / sources(1)%ud%umin
      dtlim3 = sqrt(2d0 * sources(1)%r**2 * resolution(1) / AU * &
                    (1d0 - centerpositionx) * nx / (GMsun - muR))
      idt = 1
      do while (tnow - sources(idt)%Tj > min(dtlim3, dtlim2))
         idt = idt + 1
      end do

      ! Walk forward along the trajectory, accumulating contributions
      do i_p = idt, Np-1
         ! Advance rate-map bracket if heliocentric distance decreased enough
         do while (sources(i_p)%r < rhel2 .and. mapind2 < Nmaps)
            rhel1  = rhel2
            mapind1 = mapind2
            mapind2 = mapind2 + 1
            call read_ratemap(fnames(mapind2), rhel2)
         end do
         call ratematr_interpolate(sources(i_p)%r, rhel1, rhel2)

         dt = tnow - sources(i_p)%Tj
         call runge_kutta_point_position(comet(i_p)%coords, comet(i_p)%Vastvec, muR, dt, cloudcentr)
         rMtmp = cloudcentr

!$omp parallel do default(none) private(i,j,k) collapse(3) schedule(static) &
!$omp& shared(nx,ny,nz,i_p,points,sources,tmp_res,muR,comet,dt,cloudcentr)
         do k = 1, nz
            do j = 1, ny
               do i = 1, nx
                  call hc_DUDI_simple_expansion(tmp_res(i,j,k), sources(i_p), dt, cloudcentr, points(i,j,k))
               end do
            end do
         end do
!$omp end parallel do

         density = density + tmp_res
      end do

      !--- Binary output (raw double-precision array), Fortran column-major order ---
      write(fnameout, '("results/Rg=", F5.2, "micron.bin")') Rgs(i_R)
      open(unit=77, file=trim(fnameout), access='stream', form='unformatted', status='replace')
      write(77) density
      close(77)
      write(*,*) '3D result written to ', trim(fnameout)
      if (i_R < Nrgs) write(*,*) 'continuing with next Rg...'
   end do

   deallocate(points, tmp_res, density, sources, comet)
end program phaethon_3d
