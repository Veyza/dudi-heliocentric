! This file is a part of DUDI-heliocentric, the Fortran-90 implementation 
! of the two-body model for the dynamics of dust ejected from an atmosphereless
! body moving around the Sun
! Version 1.0.1
! This is free software. You can use and redistribute it 
! under the terms of the GNU General Public License (http://www.gnu.org/licenses/)
! If you do, please cite the following paper
! Anastasiia Ershova and Jürgen Schmidt, 
! Two-body model for the spatial distribution of dust ejected from
! an atmosphereless body, 2021, A&A, 650, A186 

! Author: Anastasiia Ershova
! E-mail: vveyzaa@gmail.com

! File: phaethon_input.f90
! Description: This is the main program managing the modeling of dust
! ejection from NEA 3200 Phaethon at the orbital phase at which it is
! planned to be visited by the DESTINY+ mission using the simple 
! expansion method of DUDI-heliocentric 
! (see Sect. 4.1.1 of Ershova et. al, 2025)

program phaethon
    use const
    use help
    use distributions_fun
    use define_types
    use DUDIhc
    use phaethon_input
    use data_out
    USE OMP_LIB
    
    implicit none
    integer, parameter :: Neph = 2000  ! this many positions are in the file with ephemeridae
    integer, parameter :: Nlin = 10  ! this many positions - 1 are interpolated between the ephemeridae
  integer, parameter :: Np = (Neph-1) * Nlin + 1  ! number of points along the asteroid trajectory
    real(8), parameter :: Rast = 0d0
    real(8), parameter :: Rast_AU = Rast / AU
    real(8), parameter :: centerpositionx = 0.5d0
    real(8), parameter :: centerpositiony = 0.5d0
    integer, parameter :: Nmaps = 4
    integer, parameter :: Nrgs = 13
    logical, parameter :: pericenter = .FALSE.
    character(len = 20) outformat
    integer nt1, nt2, inds(3)
    real(8) tnow, resolution(2), beta
    integer i_p, i, ii, i_R, idt
    real(8) shift, dtlim2, dtlim3
    real(8) :: Rgs(Nrgs) = (/0.1d0, 0.2d0, 0.3d0, &
                          0.42d0, 0.55d0, 0.67d0, 0.85d0, 1d0, &
                          1.2d0, 2.5d0, 4d0, 6d0, 10d0/)
    real, allocatable, dimension(:,:) :: density, tmp_res, corr_res
    real(8) muR, Rgrain, tmp, rtmp(3), dphi, timeres, dRg, ddist, dt
    type(source_properties), allocatable, dimension(:) :: sources
    type(position_in_space), allocatable, dimension(:,:) :: points
    type(ephemeris), allocatable, dimension(:) :: comet
    real(8) cloudcentr(3), tmpvec(3)
    character(len = 5) chNlin
    character(len = 66) :: fname = 'input_data_files/Phaethon_2025-02-22_last_int=10min_ECLIPJ2000.dat'
    character(len = 50) fnameout
    character(len = 93), dimension(Nmaps) :: fnames
    real(8) rhels(Nmaps), rhel1, rhel2
    integer mapind1, mapind2
    
  ! We will compute the number density in a planar rectangular grid
    ! nt1 and nt2 are the number of the grid nodes along the vertical
    ! and horizontal directions
    nt1 = 400 ; nt2 = 400
    allocate(points(nt1,nt2), tmp_res(nt1,nt2), density(nt1,nt2), &
             corr_res(nt1, nt2), sources(Np), comet(Np))
    
    ! input sources parameters
    call get_moving_sources(fname, Np, Nlin, sources, comet)
    
    ! the moment for which we compute the number density
    tnow = sources(Np)%Tj
  
  ! resolution(1) and (2) are the distances between the grid nodes in 
  ! horizontal and vertical directions
    resolution(1) = 5d3
    resolution(2) = resolution(1)
    ! Generating the list of points where we compute number density
    call get_points(points, nt1, nt2, resolution, comet(Np)%coords, &
                centerpositionx, centerpositiony)

  ! Load the matrices with number density of impact ejecta (Szalay et al., 2019)
    call get_maps_data(rhels, fnames)
    
    density = 0.0
    ! Loop over particle radii (different values of beta and muR)
    do i_R = 0, Nrgs

         mapind1 = 1 ; mapind2 = mapind1 + 1
         call read_first_ratemap(fnames(mapind1), rhel1)
         call read_ratemap(fnames(mapind2), rhel2)
        
        density = 0d0
        if(i_R > 0) then
        ! obtain beta for the given particle radius
            call beta_from_Rg(beta, Rgs(i_R))
        else
            beta = 0d0  ! approximation for large grains of 100 um radius
        endif
        muR = GMsun * (1d0 - beta)        
        ! when the particles ejected with umin leave the FoV
        dtlim2 = resolution(1) * nt1 * (1d0 - centerpositionx) / AU / sources(1)%ud%umin
        ! when the edge of the prime cloud leaves the FoV
        dtlim3 = sqrt(2d0 * sources(1)%r**2 * resolution(1) / AU &
                * (1d0 - centerpositionx) * nt1 / (GMsun - muR))                    
        idt = 1
        ! find how far from tnow the dust that is still in the
        ! considered vicinity of the asteroid was ejected
        do while(tnow - sources(idt)%Tj > min(dtlim3, dtlim2))
            idt = idt + 1
        enddo
        write(*,*) 'start index', idt
        ! Loop over the consequently active sources along the asteroid
        ! trajectory
        do i_p = idt, Np-1
      ! use the impact-ejecta map that corresponds to the current
      ! heliocentric distance
            do while(sources(i_p)%r < rhel2)
                rhel1 = rhel2
                mapind1 = mapind2
                mapind2 = mapind2 + 1
                rmap1 = rmap2
                call read_ratemap(fnames(mapind2), rhel2)
            enddo
            call ratematr_interpolate(sources(i_p)%r, rhel1, rhel2)
            dt = tnow - sources(i_p)%Tj
      call runge_kutta_point_position(comet(i_p)%coords, &
                comet(i_p)%Vastvec, muR, dt, cloudcentr)
      rMtmp = cloudcentr
    !$omp parallel do default(none) private(i,ii) collapse(2) schedule(static) &
    !$omp& shared(nt1,nt2, i_p, points, sources, tmp_res, Rgs, muR, comet, dt, cloudcentr)
    do ii = 1, nt2
        do i = 1, nt1
            call hc_DUDI_simple_expansion(tmp_res(i,ii), sources(i_p), dt, cloudcentr, points(i,ii))
        end do
    end do
    !$omp end parallel do
      density = density + tmp_res
        enddo
        if(i_R > 0) then
      write(fnameout, '("results/Rg=", F5.2, "micron.dat")') Rgs(i_R)
        else
      write(fnameout, '("results/Rg=", F6.2, "micron.dat")') 100d0
        endif
        call matrix_out(trim(fnameout), density, nt1, nt2)
    write(*,*) 'result is in the file ', fnameout
    if(i_R < Nrgs) then
      write(*,*) 'calculations continue'
    endif
    enddo

    deallocate(points, tmp_res, density, sources, comet)
    
end
