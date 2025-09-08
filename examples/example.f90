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

! File: example.f90
! Description: A suggested template for application of the model
! This file utilizes the "data_in" module to compute the gravitational 
! parameter 'muR', reduced for solar radiation pressure, and to load asteroid
! position and velocity vectors at specified times. It inputes the properties
! of model dust sources, with 'Np' defining the count of asteroid states and
! 'Ns' the number of sources at each state. The last recorded asteroid position
! centers a grid, sized by 'nt1' and 'nt2', where dust density is calculated
! as a sum of dust from each of the Np*Ns sources by the delta-ejection method
! implemented in "DUDIhc". The 'muR' is calculated based on dust grain radius
! 'Rg' and radiation pressure efficiency 'Qpr'.
! Outputs are saved to "test_result.dat" in the "results" directory through
! a subroutine from the module "data_out".
! OpenMP library is leveraged to speed up the calculations
program example
    use const
    use define_types
    use data_in
    use data_out
    use DUDIhc
    USE OMP_LIB
    implicit none
    ! number of points along the asteroid trajectory from which dust is ejected
    integer, parameter :: Np = 41
    ! number of sources representing the asteroid at each point
    integer, parameter :: Ns = 50
    real(8), parameter :: Rast = 5d3            ! asteroid radius, meters
    real(8), parameter :: Rast_AU = Rast / AU   ! asteroid radius, AU
    ! the values of Qpr = 0.5 and Rg = 0.29 are set purposefully
    ! to obtain \beta = 0.4
    real(8), parameter :: Qpr = 0.5d0           ! radiation pressure efficiency
    real(8), parameter :: Rg = 0.29d-6          ! dust grain radius, meters
    integer, parameter :: nt1 = 200
    integer, parameter :: nt2 = 200
    ! distance between the grid nodes, meters
    real(8) :: resolution(2) = (/2d3, 2d3/)     
    integer ip, is, i, ii
    real density(nt1, nt2), tmpres(nt1, nt2)
    real(8) muR, tnow, dt
    type(source_properties) sources(Np, Ns)
    type(position_in_space) points(nt1, nt2)
    type(ephemeris) comet(Np)
    character(len = 50) :: fname = './input_data_files/ephemeridae.dat'
    
    ! calculating the parameter \mu_R = GM_sun * (1 - \beta)
    ! for the given grain radius and radiation pressure efficiency
    call reduced_gravitational_parameter(Rg, Qpr, muR)
    
    ! inputing the ephemeridae of the asteroid (`comet´)
    ! and properties of the model sources (`sources´)
    call get_sources(fname, Np, Ns, sources, comet, Rast_AU)
    
    ! a grid of points centered at the last inputed position of the asteroid
    ! `comet(Np)%coords´
    call orbital_plane_grid(points, nt1, nt2, resolution, &
                              comet(Np), comet(Np)%coords)
    
    ! sources(Np,1)%Tj is the moment at which the asteroid is at comet(Np)%coords
    tnow = sources(Np,1)%Tj
    
    density = 0.0
    
    ! the sources with the index Np were active 0 seconds before tnow,
    ! so we do not calculate number density of dust from them
    do ip = 1, Np-1
        ! the time passed from the ejection by the sources with the index `ip´
        dt = tnow - sources(ip,1)%Tj
        do is = 1, Ns
        !$OMP PARALLEL PRIVATE(i,ii) &
        !$OMP SHARED(points, sources, density, muR, comet)
        !$OMP DO
            do ii = 1, nt2
            do i = 1, nt1
        ! for our example case delta-ejection method has a sufficient accuracy
                call hc_DUDI_delta_ejection(tmpres(i,ii), points(i,ii), &
                        sources(ip,is), muR, dt, &
                        comet(ip), Rast_AU)
            enddo
            enddo
            !$OMP END DO
            !$OMP END PARALLEL
        ! adding the number density from the source with indexes `ip´ and `is´
            density = density + tmpres
        enddo
    enddo
    
    call matrix_out('./results/result.dat', &
                          density, nt1, nt2)

end
