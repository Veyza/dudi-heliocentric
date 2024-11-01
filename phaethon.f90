! This file is a part of the Fortran-95 implementation 
! of the two-body model for dust dynamics
! Version 1.0.0
! This is free software. You can use and redistribute it 
! under the terms of the GNU General Public License (http://www.gnu.org/licenses/)
! If you do, please cite the following paper
! **REFERENCE TO THE PAPER IN THE RIGHT FORMAT SO THAT IT COULD BE JUST COPIED AND PASTED TO THE REFERENCE LIST**

! File: main_program.f95
! Description: A suggested template for application of the model

! Author: Anastasiia Ershova
! E-mail: vveyzaa@gmail.com

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
    integer, parameter :: Neph = 144                                ! this many positions are in the file with ephemeridae
    integer, parameter :: Nlin = 10									! this many positions - 1 are interpolated between the ephemeridae
	integer, parameter :: Np = (Neph-1) * Nlin + 1					! number of points along the asteroid trajectory
    real(8), parameter :: Rast = 0d0
    real(8), parameter :: Rast_AU = Rast / AU
    real(8), parameter :: centerpositionx = 0.5d0
    real(8), parameter :: centerpositiony = 0.5d0
    integer, parameter :: Nmaps = 9
    integer, parameter :: Nrgs = 13
    logical, parameter :: pericenter = .FALSE.
    character(len = 20) outformat
    integer nt1, nt2, inds(3)
    real(8) tnow, resolution(2), beta
    integer i_p, i, ii, i_R, idt
    real(8) shift, dtlim2, dtlim3
    real(8) :: Rgs(Nrgs) = (/0.1d0, 0.2d0, 0.3d0, &
                          0.4d0, 0.5d0, 0.65d0, 0.85d0, 1d0, &
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
    
    nt1 = 400 ; nt2 = 400
    allocate(points(nt1,nt2), tmp_res(nt1,nt2), density(nt1,nt2), &
             corr_res(nt1, nt2), sources(Np), comet(Np))
    
    call get_moving_sources(fname, Np, Neph, Nlin, sources, comet)
    
    tnow = sources(Np)%Tj

    resolution(1) = 5d3
    resolution(2) = resolution(1)
    
    call get_points(points, nt1, nt2, resolution, comet(Np)%coords, &
                comet(Np)%Vastvec, comet(Np)%Vast, &
                centerpositionx, centerpositiony)


    call get_maps_data(rhels, fnames)
    density = 0.0
    do i_R = 0, Nrgs

         mapind1 = 1 ; mapind2 = mapind1 + 1
         call read_first_ratemap(fnames(mapind1), rhel1)
         call read_ratemap(fnames(mapind2), rhel2)
        
        density = 0d0
        if(i_R > 0) then
            call beta_from_Rg(beta, Rgs(i_R))
        else
            beta = 0d0
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
        do i_p = idt, Np-1
            do while(sources(i_p)%r < rhel2)
                rhel1 = rhel2
                mapind1 = mapind2
                mapind2 = mapind2 + 1
                rmap1 = rmap2
                call read_ratemap(fnames(mapind2), rhel2)
                write(*,*) 'rhel', rhel1, sources(i_p)%r, rhel2
            enddo
            call ratematr_interpolate(sources(i_p)%r, rhel1, rhel2)
            dt = tnow - sources(i_p)%Tj
			call runge_kutta_point_position(comet(i_p)%coords, &
								comet(i_p)%Vastvec, muR, dt, cloudcentr)
			rMtmp = cloudcentr
			!$OMP PARALLEL PRIVATE(i,ii) &
			!$OMP SHARED(i_p, i_R, points, sources, tmp_res, Rgs, muR, comet)
			!$OMP DO
			do ii = 1, nt2
			do i = 1, nt1
				
				call hc_DUDI_simple_expansion(tmp_res(i,ii), &
							sources(i_p), dt, cloudcentr, points(i,ii))

			enddo
			enddo
			!$OMP END DO
			!$OMP END PARALLEL
			density = density + tmp_res
        enddo
        write(*,*) minval(density), maxval(density)
        if(i_R > 0) then
			write(fnameout, '("Rg=", F4.2, "micron.dat")') Rgs(i_R)
        else
			write(fnameout, '("Rg=", F4.2, "micron.dat")') 100d0
        endif
        call matrix_out(trim(fnameout), density, nt1, nt2)
    enddo

    deallocate(points, tmp_res, density, sources, comet)
    
end
