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

! File: select_method.f90
! Description: A routine that calculates dust number density using three
! different methods, compares the results, and provides recommendations on
! the applicability of each method.

! Author: Anastasiia Ershova
! E-mail: vveyzaa@gmail.com

program select_method
    use const
    use help
    use define_types
    use data_in
    use DUDIhc
    use data_out
    USE OMP_LIB
    
    implicit none
    real(8), parameter :: Rast = 0d0
    real(8), parameter :: Rast_AU = Rast / AU
    real(8), parameter :: Rgrain = 2d-6
    real, parameter :: accuracy = 5.0 ! %
    logical, parameter :: pericenter = .FALSE.
    character(len = 20) outformat
    integer nt1, nt2
    real(8) tnow, resolution(2), dtau
    integer i, ii, k
    real, allocatable, dimension(:,:) :: density_s, density_d, density_v
    real, allocatable, dimension(:,:) :: test_d, test_s
    real(8) muR, rtmp(3)
    type(source_properties) source
    type(position_in_space), allocatable, dimension(:,:) :: points
    type(ephemeris) comet
    character(len = 61) :: fname = './input_data_files/orbit_and_time_test.dat'
    
    nt1 = 200 ; nt2 = 200
    allocate(points(nt1,nt2), density_s(nt1, nt2), &
    density_d(nt1,nt2), density_v(nt1,nt2), test_s(nt1,nt2), test_d(nt1,nt2))
    
    ! inputing the parameters of the test case
    open(200, file = fname, status = 'old')
        read(200,*) comet%coords
        read(200,*) comet%Vastvec
        comet%Vast = norma3d(comet%Vastvec)
        read(200,*) tnow
        read(200,*) dtau
    close(200)
    
! define a point source ejecting dust uniformly to all the directions
! with the ejection speed distributed uniformly between 0.5 m/s and 50 m/s
    source%rrM = comet%coords
    source%r = norma3d(source%rrM)
    source%alphaM = acos(source%rrM(3) / source%r)
    source%betaM = atan(source%rrM(2), source%rrM(1))
    source%symmetry_axis = source%rrM / source%r
    source%zeta = 0d0
    source%eta = 0d0
    source%ud%ud_shape = 0
    source%ud%umin = 0.5d0 / AUdays2SI
    source%ud%umax = 50d0 / AUdays2SI
    source%ejection_angle_distr = 1
    source%Nparticles = 1d10
    source%Tj = 0d0
    source%dtau = 2d0 * dtau / s_in_day
    
    density_s = 0.0
    density_d = 0.0
    density_v = 0.0
    
    muR = GMsun
    
    ! the points where we conpute the number density are chosen so that
    ! they form a plain section through the middle of the dust cloud
    ! points with zero-density must present only in a small region in the middle
    resolution(1) = source%ud%umax * tnow * AU / dble(nt1)
    resolution(2) = resolution(1)
    ! find the coordinates of the cloud center
    call runge_kutta_point_position(comet%coords, comet%Vastvec, muR, tnow, rtmp)
    
    call orbital_plane_grid(points, nt1, nt2, resolution, comet, rtmp)
    ! compute the number density distribution in the defined plane
    ! using two different solutions
    !$OMP PARALLEL PRIVATE(i,ii) &
    !$OMP SHARED(points, source, density_s, density_d, density_v, muR, comet)
    !$OMP DO
    do i = 1, nt1
    do ii = 1, nt2
        call hc_DUDI_delta_ejection(density_d(i,ii), points(i,ii), &
                        source, muR, tnow, comet, Rast_AU)
                        
        call hc_DUDI_v_integration(density_v(i,ii), points(i,ii), &
                        source, muR, tnow, comet, Rast_AU, pericenter)

        call hc_DUDI_simple_expansion(density_s(i,ii), &
                        source, tnow, rtmp, points(i,ii))
        
    enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    
    ! exclude from the consideration the center part of the cloud
    ! where the methods may provide different results due to
    ! their different resolution
    k = int(source%ud%umin * tnow / (resolution(1) / AU)) + 1
    forall(i = (nt1/2-k):(nt1/2+k))
        forall(ii = (nt2/2-k):(nt2/2+k))
            density_d(i,ii) = 1.0
            density_v(i,ii) = 1.0
            density_s(i,ii) = 1.0
        endforall
    endforall
    
    call matrix_out('./results/test_simple_exp_meth.dat', density_s, nt1, nt2)
    call matrix_out('./results/test_delta-eject_meth.dat', density_d, nt1, nt2)
	call matrix_out('./results/test_v-integr_meth.dat', density_v, nt1, nt2)

    ! compute the matrix of delta-ejection solution deviations
    ! from the v-integration solution
    forall(i = 1:nt1)
        forall(ii = 1:nt2) test_d(i,ii) = density_d(i,ii) / density_v(i,ii)
       forall(ii = 1:nt2) test_s(i,ii) = density_s(i,ii) / density_d(i,ii)
    endforall
    test_s = test_s - 1.0
    test_d = test_d - 1.0
    call matrix_out('./results/test_simp_exp_vs_delta-eject.dat', test_s, nt1, nt2)
    call matrix_out('./results/test_delta-eject_vs_v-integr.dat', test_d, nt1, nt2)
    
    ! estimate the minimum and maximum deviations
    write(*,*) 'difference between the delta-ejection solution and v-integration solution'
    write(*,'(A3, x, f7.1, A1)') 'min', minval(test_d) * 100, '%'
    write(*,'(A3, x, f7.1, A1)') 'max', maxval(test_d) * 100, '%'
    
    write(*,*) 'difference between the simple expansion solution and delta-ejection solution'
    write(*,'(A3, x, f7.1, A1)') 'min', minval(test_s) * 100, '%'
    write(*,'(A3, x, f7.1, A1)') 'max', maxval(test_s) * 100, '%'

    ! recommend the delta ejection method if the deviations are small
    if(sum(abs(test_d)) / dble(nt1 * nt2) < accuracy * 1e-2 &
    .and. max(abs(minval(test_d) * 100), abs(maxval(test_d) * 100)) < accuracy) then
        write(*,*) 'delta-ejection method is applicable'
    ! otherwise, the v-integration solution is to be preferred
    else
        write(*,*) 'v-integration method is recommended'
    endif
	! recommend the simple expansion method if the deviations are small
    if(sum(abs(test_s)) / dble(nt1 * nt2) < accuracy * 1e-2 &
    .and. max(abs(minval(test_s) * 100), abs(maxval(test_s) * 100)) < accuracy) then
        write(*,*) 'simple expansion method is applicable too'
    ! otherwise, the v-integration solution is to be preferred
    else
        write(*,*) 'simple expansion method is NOT recommended'
    endif
    
    deallocate(points, test_s, test_d, density_d, density_s, density_v)
    
end
