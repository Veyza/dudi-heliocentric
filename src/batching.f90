! This file is a part of DUDI-heliocentric, the Fortran-90 implementation 
! of the two-body model for the dynamics of dust ejected from an atmosphereless
! body moving around the Sun
! Version 1.1.0
! This is free software. You can use and redistribute it 
! under the terms of the GNU General Public License (http://www.gnu.org/licenses/)
! If you do, please cite the following papers
!
! Anastasiia Ershova and Juergen Schmidt, 
! Two-body model for the spatial distribution of dust ejected from
! an atmosphereless body, 2021, A&A, 650, A186 
! and Ershova, A., Schmidt, J., Liu, X., Szalay, J., Kimura, H., Hirai,
! T., Arai, T., and Kobayashi, M.,
! A computationally efficient semi-analytical model for the dust
! environment of comets and asteroids, A&A 693, A80 (2025).

! Module: batching.f90
! Description:
!   Batched wrappers over the three DUDI methods:
!     - hc_DUDI_simple_expansion
!     - hc_DUDI_delta_ejection
!     - hc_DUDI_v_integration
!   providing:
!     * parallelization over an array of points (fixed source),
!     * parallelization over an array of sources (fixed point),
!   with an integer parameter selecting the method.

module batching
    use define_types
    use help
    use DUDIhc
    implicit none

    ! Method selector codes
    integer, parameter :: METHOD_SIMPLE_EXPANSION = 1
    integer, parameter :: METHOD_DELTA_EJECTION   = 2
    integer, parameter :: METHOD_V_INTEGRATION    = 3

contains
    !==================================================================
    !> General batched computation over:
    !>   - n_points observation points
    !>   - Nt time steps
    !>   - Ns sources active at each time
    !>
    !>  sources(Nt, Ns) : all sources, grouped by time
    !>  comets(Nt)      : comet ephemeris at each time
    !>
    !>  For each time i_t:
    !>     dt = tnow - sources(i_t,1)%Tj
    !>     (simple-expansion only) cloudcentr = propagated comet position
    !>
    !>  Then contributions from all sources at all times are summed:
    !>
    !>      density(i_point) = sum_{t,s} method(point_i, source_{t,s}, ...)
    !>
    !>  method_id:
    !>    = METHOD_SIMPLE_EXPANSION : hc_DUDI_simple_expansion
    !>    = METHOD_DELTA_EJECTION   : hc_DUDI_delta_ejection
    !>    = METHOD_V_INTEGRATION    : hc_DUDI_v_integration
    !==================================================================
    subroutine hc_DUDI_batch_sources_points( &
        n_points, Nt, Ns, density, points, sources, &
        muR, tnow, comets, Rast_AU, pericenter, method_id )
		use distributions_fun
        implicit none

        integer, intent(in) :: n_points, Nt, Ns
        real,    intent(out) :: density(n_points)
        type(position_in_space), intent(in) :: points(n_points)
        type(source_properties), intent(in) :: sources(Nt, Ns)
        real(8), intent(in) :: muR, tnow, Rast_AU
        type(ephemeris), intent(in) :: comets(Nt)
        logical, intent(in) :: pericenter
        integer, intent(in) :: method_id

        integer :: i_t, i_s, i
        real(8) :: dt
        real(8) :: cloudcentr(3)
        real    :: tmp(n_points)
		
        density(:) = 0.0

        do i_t = 1, Nt
           ! dt for this time slice (all Ns sources share same ejection time)
           dt = tnow - sources(i_t, 1)%Tj

           select case (method_id)

           case (METHOD_SIMPLE_EXPANSION)
              ! Compute cloud center for this time
              call runge_kutta_point_position( comets(i_t)%coords, &
                                               comets(i_t)%Vastvec, &
                                               muR, dt, cloudcentr )
			  rMtmp = cloudcentr
              do i_s = 1, Ns
              !$OMP PARALLEL PRIVATE(i) &
			  !$OMP SHARED(points, sources, density, muR, comets, dt, i_t, i_s, tmp, cloudcentr)
			  !$OMP DO
                 do i = 1, n_points
                    call hc_DUDI_simple_expansion(tmp(i), sources(i_t, i_s), dt, &
                                                   cloudcentr, points(i))
                 end do
              !$OMP END DO
			  !$OMP END PARALLEL
              density = density + tmp
!~               write(*,*) tmp(1), density(1)
!~               write(*,*) sources(i_t, i_s)
              end do

           case (METHOD_DELTA_EJECTION)
              do i_s = 1, Ns
              !$OMP PARALLEL PRIVATE(i) &
			  !$OMP SHARED(points, sources, density, muR, comets, dt, i_t, tmp)
		  	  !$OMP DO
                 do i = 1, n_points
                    call hc_DUDI_delta_ejection( tmp(i), points(i), sources(i_t, i_s), &
                                                 muR, dt, comets(i_t), Rast_AU )
                 end do
              !$OMP END DO
			  !$OMP END PARALLEL
              density = density + tmp
              end do

           case (METHOD_V_INTEGRATION)
              do i_s = 1, Ns
              !$OMP PARALLEL PRIVATE(i) &
			  !$OMP SHARED(points, sources, density, muR, comets, dt, i_t, tmp, pericenter)
			  !$OMP DO
                 do i = 1, n_points
                    call hc_DUDI_v_integration( tmp(i), points(i), sources(i_t, i_s), &
                                                muR, tnow, comets(i_t), Rast_AU, pericenter )
                 end do
              !$OMP END DO
			  !$OMP END PARALLEL
              density = density + tmp
              end do

           case default
            ! Unknown method: set all densities to negative
            density(:) = -10.0

           end select
           
        end do

    end subroutine hc_DUDI_batch_sources_points



    !==================================================================
    !> Batched computation over an array of points for a single source.
    !!
    !!  For each i = 1..n_points, computes density(i) at points(i)
    !!  due to the same source, using the selected method.
    !!
    !!  method_id:
    !!    = METHOD_SIMPLE_EXPANSION : hc_DUDI_simple_expansion
    !!    = METHOD_DELTA_EJECTION   : hc_DUDI_delta_ejection
    !!    = METHOD_V_INTEGRATION    : hc_DUDI_v_integration
    !==================================================================
    subroutine hc_DUDI_batch_points( n_points, density, points, source, &
                                     muR, tnow, dt, comet, Rast_AU,      &
                                     pericenter, cloudcentr, method_id )

        implicit none

        integer, intent(in) :: n_points
        real,    intent(out) :: density(n_points)
        type(position_in_space), intent(in) :: points(n_points)
        type(source_properties), intent(in) :: source
        real(8), intent(in) :: muR, tnow, dt, Rast_AU
        type(ephemeris), intent(in) :: comet
        logical, intent(in) :: pericenter
        real(8), intent(in) :: cloudcentr(3)
        integer, intent(in) :: method_id

        integer :: i

        select case (method_id)

        case (METHOD_SIMPLE_EXPANSION)
            !$OMP PARALLEL PRIVATE(i) &
			!$OMP SHARED(points, source, density, muR, comet, dt, cloudcentr)
			!$OMP DO
            do i = 1, n_points
                ! hc_DUDI_simple_expansion(density, source, dt, cloudcentr, point)
                call hc_DUDI_simple_expansion( density(i), source, dt, &
                                               cloudcentr, points(i) )
            end do
			!$OMP END DO
			!$OMP END PARALLEL

        case (METHOD_DELTA_EJECTION)
            !$OMP PARALLEL PRIVATE(i) &
			!$OMP SHARED(points, source, density, muR, comet, dt, Rast_AU)
			!$OMP DO
            do i = 1, n_points
                ! hc_DUDI_delta_ejection(density, point, source, muR, dt, comet, Rast_AU)
                call hc_DUDI_delta_ejection( density(i), points(i), source, &
                                             muR, dt, comet, Rast_AU )
            end do
			!$OMP END DO
			!$OMP END PARALLEL

        case (METHOD_V_INTEGRATION)
            !$OMP PARALLEL PRIVATE(i) &
			!$OMP SHARED(points, source, density, muR, comet, dt, tnow, Rast_AU, pericenter)
			!$OMP DO
            do i = 1, n_points
                ! hc_DUDI_v_integration(density, point, source, muR, tnow, comet, Rast_AU, pericenter)
                call hc_DUDI_v_integration( density(i), points(i), source, &
                                            muR, tnow, comet, Rast_AU, pericenter )
            end do
			!$OMP END DO
			!$OMP END PARALLEL

        case default
            ! Unknown method: set all densities to negative
            density(:) = -10.0

        end select

    end subroutine hc_DUDI_batch_points

    !==================================================================
    !> Batched computation over an array of sources for a single point.
    !!
    !!  For each i = 1..n_sources, computes density(i) at the same point
    !!  due to sources(i), using the selected method.
    !!
    !!  method_id:
    !!    = METHOD_SIMPLE_EXPANSION : hc_DUDI_simple_expansion
    !!    = METHOD_DELTA_EJECTION   : hc_DUDI_delta_ejection
    !!    = METHOD_V_INTEGRATION    : hc_DUDI_v_integration
    !==================================================================
    subroutine hc_DUDI_batch_sources( n_sources, density, point, sources, &
                                      muR, tnow, dt, comet, Rast_AU,      &
                                      pericenter, cloudcentr, method_id )

        implicit none

        integer, intent(in) :: n_sources
        real,    intent(out) :: density(n_sources)
        type(position_in_space), intent(in) :: point
        type(source_properties), intent(in) :: sources(n_sources)
        real(8), intent(in) :: muR, tnow, dt, Rast_AU
        type(ephemeris), intent(in) :: comet
        logical, intent(in) :: pericenter
        real(8), intent(in) :: cloudcentr(3)
        integer, intent(in) :: method_id

        integer :: i

        select case (method_id)

        case (METHOD_SIMPLE_EXPANSION)
            !$omp parallel do default(shared) private(i) schedule(dynamic)
            do i = 1, n_sources
                ! hc_DUDI_simple_expansion(density, source, dt, cloudcentr, point)
                call hc_DUDI_simple_expansion( density(i), sources(i), dt, &
                                               cloudcentr, point )
            end do
            !$omp end parallel do

        case (METHOD_DELTA_EJECTION)
            !$omp parallel do default(shared) private(i) schedule(dynamic)
            do i = 1, n_sources
                ! hc_DUDI_delta_ejection(density, point, source, muR, dt, comet, Rast_AU)
                call hc_DUDI_delta_ejection( density(i), point, sources(i), &
                                             muR, dt, comet, Rast_AU )
            end do
            !$omp end parallel do

        case (METHOD_V_INTEGRATION)
            !$omp parallel do default(shared) private(i) schedule(dynamic)
            do i = 1, n_sources
                ! hc_DUDI_v_integration(density, point, source, muR, tnow, comet, Rast_AU, pericenter)
                call hc_DUDI_v_integration( density(i), point, sources(i), &
                                            muR, tnow, comet, Rast_AU, pericenter )
            end do
            !$omp end parallel do

        case default
            density(:) = 0.0

        end select

    end subroutine hc_DUDI_batch_sources

end module batching
