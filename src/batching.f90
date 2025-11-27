! This file is a part of DUDI-heliocentric, the Fortran-90 implementation 
! of the two-body model for the dynamics of dust ejected from an atmosphereless
! body moving around the Sun
! Version 1.0.2
! This is free software. You can use and redistribute it 
! under the terms of the GNU General Public License (http://www.gnu.org/licenses/)
! If you do, please cite the following paper
! Anastasiia Ershova and Juergen Schmidt, 
! Two-body model for the spatial distribution of dust ejected from
! an atmosphereless body, 2021, A&A, 650, A186 
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
    use DUDIhc
    implicit none

    ! Method selector codes
    integer, parameter :: METHOD_SIMPLE_EXPANSION = 1
    integer, parameter :: METHOD_DELTA_EJECTION   = 2
    integer, parameter :: METHOD_V_INTEGRATION    = 3

contains

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
            !!$omp parallel do default(shared) private(i) schedule(dynamic)
            do i = 1, n_points
                ! hc_DUDI_simple_expansion(density, source, dt, cloudcentr, point)
                call hc_DUDI_simple_expansion( density(i), source, dt, &
                                               cloudcentr, points(i) )
            end do
            !!$omp end parallel do

        case (METHOD_DELTA_EJECTION)
            !!$omp parallel do default(shared) private(i) schedule(dynamic)
            do i = 1, n_points
                ! hc_DUDI_delta_ejection(density, point, source, muR, dt, comet, Rast_AU)
                call hc_DUDI_delta_ejection( density(i), points(i), source, &
                                             muR, dt, comet, Rast_AU )
            end do
            !!$omp end parallel do

        case (METHOD_V_INTEGRATION)
            !!$omp parallel do default(shared) private(i) schedule(dynamic)
            do i = 1, n_points
                ! hc_DUDI_v_integration(density, point, source, muR, tnow, comet, Rast_AU, pericenter)
                call hc_DUDI_v_integration( density(i), points(i), source, &
                                            muR, tnow, comet, Rast_AU, pericenter )
            end do
            !!$omp end parallel do

        case default
            ! Unknown method: set all densities to zero (or handle differently if you prefer)
            density(:) = 0.0

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
            !!$omp parallel do default(shared) private(i) schedule(dynamic)
            do i = 1, n_sources
                ! hc_DUDI_simple_expansion(density, source, dt, cloudcentr, point)
                call hc_DUDI_simple_expansion( density(i), sources(i), dt, &
                                               cloudcentr, point )
            end do
            !!$omp end parallel do

        case (METHOD_DELTA_EJECTION)
            !!$omp parallel do default(shared) private(i) schedule(dynamic)
            do i = 1, n_sources
                ! hc_DUDI_delta_ejection(density, point, source, muR, dt, comet, Rast_AU)
                call hc_DUDI_delta_ejection( density(i), point, sources(i), &
                                             muR, dt, comet, Rast_AU )
            end do
            !!$omp end parallel do

        case (METHOD_V_INTEGRATION)
            !!$omp parallel do default(shared) private(i) schedule(dynamic)
            do i = 1, n_sources
                ! hc_DUDI_v_integration(density, point, source, muR, tnow, comet, Rast_AU, pericenter)
                call hc_DUDI_v_integration( density(i), point, sources(i), &
                                            muR, tnow, comet, Rast_AU, pericenter )
            end do
            !!$omp end parallel do

        case default
            density(:) = 0.0

        end select

    end subroutine hc_DUDI_batch_sources

end module batching
