! python_interface/fortran_bridge/py_bridge.f90
!
! C-interoperable wrappers around DUDIhc.
! - Only interoperable scalars/1D arrays cross the boundary.
! - We reconstruct derived types locally and call the real routines.
! - Fortran LOGICAL is replaced by INTEGER(C_INT) on the C boundary.
! - DUDIhc returns density as REAL(4); we convert to REAL(8) on output.

module py_dudihc_bridge
  use iso_c_binding, only: c_int, c_double
  use iso_fortran_env, only: real32, real64
  use define_types, only: &
       position_in_space, &
       source_properties, &
       ejection_speed_properties, &
       ephemeris
  use DUDIhc, only: &
       hc_DUDI_v_integration, &
       hc_DUDI_delta_ejection, &
       hc_DUDI_simple_expansion
  implicit none
contains

  subroutine py_hc_v_integration( &
       point_r, point_alpha, point_beta, point_rvector,             &
       src_r, src_alphaM, src_betaM, src_rrM, src_zeta, src_eta,    &
       src_axis, src_eject_distr, src_ud_shape, src_umin, src_umax, &
       comet_coords, comet_vastvec, comet_vast,                     &
       muR, tnow, Rast_AU, pericenter_c, density) bind(C, name="py_hc_v_integration")

    ! C boundary types
    real(c_double), intent(in) :: point_r, point_alpha, point_beta
    real(c_double), intent(in) :: point_rvector(3)

    real(c_double), intent(in) :: src_r, src_alphaM, src_betaM, src_zeta, src_eta
    real(c_double), intent(in) :: src_rrM(3), src_axis(3)
    integer(c_int), intent(in) :: src_eject_distr, src_ud_shape
    real(c_double), intent(in) :: src_umin, src_umax

    real(c_double), intent(in) :: comet_coords(3), comet_vastvec(3)
    real(c_double), intent(in) :: comet_vast

    real(c_double), intent(in) :: muR, tnow, Rast_AU
    integer(c_int), intent(in) :: pericenter_c   ! 0/1 instead of LOGICAL

    real(c_double), intent(out) :: density       ! Python sees a double

    ! Local Fortran-side types
    type(position_in_space)          :: point
    type(source_properties)                :: source
    type(ejection_speed_properties)  :: ud
    type(ephemeris)                  :: comet
    logical                          :: pericenter
    real(real32)                     :: density_sp

    ! Build POINT
    point%r       = real(point_r,   kind=real64)
    point%alpha   = real(point_alpha, kind=real64)
    point%beta    = real(point_beta,  kind=real64)
    point%rvector = real(point_rvector, kind=real64)

    ! Build SOURCE
    source%r      = real(src_r,      kind=real64)
    source%alphaM = real(src_alphaM, kind=real64)
    source%betaM  = real(src_betaM,  kind=real64)
    source%rrM    = real(src_rrM,    kind=real64)
    source%zeta   = real(src_zeta,   kind=real64)
    source%eta    = real(src_eta,    kind=real64)
    source%symmetry_axis        = real(src_axis, kind=real64)
    source%ejection_angle_distr = int(src_eject_distr, kind=kind(source%ejection_angle_distr))
    ud%ud_shape   = int(src_ud_shape, kind=kind(ud%ud_shape))
    ud%umin       = real(src_umin, kind=real64)
    ud%umax       = real(src_umax, kind=real64)
    source%ud     = ud

    ! Build COMET
    comet%coords  = real(comet_coords,  kind=real64)
    comet%Vastvec = real(comet_vastvec, kind=real64)
    comet%Vast    = real(comet_vast,    kind=real64)

    pericenter = (pericenter_c /= 0)

    call hc_DUDI_v_integration(density_sp, point, source, real(muR,kind=real64), real(tnow,kind=real64), &
                               comet, real(Rast_AU,kind=real64), pericenter)

    density = real(density_sp, kind=real64)
  end subroutine py_hc_v_integration


  subroutine py_hc_delta_ejection( &
       point_r, point_alpha, point_beta, point_rvector,             &
       src_r, src_alphaM, src_betaM, src_rrM, src_zeta, src_eta,    &
       src_axis, src_eject_distr, src_ud_shape, src_umin, src_umax, &
       comet_coords, comet_vastvec, comet_vast,                     &
       muR, dt, Rast_AU, density) bind(C, name="py_hc_delta_ejection")

    real(c_double), intent(in) :: point_r, point_alpha, point_beta
    real(c_double), intent(in) :: point_rvector(3)

    real(c_double), intent(in) :: src_r, src_alphaM, src_betaM, src_zeta, src_eta
    real(c_double), intent(in) :: src_rrM(3), src_axis(3)
    integer(c_int), intent(in) :: src_eject_distr, src_ud_shape
    real(c_double), intent(in) :: src_umin, src_umax

    real(c_double), intent(in) :: comet_coords(3), comet_vastvec(3)
    real(c_double), intent(in) :: comet_vast

    real(c_double), intent(in) :: muR, dt, Rast_AU
    real(c_double), intent(out) :: density

    type(position_in_space)          :: point
    type(source_properties)                :: source
    type(ejection_speed_properties)  :: ud
    type(ephemeris)                  :: comet
    real(real32)                     :: density_sp

    point%r       = real(point_r,   kind=real64)
    point%alpha   = real(point_alpha, kind=real64)
    point%beta    = real(point_beta,  kind=real64)
    point%rvector = real(point_rvector, kind=real64)

    source%r      = real(src_r,      kind=real64)
    source%alphaM = real(src_alphaM, kind=real64)
    source%betaM  = real(src_betaM,  kind=real64)
    source%rrM    = real(src_rrM,    kind=real64)
    source%zeta   = real(src_zeta,   kind=real64)
    source%eta    = real(src_eta,    kind=real64)
    source%symmetry_axis        = real(src_axis, kind=real64)
    source%ejection_angle_distr = int(src_eject_distr, kind=kind(source%ejection_angle_distr))
    ud%ud_shape   = int(src_ud_shape, kind=kind(ud%ud_shape))
    ud%umin       = real(src_umin, kind=real64)
    ud%umax       = real(src_umax, kind=real64)
    source%ud     = ud

    comet%coords  = real(comet_coords,  kind=real64)
    comet%Vastvec = real(comet_vastvec, kind=real64)
    comet%Vast    = real(comet_vast,    kind=real64)

    call hc_DUDI_delta_ejection(density_sp, point, source, real(muR,kind=real64), real(dt,kind=real64), &
                                comet, real(Rast_AU,kind=real64))

    density = real(density_sp, kind=real64)
  end subroutine py_hc_delta_ejection


  subroutine py_hc_simple_expansion( &
       point_r, point_alpha, point_beta, point_rvector,             &
       src_r, src_alphaM, src_betaM, src_rrM, src_zeta, src_eta,    &
       src_axis, src_eject_distr, src_ud_shape, src_umin, src_umax, &
       cloudcentr, dt, density) bind(C, name="py_hc_simple_expansion")

    real(c_double), intent(in) :: point_r, point_alpha, point_beta
    real(c_double), intent(in) :: point_rvector(3)

    real(c_double), intent(in) :: src_r, src_alphaM, src_betaM, src_zeta, src_eta
    real(c_double), intent(in) :: src_rrM(3), src_axis(3)
    integer(c_int), intent(in) :: src_eject_distr, src_ud_shape
    real(c_double), intent(in) :: src_umin, src_umax

    real(c_double), intent(in) :: cloudcentr(3)
    real(c_double), intent(in) :: dt
    real(c_double), intent(out):: density

    type(position_in_space)          :: point
    type(source_properties)                :: source
    type(ejection_speed_properties)  :: ud
    real(real32)                     :: density_sp

    point%r       = real(point_r,   kind=real64)
    point%alpha   = real(point_alpha, kind=real64)
    point%beta    = real(point_beta,  kind=real64)
    point%rvector = real(point_rvector, kind=real64)

    source%r      = real(src_r,      kind=real64)
    source%alphaM = real(src_alphaM, kind=real64)
    source%betaM  = real(src_betaM,  kind=real64)
    source%rrM    = real(src_rrM,    kind=real64)
    source%zeta   = real(src_zeta,   kind=real64)
    source%eta    = real(src_eta,    kind=real64)
    source%symmetry_axis        = real(src_axis, kind=real64)
    source%ejection_angle_distr = int(src_eject_distr, kind=kind(source%ejection_angle_distr))
    ud%ud_shape   = int(src_ud_shape, kind=kind(ud%ud_shape))
    ud%umin       = real(src_umin, kind=real64)
    ud%umax       = real(src_umax, kind=real64)
    source%ud     = ud

    call hc_DUDI_simple_expansion(density_sp, source, real(dt,kind=real64), real(cloudcentr,kind=real64), point)

    density = real(density_sp, kind=real64)
  end subroutine py_hc_simple_expansion

end module py_dudihc_bridge
