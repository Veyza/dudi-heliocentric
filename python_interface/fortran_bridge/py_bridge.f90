! This file is a part of DUDI-heliocentric, the Fortran-95 implementation 
! of the two-body model for the dynamics of dust ejected from an atmosphereless
! body moving around the Sun
! Version 1.1.1
! This is free software. You can use and redistribute it 
! under the terms of the GNU General Public License (http://www.gnu.org/licenses/)
! If you do, please cite the following paper
! Anastasiia Ershova and Jürgen Schmidt, 
! Two-body model for the spatial distribution of dust ejected from
! an atmosphereless body, 2021, A&A, 650, A186 

! Author: Anastasiia Ershova
! E-mail: vveyzaa@gmail.com

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
  use batching, only: hc_DUDI_batch_points, hc_DUDI_batch_sources, &
                      hc_DUDI_batch_sources_points
  use define_types, only: &
       position_in_space, &
       source_properties, &
       ejection_speed_properties, &
       ephemeris
  use DUDIhc, only: &
       hc_DUDI_v_integration, &
       hc_DUDI_delta_ejection, &
       hc_DUDI_simple_expansion
  use distributions_fun, only: nlats, nlons, lonmax, lonmin, &
                         lats, lons, rmap1, rmap2, ratemap, rMtmp, &
                         ratematr_interpolate, read_ratemap, &
                         set_tabulated_fu, set_tabulated_fpsi
  implicit none


contains
    logical function env_enabled(name)
    use iso_c_binding, only: c_int
    implicit none
    character(len=*), intent(in) :: name
    character(len=8) :: val
    integer :: stat, lenval
    env_enabled = .false.
    call get_environment_variable(name, value=val, length=lenval, status=stat, trim_name=.true.)
    if (stat == 0) then
       if (lenval > 0) then
          ! treat any non-'0' first char as enabled
          env_enabled = (val(1:1) /= '0')
       end if
    end if
  end function env_enabled

  subroutine py_hc_v_integration( &
        point_r, point_alpha, point_beta, point_rvector,             &
        src_r, src_alphaM, src_betaM, src_rrM, src_zeta, src_eta,    &
        src_axis, src_eject_distr, src_ud_shape, src_umin, src_umax, &
        src_Nparticles, src_Tj, src_dtau,                            &
        comet_coords, comet_vastvec, comet_vast,                     &
        muR, tnow, Rast_AU, pericenter_c, density) bind(C, name="py_hc_v_integration")
    use, intrinsic :: iso_fortran_env, only: real64, real32, output_unit
    ! C boundary types
    real(c_double), value, intent(in) :: point_r, point_alpha, point_beta
    real(c_double), intent(in) :: point_rvector(3)

    real(c_double), value, intent(in) :: src_r, src_alphaM, src_betaM, src_zeta, src_eta
    real(c_double), intent(in) :: src_rrM(3), src_axis(3)
    integer(c_int), value, intent(in) :: src_eject_distr, src_ud_shape
    real(c_double), value, intent(in) :: src_umin, src_umax
    real(c_double), value, intent(in) :: src_Nparticles, src_Tj, src_dtau

    real(c_double), intent(in) :: comet_coords(3), comet_vastvec(3)
    real(c_double), value, intent(in) :: comet_vast

    real(c_double), value, intent(in) :: muR, tnow, Rast_AU
    integer(c_int), value, intent(in) :: pericenter_c   ! 0/1 instead of LOGICAL

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
    source%Nparticles = real(src_Nparticles, kind=real64)
    source%Tj         = real(src_Tj,         kind=real64)
    source%dtau       = real(src_dtau,       kind=real64)

    ! Build COMET
    comet%coords  = real(comet_coords,  kind=real64)
    comet%Vastvec = real(comet_vastvec, kind=real64)
    comet%Vast    = real(comet_vast,    kind=real64)

    pericenter = (pericenter_c /= 0)

    if (env_enabled("HC_BRIDGE_DIAG")) then
       print *, "[py_bridge] DIAG enabled in py_hc_v_integration"
       print *, "  point: r, alpha, beta =", point%r, point%alpha, point%beta
       print *, "  point: rvector        =", point%rvector
       print *, "  source: r, alphaM, betaM =", source%r, source%alphaM, source%betaM
       print *, "  source: rrM, zeta, eta   =", source%rrM, source%zeta, source%eta
       print *, "  source: axis              =", source%symmetry_axis
       print *, "  source: ang_distr, ud_shape =", source%ejection_angle_distr, source%ud%ud_shape
       print *, "  source: umin, umax          =", source%ud%umin, source%ud%umax
       print *, "  source: Nparticles, Tj, dtau =", source%Nparticles, source%Tj, source%dtau
       ! comet/cloud depending on routine
       ! v_integration & delta_ejection:
       !   print *, "  comet: coords, Vastvec, Vast =", comet%coords, comet%Vastvec, comet%Vast
       ! simple_expansion:
       !   print *, "  cloudcentr, dt =", cloudcentr, dt
       density_sp = 0.0_real32
       ! return early, skip kernel call
       call flush(output_unit)   ! <<--- add this
        density = 0.0_c_double
        return
    end if
    
    call hc_DUDI_v_integration(density_sp, point, source, real(muR,kind=real64), real(tnow,kind=real64), &
                               comet, real(Rast_AU,kind=real64), pericenter)

    density = real(density_sp, kind=real64)
  end subroutine py_hc_v_integration


  subroutine py_hc_delta_ejection( &
        point_r, point_alpha, point_beta, point_rvector,             &
        src_r, src_alphaM, src_betaM, src_rrM, src_zeta, src_eta,    &
        src_axis, src_eject_distr, src_ud_shape, src_umin, src_umax, &
        src_Nparticles, src_Tj, src_dtau,                            &
        comet_coords, comet_vastvec, comet_vast,                     &
        muR, dt, Rast_AU, density) bind(C, name="py_hc_delta_ejection")
    use, intrinsic :: iso_fortran_env, only: real64, real32, output_unit

    real(c_double), value, intent(in) :: point_r, point_alpha, point_beta
    real(c_double), intent(in) :: point_rvector(3)

    real(c_double), value, intent(in) :: src_r, src_alphaM, src_betaM, src_zeta, src_eta
    real(c_double), intent(in) :: src_rrM(3), src_axis(3)
    integer(c_int), value, intent(in) :: src_eject_distr, src_ud_shape
    real(c_double), value, intent(in) :: src_umin, src_umax
    real(c_double), value, intent(in) :: src_Nparticles, src_Tj, src_dtau

    real(c_double), intent(in) :: comet_coords(3), comet_vastvec(3)
    real(c_double), value, intent(in) :: comet_vast

    real(c_double), value, intent(in) :: muR, dt, Rast_AU
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
    source%Nparticles = real(src_Nparticles, kind=real64)
    source%Tj         = real(src_Tj,         kind=real64)
    source%dtau       = real(src_dtau,       kind=real64)

    comet%coords  = real(comet_coords,  kind=real64)
    comet%Vastvec = real(comet_vastvec, kind=real64)
    comet%Vast    = real(comet_vast,    kind=real64)

    if (env_enabled("HC_BRIDGE_DIAG")) then
       print *, "[py_bridge] DIAG enabled in py_hc_delta_ejection"
       print *, "  point: r, alpha, beta =", point%r, point%alpha, point%beta
       print *, "  point: rvector        =", point%rvector
       print *, "  source: r, alphaM, betaM =", source%r, source%alphaM, source%betaM
       print *, "  source: rrM, zeta, eta   =", source%rrM, source%zeta, source%eta
       print *, "  source: axis              =", source%symmetry_axis
       print *, "  source: ang_distr, ud_shape =", source%ejection_angle_distr, source%ud%ud_shape
       print *, "  source: umin, umax          =", source%ud%umin, source%ud%umax
       print *, "  source: Nparticles, Tj, dtau =", source%Nparticles, source%Tj, source%dtau
       ! comet/cloud depending on routine
       ! v_integration & delta_ejection:
       !   print *, "  comet: coords, Vastvec, Vast =", comet%coords, comet%Vastvec, comet%Vast
       ! simple_expansion:
       !   print *, "  cloudcentr, dt =", cloudcentr, dt
       density_sp = 0.0_real32
       ! return early, skip kernel call
       call flush(output_unit)   ! <<--- add this
        density = 0.0_c_double
        return
    end if
    
    call hc_DUDI_delta_ejection(density_sp, point, source, real(muR,kind=real64), real(dt,kind=real64), &
                                comet, real(Rast_AU,kind=real64))

    density = real(density_sp, kind=real64)
  end subroutine py_hc_delta_ejection


  subroutine py_hc_simple_expansion( &
        point_r, point_alpha, point_beta, point_rvector,             &
        src_r, src_alphaM, src_betaM, src_rrM, src_zeta, src_eta,    &
        src_axis, src_eject_distr, src_ud_shape, src_umin, src_umax, &
        src_Nparticles, src_Tj, src_dtau,                            &
        cloudcentr, dt, density) bind(C, name="py_hc_simple_expansion")
    use, intrinsic :: iso_fortran_env, only: real64, real32, output_unit

    real(c_double), value, intent(in) :: point_r, point_alpha, point_beta
    real(c_double), intent(in) :: point_rvector(3)

    real(c_double), value, intent(in) :: src_r, src_alphaM, src_betaM, src_zeta, src_eta
    real(c_double), intent(in) :: src_rrM(3), src_axis(3)
    integer(c_int), value, intent(in) :: src_eject_distr, src_ud_shape
    real(c_double), value, intent(in) :: src_umin, src_umax
    real(c_double), value, intent(in) :: src_Nparticles, src_Tj, src_dtau

    real(c_double), intent(in) :: cloudcentr(3)
    real(c_double), value, intent(in) :: dt
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
    source%Nparticles = real(src_Nparticles, kind=real64)
    source%Tj         = real(src_Tj,         kind=real64)
    source%dtau       = real(src_dtau,       kind=real64)
 
    if (env_enabled("HC_BRIDGE_DIAG")) then
       print *, "[py_bridge] DIAG enabled in py_hc_simple_expansion"
       print *, "  point: r, alpha, beta =", point%r, point%alpha, point%beta
       print *, "  point: rvector        =", point%rvector
       print *, "  source: r, alphaM, betaM =", source%r, source%alphaM, source%betaM
       print *, "  source: rrM, zeta, eta   =", source%rrM, source%zeta, source%eta
       print *, "  source: axis              =", source%symmetry_axis
       print *, "  source: ang_distr, ud_shape =", source%ejection_angle_distr, source%ud%ud_shape
       print *, "  source: umin, umax          =", source%ud%umin, source%ud%umax
       print *, "  source: Nparticles, Tj, dtau =", source%Nparticles, source%Tj, source%dtau
       ! comet/cloud depending on routine
       ! v_integration & delta_ejection:
       !   print *, "  comet: coords, Vastvec, Vast =", comet%coords, comet%Vastvec, comet%Vast
       ! simple_expansion:
       !   print *, "  cloudcentr, dt =", cloudcentr, dt
       density_sp = 0.0_real32
       ! return early, skip kernel call
        call flush(output_unit)   ! <<--- add this
        density = 0.0_c_double
        return
    end if   
    call hc_DUDI_simple_expansion(density_sp, source, real(dt,kind=real64), real(cloudcentr,kind=real64), point)

    density = real(density_sp, kind=real64)
  end subroutine py_hc_simple_expansion


  !===================================================================
  !  Batched wrappers: over points (fixed source) and over sources
  !  (fixed point). These are C-interoperable entry points that
  !  reconstruct derived types and call hc_DUDI_batch_*.
  !
  !  NOTE:
  !    - Fortran LOGICAL is replaced by INTEGER(C_INT) at C boundary.
  !    - DUDIhc returns REAL (single); we convert to REAL(C_DOUBLE).
  !===================================================================

  subroutine py_hc_batch_points( &
       n_points, density,                                           &
       point_r, point_alpha, point_beta, point_rvector,             &
       src_r, src_alphaM, src_betaM, src_rrM, src_zeta, src_eta,    &
       src_axis, src_eject_distr, src_ud_shape, src_umin, src_umax, &
       src_Nparticles, src_Tj, src_dtau,                            &
       comet_coords, comet_vastvec, comet_vast,                     &
       muR, tnow, dt, Rast_AU, pericenter_c,                        &
       cloudcentr, method_id)                                       &
       bind(C, name="py_hc_batch_points")

    use, intrinsic :: iso_fortran_env, only: real64, real32, output_unit
    use iso_c_binding, only: c_int, c_double
    use define_types, only: position_in_space, source_properties, ejection_speed_properties, ephemeris
    use batching,      only: hc_DUDI_batch_points

    ! sizes
    integer(c_int), value, intent(in) :: n_points

    ! C-side outputs
    real(c_double), intent(out) :: density(n_points)

    ! C-side point arrays
    real(c_double), intent(in) :: point_r(n_points)
    real(c_double), intent(in) :: point_alpha(n_points)
    real(c_double), intent(in) :: point_beta(n_points)
    ! flattened or 2D; here we assume (3, n_points) layout
    real(c_double), intent(in) :: point_rvector(3, n_points)

    ! C-side source (scalars, same for all points)
    real(c_double), value, intent(in) :: src_r, src_alphaM, src_betaM, src_zeta, src_eta
    real(c_double),        intent(in) :: src_rrM(3), src_axis(3)
    integer(c_int), value, intent(in) :: src_eject_distr, src_ud_shape
    real(c_double), value, intent(in) :: src_umin, src_umax
    real(c_double), value, intent(in) :: src_Nparticles, src_Tj, src_dtau

    ! C-side comet (same for all points)
    real(c_double), intent(in) :: comet_coords(3), comet_vastvec(3)
    real(c_double), value, intent(in) :: comet_vast

    ! scalars / flags
    real(c_double), value, intent(in) :: muR, tnow, dt, Rast_AU
    integer(c_int), value, intent(in) :: pericenter_c
    real(c_double), intent(in) :: cloudcentr(3)
    integer(c_int), value, intent(in) :: method_id

    ! Local derived types
    type(position_in_space)        :: points(n_points)
    type(source_properties)        :: source
    type(ejection_speed_properties):: ud
    type(ephemeris)                :: comet
    logical                        :: pericenter
    integer                        :: i, method_f

    ! internal single-precision densities
    real(real32) :: density_sp(n_points)

    !--------------------------------------------
    ! diagnostic mode: just print and return zero
    !--------------------------------------------
    if (env_enabled("HC_BRIDGE_DIAG")) then
       print *, "[py_bridge] DIAG enabled in py_hc_batch_points"
       print *, "  n_points =", n_points
       print *, "  method_id =", method_id
       call flush(output_unit)
       density(:) = 0.0_c_double
       return
    end if

    !--------------------------------------------
    ! Build POINT array
    !--------------------------------------------
    do i = 1, n_points
       points(i)%r       = real(point_r(i),      kind=real64)
       points(i)%alpha   = real(point_alpha(i),  kind=real64)
       points(i)%beta    = real(point_beta(i),   kind=real64)
       points(i)%rvector = real(point_rvector(:, i), kind=real64)
    end do

    !--------------------------------------------
    ! Build SOURCE (same for all points)
    !--------------------------------------------
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
    source%Nparticles = real(src_Nparticles, kind=real64)
    source%Tj         = real(src_Tj,         kind=real64)
    source%dtau       = real(src_dtau,       kind=real64)

    !--------------------------------------------
    ! Build COMET (same for all points)
    !--------------------------------------------
    comet%coords  = real(comet_coords,  kind=real64)
    comet%Vastvec = real(comet_vastvec, kind=real64)
    comet%Vast    = real(comet_vast,    kind=real64)

    pericenter = (pericenter_c /= 0_c_int)
    method_f   = int(method_id, kind=kind(method_f))

    !--------------------------------------------
    ! Call batched kernel (single precision density)
    !--------------------------------------------
    call hc_DUDI_batch_points( n_points, density_sp, points, source, &
                               real(muR,    kind=real64),           &
                               real(tnow,   kind=real64),           &
                               real(dt,     kind=real64),           &
                               comet,                                &
                               real(Rast_AU, kind=real64),          &
                               pericenter,                           &
                               real(cloudcentr, kind=real64),       &
                               method_f )

    !--------------------------------------------
    ! Convert to REAL(C_DOUBLE) for C boundary
    !--------------------------------------------
    do i = 1, n_points
       density(i) = real(density_sp(i), kind=real64)
    end do

  end subroutine py_hc_batch_points


  subroutine py_hc_batch_sources( &
       n_sources, density,                                           &
       point_r, point_alpha, point_beta, point_rvector,              &
       src_r, src_alphaM, src_betaM, src_rrM, src_zeta, src_eta,     &
       src_axis, src_eject_distr, src_ud_shape, src_umin, src_umax,  &
       src_Nparticles, src_Tj, src_dtau,                             &
       comet_coords, comet_vastvec, comet_vast,                      &
       muR, tnow, dt, Rast_AU, pericenter_c,                         &
       cloudcentr, method_id)                                        &
       bind(C, name="py_hc_batch_sources")

    use, intrinsic :: iso_fortran_env, only: real64, real32, output_unit
    use iso_c_binding, only: c_int, c_double
    use define_types, only: position_in_space, source_properties, ejection_speed_properties, ephemeris
    use batching,      only: hc_DUDI_batch_sources

    ! sizes
    integer(c_int), value, intent(in) :: n_sources

    ! C-side outputs
    real(c_double), intent(out) :: density(n_sources)

    ! C-side point (single)
    real(c_double), value, intent(in) :: point_r, point_alpha, point_beta
    real(c_double),       intent(in) :: point_rvector(3)

    ! C-side source arrays (vary over sources)
    real(c_double), intent(in) :: src_r(n_sources)
    real(c_double), intent(in) :: src_alphaM(n_sources)
    real(c_double), intent(in) :: src_betaM(n_sources)
    real(c_double), intent(in) :: src_rrM(3, n_sources)
    real(c_double), intent(in) :: src_zeta(n_sources)
    real(c_double), intent(in) :: src_eta(n_sources)
    real(c_double), intent(in) :: src_axis(3, n_sources)
    integer(c_int), intent(in) :: src_eject_distr(n_sources)
    integer(c_int), intent(in) :: src_ud_shape(n_sources)
    real(c_double), intent(in) :: src_umin(n_sources)
    real(c_double), intent(in) :: src_umax(n_sources)
    real(c_double), intent(in) :: src_Nparticles(n_sources)
    real(c_double), intent(in) :: src_Tj(n_sources)
    real(c_double), intent(in) :: src_dtau(n_sources)

    ! C-side comet (same for all sources)
    real(c_double), intent(in) :: comet_coords(3), comet_vastvec(3)
    real(c_double), value, intent(in) :: comet_vast

    ! scalars / flags
    real(c_double), value, intent(in) :: muR, tnow, dt, Rast_AU
    integer(c_int), value, intent(in) :: pericenter_c
    real(c_double), intent(in) :: cloudcentr(3)
    integer(c_int), value, intent(in) :: method_id

    ! Local derived types
    type(position_in_space)         :: point
    type(source_properties)         :: sources(n_sources)
    type(ejection_speed_properties) :: ud
    type(ephemeris)                 :: comet
    logical                         :: pericenter
    integer                         :: i, method_f

    ! internal single-precision densities
    real(real32) :: density_sp(n_sources)

    !--------------------------------------------
    ! diagnostic mode: just print and return zero
    !--------------------------------------------
    if (env_enabled("HC_BRIDGE_DIAG")) then
       print *, "[py_bridge] DIAG enabled in py_hc_batch_sources"
       print *, "  n_sources =", n_sources
       print *, "  method_id =", method_id
       call flush(output_unit)
       density(:) = 0.0_c_double
       return
    end if

    !--------------------------------------------
    ! Build POINT (single)
    !--------------------------------------------
    point%r       = real(point_r,      kind=real64)
    point%alpha   = real(point_alpha,  kind=real64)
    point%beta    = real(point_beta,   kind=real64)
    point%rvector = real(point_rvector, kind=real64)

    !--------------------------------------------
    ! Build SOURCES array
    !--------------------------------------------
    do i = 1, n_sources
       sources(i)%r      = real(src_r(i),      kind=real64)
       sources(i)%alphaM = real(src_alphaM(i), kind=real64)
       sources(i)%betaM  = real(src_betaM(i),  kind=real64)
       sources(i)%rrM    = real(src_rrM(:, i), kind=real64)
       sources(i)%zeta   = real(src_zeta(i),   kind=real64)
       sources(i)%eta    = real(src_eta(i),    kind=real64)
       sources(i)%symmetry_axis        = real(src_axis(:, i), kind=real64)
       sources(i)%ejection_angle_distr = int(src_eject_distr(i), kind=kind(sources(i)%ejection_angle_distr))
       ud%ud_shape   = int(src_ud_shape(i), kind=kind(ud%ud_shape))
       ud%umin       = real(src_umin(i), kind=real64)
       ud%umax       = real(src_umax(i), kind=real64)
       sources(i)%ud  = ud
       sources(i)%Nparticles = real(src_Nparticles(i), kind=real64)
       sources(i)%Tj         = real(src_Tj(i),         kind=real64)
       sources(i)%dtau       = real(src_dtau(i),       kind=real64)
    end do

    !--------------------------------------------
    ! Build COMET (same for all sources)
    !--------------------------------------------
    comet%coords  = real(comet_coords,  kind=real64)
    comet%Vastvec = real(comet_vastvec, kind=real64)
    comet%Vast    = real(comet_vast,    kind=real64)

    pericenter = (pericenter_c /= 0_c_int)
    method_f   = int(method_id, kind=kind(method_f))

    !--------------------------------------------
    ! Call batched kernel
    !--------------------------------------------
    call hc_DUDI_batch_sources( n_sources, density_sp, point, sources, &
                                real(muR,    kind=real64),           &
                                real(tnow,   kind=real64),           &
                                real(dt,     kind=real64),           &
                                comet,                                &
                                real(Rast_AU, kind=real64),          &
                                pericenter,                           &
                                real(cloudcentr, kind=real64),       &
                                method_f )

    !--------------------------------------------
    ! Convert to REAL(C_DOUBLE) for C boundary
    !--------------------------------------------
    do i = 1, n_sources
       density(i) = real(density_sp(i), kind=real64)
    end do

  end subroutine py_hc_batch_sources
  
  
  !===================================================================
  !  C-interoperable wrapper for the general time×sources×points batch.
  !
  !  Layout expectations at C/Python side:
  !
  !   - n_points, Nt, Ns: sizes
  !   - point_*: length n_points
  !   - point_rvector: length 3*n_points   (flattened [x0,y0,z0,x1,y1,z1,...])
  !
  !   - src_*: length Nt*Ns, flattened with time-major order:
  !       idx = i_t * Ns + i_s   (0-based)
  !
  !   - src_rrM, src_axis: length 3*Nt*Ns, layout like:
  !       [x(t0,s0),y(t0,s0),z(t0,s0), x(t0,s1),..., x(t1,s0),...]
  !
  !   - comet_coords, comet_vastvec: length 3*Nt (coords per time)
  !   - comet_vast: length Nt
  !
  !===================================================================
  subroutine py_hc_batch_sources_points( &
       n_points, Nt, Ns, density,                    &
       point_r, point_alpha, point_beta,             &
       point_rvector,                               &
       src_r, src_alphaM, src_betaM,                &
       src_rrM, src_zeta, src_eta,                  &
       src_axis, src_eject_distr, src_ud_shape,     &
       src_umin, src_umax,                          &
       src_Nparticles, src_Tj, src_dtau,            &
       comet_coords, comet_vastvec, comet_vast,     &
       muR, tnow, Rast_AU, pericenter_c,            &
       method_id)                                    &
       bind(C, name="py_hc_batch_sources_points")

    use, intrinsic :: iso_fortran_env, only: real64
    use iso_c_binding, only: c_int, c_double
    use define_types, only: position_in_space, source_properties, ejection_speed_properties, ephemeris
    use batching,      only: hc_DUDI_batch_sources_points

    integer(c_int), value, intent(in) :: n_points, Nt, Ns
    real(c_double), intent(out) :: density(n_points)

    ! points
    real(c_double), intent(in) :: point_r(n_points)
    real(c_double), intent(in) :: point_alpha(n_points)
    real(c_double), intent(in) :: point_beta(n_points)
    real(c_double), intent(in) :: point_rvector(3*n_points)

    ! sources (flattened Nt*Ns and 3*Nt*Ns)
    real(c_double), intent(in) :: src_r(Nt*Ns)
    real(c_double), intent(in) :: src_alphaM(Nt*Ns)
    real(c_double), intent(in) :: src_betaM(Nt*Ns)
    real(c_double), intent(in) :: src_rrM(3*Nt*Ns)
    real(c_double), intent(in) :: src_zeta(Nt*Ns)
    real(c_double), intent(in) :: src_eta(Nt*Ns)
    real(c_double), intent(in) :: src_axis(3*Nt*Ns)
    integer(c_int), intent(in) :: src_eject_distr(Nt*Ns)
    integer(c_int), intent(in) :: src_ud_shape(Nt*Ns)
    real(c_double), intent(in) :: src_umin(Nt*Ns)
    real(c_double), intent(in) :: src_umax(Nt*Ns)
    real(c_double), intent(in) :: src_Nparticles(Nt*Ns)
    real(c_double), intent(in) :: src_Tj(Nt*Ns)
    real(c_double), intent(in) :: src_dtau(Nt*Ns)

    ! comets: coords, vvec, Vast per time
    real(c_double), intent(in) :: comet_coords(3*Nt)
    real(c_double), intent(in) :: comet_vastvec(3*Nt)
    real(c_double), intent(in) :: comet_vast(Nt)

    ! scalars
    real(c_double), value, intent(in) :: muR, tnow, Rast_AU
    integer(c_int), value, intent(in) :: pericenter_c
    integer(c_int), value, intent(in) :: method_id

    ! local derived types
    type(position_in_space)         :: points(n_points)
    type(source_properties)         :: sources(Nt, Ns)
    type(ejection_speed_properties) :: ud
    type(ephemeris)                 :: comets(Nt)
    logical                         :: pericenter
    integer                         :: i, i_t, i_s, idx
    integer                         :: method_f
    integer                         :: k0

    real(real32) :: dens_sp(n_points)

    !----- build points -----
    do i = 1, n_points
       points(i)%r     = real(point_r(i),     kind=real64)
       points(i)%alpha = real(point_alpha(i), kind=real64)
       points(i)%beta  = real(point_beta(i),  kind=real64)
       k0 = 3*(i-1)
       points(i)%rvector(1) = real(point_rvector(k0+1), kind=real64)
       points(i)%rvector(2) = real(point_rvector(k0+2), kind=real64)
       points(i)%rvector(3) = real(point_rvector(k0+3), kind=real64)
    end do

    !----- build sources(Nt,Ns) -----
    do i_t = 1, Nt
       do i_s = 1, Ns
          idx = (i_t-1)*Ns + i_s   ! matches C-order flattening
          sources(i_t, i_s)%r      = real(src_r(idx),      kind=real64)
          sources(i_t, i_s)%alphaM = real(src_alphaM(idx), kind=real64)
          sources(i_t, i_s)%betaM  = real(src_betaM(idx),  kind=real64)

          k0 = 3*(idx-1)
          sources(i_t, i_s)%rrM(1) = real(src_rrM(k0+1), kind=real64)
          sources(i_t, i_s)%rrM(2) = real(src_rrM(k0+2), kind=real64)
          sources(i_t, i_s)%rrM(3) = real(src_rrM(k0+3), kind=real64)

          sources(i_t, i_s)%zeta = real(src_zeta(idx), kind=real64)
          sources(i_t, i_s)%eta  = real(src_eta(idx),  kind=real64)

          sources(i_t, i_s)%symmetry_axis(1) = real(src_axis(k0+1), kind=real64)
          sources(i_t, i_s)%symmetry_axis(2) = real(src_axis(k0+2), kind=real64)
          sources(i_t, i_s)%symmetry_axis(3) = real(src_axis(k0+3), kind=real64)

          sources(i_t, i_s)%ejection_angle_distr = int(src_eject_distr(idx), kind=kind(sources(i_t,i_s)%ejection_angle_distr))

          ud%ud_shape = int(src_ud_shape(idx), kind=kind(ud%ud_shape))
          ud%umin     = real(src_umin(idx), kind=real64)
          ud%umax     = real(src_umax(idx), kind=real64)
          sources(i_t, i_s)%ud = ud

          sources(i_t, i_s)%Nparticles = real(src_Nparticles(idx), kind=real64)
          sources(i_t, i_s)%Tj         = real(src_Tj(idx),         kind=real64)
          sources(i_t, i_s)%dtau       = real(src_dtau(idx),       kind=real64)
       end do
    end do

    !----- build comets(Nt) -----
    do i_t = 1, Nt
       k0 = 3*(i_t-1)
       comets(i_t)%coords(1)  = real(comet_coords(k0+1),  kind=real64)
       comets(i_t)%coords(2)  = real(comet_coords(k0+2),  kind=real64)
       comets(i_t)%coords(3)  = real(comet_coords(k0+3),  kind=real64)
       comets(i_t)%Vastvec(1) = real(comet_vastvec(k0+1), kind=real64)
       comets(i_t)%Vastvec(2) = real(comet_vastvec(k0+2), kind=real64)
       comets(i_t)%Vastvec(3) = real(comet_vastvec(k0+3), kind=real64)
       comets(i_t)%Vast       = real(comet_vast(i_t),     kind=real64)
    end do

    pericenter = (pericenter_c /= 0_c_int)
    method_f   = int(method_id, kind=kind(method_f))

    call hc_DUDI_batch_sources_points( n_points, Nt, Ns, dens_sp, points, sources, &
                                       real(muR,   kind=real64), &
                                       real(tnow,  kind=real64), &
                                       comets,                    &
                                       real(Rast_AU, kind=real64), &
                                       pericenter, method_f )

    do i = 1, n_points
       density(i) = real(dens_sp(i), kind=real64)
    end do

  end subroutine py_hc_batch_sources_points


	 !--- tabulated ejection speed distribution ---------------------------------

	  subroutine py_set_tabulated_fu(nu_in, u_in, fu_in) &
					   bind(C, name="py_set_tabulated_fu")
		use iso_c_binding, only: c_int, c_double
		implicit none
		integer(c_int), value, intent(in) :: nu_in
		real(c_double), intent(in) :: u_in(nu_in)
		real(c_double), intent(in) :: fu_in(nu_in)

		call set_tabulated_fu(int(nu_in, kind=kind(1)), u_in, fu_in)
	  end subroutine py_set_tabulated_fu


	  !--- tabulated ejection direction distribution -----------------------------

	  subroutine py_set_tabulated_fpsi(Npsi_in, NlambdaM_in, psi_in, lambdaM_in, fpsi_in) &
					   bind(C, name="py_set_tabulated_fpsi")
		use iso_c_binding, only: c_int, c_double
		implicit none
		integer(c_int), value, intent(in) :: Npsi_in, NlambdaM_in
		real(c_double), intent(in) :: psi_in(Npsi_in)
		real(c_double), intent(in) :: lambdaM_in(NlambdaM_in)
		real(c_double), intent(in) :: fpsi_in(Npsi_in, NlambdaM_in)

		call set_tabulated_fpsi( int(Npsi_in,     kind=kind(1)), &
								 int(NlambdaM_in, kind=kind(1)), &
								 psi_in, lambdaM_in, fpsi_in )
	  end subroutine py_set_tabulated_fpsi
  !--- dimensions and limits --------------------------------------------------

  subroutine get_ratemap_dims(nlats_out, nlons_out) & 
                   bind(C, name="py_get_ratemap_dims")
    integer, intent(out) :: nlats_out, nlons_out
    nlats_out = nlats
    nlons_out = nlons
  end subroutine get_ratemap_dims


  subroutine get_lon_limits(lonmin_out, lonmax_out) & 
                   bind(C, name="py_get_lon_limits")
    real(8), intent(out) :: lonmin_out, lonmax_out
    lonmin_out = lonmin
    lonmax_out = lonmax
  end subroutine get_lon_limits


  subroutine set_lon_limits(lonmin_in, lonmax_in) & 
                   bind(C, name="py_set_lon_limits")
    real(8), intent(in) :: lonmin_in, lonmax_in
    lonmin = lonmin_in
    lonmax = lonmax_in
  end subroutine set_lon_limits


  !--- latitude / longitude grids ---------------------------------------------

  subroutine get_lats(lats_out) & 
                   bind(C, name="py_get_lats")
    real, intent(out) :: lats_out(nlats)
    lats_out = lats
  end subroutine get_lats


  subroutine get_lons(lons_out) &
                   bind(C, name="py_get_lons")
    real, intent(out) :: lons_out(nlons)
    lons_out = lons
  end subroutine get_lons


  subroutine set_lats(lats_in) &
                   bind(C, name="py_set_lats")
    real, intent(in) :: lats_in(nlats)
    lats = lats_in
  end subroutine set_lats


  subroutine set_lons(lons_in) &
                   bind(C, name="py_set_lons")
    real, intent(in) :: lons_in(nlons)
    lons = lons_in
  end subroutine set_lons


  !--- ratemap matrices -------------------------------------------------------

  subroutine get_ratemap(ratemap_out) &
                   bind(C, name="py_get_ratemap")
    real(8), intent(out) :: ratemap_out(nlons, nlats)
    ratemap_out = ratemap
  end subroutine get_ratemap


  subroutine get_rmap1(rmap1_out) &
                   bind(C, name="py_get_rmap1")
    real(8), intent(out) :: rmap1_out(nlons, nlats)
    rmap1_out = rmap1
  end subroutine get_rmap1


  subroutine get_rmap2(rmap2_out) &
                   bind(C, name="py_get_rmap2")
    real(8), intent(out) :: rmap2_out(nlons, nlats)
    rmap2_out = rmap2
  end subroutine get_rmap2


  subroutine set_ratemap(ratemap_in) &
                   bind(C, name="py_set_ratemap")
    real(8), intent(in) :: ratemap_in(nlons, nlats)
    ratemap = ratemap_in
  end subroutine set_ratemap


  subroutine set_rmap1(rmap1_in) &
                   bind(C, name="py_set_rmap1")
    real(8), intent(in) :: rmap1_in(nlons, nlats)
    rmap1 = rmap1_in
  end subroutine set_rmap1


  subroutine set_rmap2(rmap2_in) &
                   bind(C, name="py_set_rmap2")
    real(8), intent(in) :: rmap2_in(nlons, nlats)
    rmap2 = rmap2_in
  end subroutine set_rmap2


  !--- temporary vector rMtmp -------------------------------------------------

  subroutine get_rMtmp(rMtmp_out) &
                   bind(C, name="py_get_rMtmp")
    real(8), intent(out) :: rMtmp_out(3)
    rMtmp_out = rMtmp
  end subroutine get_rMtmp


  subroutine set_rMtmp(rMtmp_in) &
                   bind(C, name="py_set_rMtmp")
    real(8), intent(in) :: rMtmp_in(3)
    rMtmp = rMtmp_in
  end subroutine set_rMtmp
  


  ! Wrapper for:
  !   subroutine read_ratemap(fname, rhel)
  !     character(*), intent(in) :: fname
  !     real(8),      intent(out):: rhel
  !
  ! C interface:
  !   void py_read_ratemap(const char *fname, double *rhel);
  !
  subroutine read_ratemap_get_rhel(fname_c, rhel) bind(C, name="py_read_ratemap_get_rhel")
    use iso_c_binding
    implicit none
    character(kind=c_char), intent(in) :: fname_c(*)   ! C string (null-terminated)
    real(c_double),          intent(out) :: rhel       ! C double* (by reference)

    character(len=512) :: fname_f
    real(8)            :: rhel_f
    integer            :: i

    ! Convert C null-terminated string to Fortran CHARACTER(*)
    fname_f = ' '
    do i = 1, len(fname_f)
       if (fname_c(i) == c_null_char) exit
       fname_f(i:i) = achar(iachar(fname_c(i)))
    end do

    call read_ratemap(trim(fname_f), rhel_f)
    rhel = rhel_f
  end subroutine read_ratemap_get_rhel



  ! Wrapper for:
  !   subroutine ratematr_interpolate(rhel, rhel1, rhel2)
  !     real(8), intent(in) :: rhel, rhel1, rhel2
  !
  ! C interface:
  !   void py_ratematr_interpolate(double rhel,
  !                                double rhel1,
  !                                double rhel2);
  !
  subroutine py_ratematr_interpolate(rhel, rhel1, rhel2) &
       bind(C, name="py_ratematr_interpolate")
    use iso_c_binding
    implicit none
    real(c_double), value :: rhel, rhel1, rhel2   ! C doubles passed by value

    real(8) :: rhel_f, rhel1_f, rhel2_f

    rhel_f  = real(rhel,  kind=8)
    rhel1_f = real(rhel1, kind=8)
    rhel2_f = real(rhel2, kind=8)

    call ratematr_interpolate(rhel_f, rhel1_f, rhel2_f)
  end subroutine py_ratematr_interpolate


end module py_dudihc_bridge
