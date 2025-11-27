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
  use batching, only: hc_DUDI_batch_points, hc_DUDI_batch_sources
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
                         lats, lons, rmap1, rmap2, ratemap, rMtmp
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
  
      !===================== getters for sizes =====================
    integer(c_int) function py_get_nlats() bind(C, name="py_get_nlats")
      use, intrinsic :: iso_c_binding, only: c_int
      py_get_nlats = nlats
    end function

    integer(c_int) function py_get_nlons() bind(C, name="py_get_nlons")
      use, intrinsic :: iso_c_binding, only: c_int
      py_get_nlons = nlons
    end function

    !===================== lon bounds ============================
    subroutine py_set_lon_bounds(lonmin_in, lonmax_in) bind(C, name="py_set_lon_bounds")
      use, intrinsic :: iso_c_binding, only: c_double
      real(c_double), value :: lonmin_in, lonmax_in
      lonmin = lonmin_in
      lonmax = lonmax_in
    end subroutine

    subroutine py_get_lon_bounds(lonmin_out, lonmax_out) bind(C, name="py_get_lon_bounds")
      use, intrinsic :: iso_c_binding, only: c_double
      real(c_double) :: lonmin_out, lonmax_out
      lonmin_out = lonmin
      lonmax_out = lonmax
    end subroutine

    !===================== 1D arrays (REAL(4)) ===================
    subroutine py_set_lats(n, arr) bind(C, name="py_set_lats")
      use, intrinsic :: iso_c_binding, only: c_int, c_float
      integer(c_int), value :: n
      real(c_float)         :: arr(n)
      if (n /= nlats) return
      lats(1:n) = arr(1:n)
    end subroutine

    subroutine py_set_lons(n, arr) bind(C, name="py_set_lons")
      use, intrinsic :: iso_c_binding, only: c_int, c_float
      integer(c_int), value :: n
      real(c_float)         :: arr(n)
      if (n /= nlons) return
      lons(1:n) = arr(1:n)
    end subroutine

    ! Optional getters (handy for tests/validation)
    subroutine py_get_lats(n, arr) bind(C, name="py_get_lats")
      use, intrinsic :: iso_c_binding, only: c_int, c_float
      integer(c_int), value :: n
      real(c_float)         :: arr(n)
      integer               :: k, m
      m = min(n, nlats)
      do k = 1, m
        arr(k) = lats(k)
      end do
    end subroutine

    subroutine py_get_lons(n, arr) bind(C, name="py_get_lons")
      use, intrinsic :: iso_c_binding, only: c_int, c_float
      integer(c_int), value :: n
      real(c_float)         :: arr(n)
      integer               :: k, m
      m = min(n, nlons)
      do k = 1, m
        arr(k) = lons(k)
      end do
    end subroutine

    !===================== 2D maps (REAL(8)) =====================
    ! NOTE: arrays are (nlons, nlats) in Fortran column-major.
    subroutine py_set_rmap1(nx, ny, A) bind(C, name="py_set_rmap1")
      use, intrinsic :: iso_c_binding, only: c_int, c_double
      integer(c_int), value :: nx, ny
      real(c_double)        :: A(nx, ny)
      if (nx==nlons .and. ny==nlats) rmap1(:,:) = A(:,:)
    end subroutine

    subroutine py_set_rmap2(nx, ny, A) bind(C, name="py_set_rmap2")
      use, intrinsic :: iso_c_binding, only: c_int, c_double
      integer(c_int), value :: nx, ny
      real(c_double)        :: A(nx, ny)
      if (nx==nlons .and. ny==nlats) rmap2(:,:) = A(:,:)
    end subroutine

    subroutine py_set_ratemap(nx, ny, A) bind(C, name="py_set_ratemap")
      use, intrinsic :: iso_c_binding, only: c_int, c_double
      integer(c_int), value :: nx, ny
      real(c_double)        :: A(nx, ny)
      if (nx==nlons .and. ny==nlats) ratemap(:,:) = A(:,:)
    end subroutine

    !===================== rMtmp (REAL(8), len=3) ================
    subroutine py_set_rmtmp(v) bind(C, name="py_set_rmtmp")
      use, intrinsic :: iso_c_binding, only: c_double
      real(c_double) :: v(3)
      rMtmp(1:3) = v(1:3)
    end subroutine

    subroutine py_get_rmtmp(v) bind(C, name="py_get_rmtmp")
      use, intrinsic :: iso_c_binding, only: c_double
      real(c_double) :: v(3)
      v(1:3) = rMtmp(1:3)
    end subroutine


end module py_dudihc_bridge
