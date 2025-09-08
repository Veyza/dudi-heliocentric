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

! File: twobody_fun.f90
! Description: Contains subroutines that calculate integrands
! for the delta-ejection and v-integration methods, along with auxiliary
! subroutines that support these calculations.


! Author: Anastasiia Ershova
! E-mail: vveyzaa@gmail.com

module twobody_fun
    implicit none
        contains
            
            ! comet velocity in system \tilde{S}
            subroutine find_astapex(rrM, r, VVast, Vast, VVVast)
                use help
                implicit none
                real(8), intent(in) :: r, rrM(3), VVast(3), Vast
                real(8), intent(out) :: VVVast(3)
                real(8) zetaast, etaast
                
                ! comet apex in the local horizontal SC
                zetaast = acos(dot_product(VVast, rrM) / r / Vast)
                etaast = azimuth(VVast, rrM)
                
                ! comet velocity vector in the local horizontal CS
                VVVast(1) = sin(zetaast) * cos(etaast)
                VVVast(2) = sin(zetaast) * sin(etaast)
                VVVast(3) = cos(zetaast)
                VVVast = VVVast * Vast
                
            end subroutine find_astapex
            



            function hc_jacobian_ueject(uej, usinpsieject, u, usinpsi) result(Jmotion)
                use const
                use help
                implicit none
                real(8) uej, u, usinpsi, usinpsieject
                real(8) Jmotion
                
                Jmotion = u * usinpsi / uej / usinpsieject
            
            end function hc_jacobian_ueject



            ! tests if the particle  re-impacted the asteroid on its way
            ! from the source to the point of interest
            subroutine collision_check(muR, dt, coords, s, point, &
                       uejectvec, collision, Rast_AU)
                use const
                use define_types
                use help
                implicit none
                type(position_in_space), intent(in) :: point
                type(source_properties), intent(in) :: s
                integer i, ind
                real(8), intent(in) :: muR, dt, Rast_AU
                real(8), intent(in) :: coords(3), uejectvec(3)
                real(8) rdif(3), polan, az, radast
                real(8) da(3), dv(3), dx(3), rero, imro, tmpa(3), dd(3)
                real(8) dda(3), ddda(3)
                real(8) coefs(0:3), rtmp, tmp
                complex(8) roots(3), ctmp
                logical, intent(out) :: collision
                
                collision = .False.
                ! if the grain was ejected toward the Sun
                if(uejectvec(3) < 0d0) then
                ! comet heliocentric distance at the moment of ejection
                    radast = norma3d(coords) 
                ! vector from the source position to the comet center 
                    rdif =  s%rrM - coords
                ! only direction of this vector  
                    rdif = rdif / Rast_AU
                ! polar angle and longitude of the comet center in local 
                ! (defined at the source position) horizontal CS
                    polan = acos(dot_product(rdif, s%rrM) / s%r)
                    az = azimuth(rdif, s%rrM)
                ! vector from the source position to the comet center
                ! in local horizontal CS
                    dx = Rast_AU * (/sin(polan) * cos(az), sin(polan) * sin(az), cos(polan)/)
                ! angle between heliocentric vectors
                ! of the source and the comet center
                    tmp = dot_product(s%rrM, -coords) / s%r / radast
                    if(abs(tmp) > 1d0) tmp = sign(1d0, tmp)
                    polan = acos(tmp)
                    az = azimuth(-coords/radast, s%rrM)
                ! acceleration acting at the comet in the horizontal CS
                ! (z pointing along the source heliocentric vector)
                    tmpa = -GMsun / s%r**2 * (/0d0, 0d0, 1d0/)
                ! differences
                    da = (/0d0, 0d0, -muR / s%r**2/) - tmpa
                    dv = uejectvec
                    dda = (/-dv(1), -dv(2), dv(3) * 2d0/) * muR / s%r**3

                    coefs(0) = dot_product(da, da) + 4d0 * dot_product(dv, dda) / 3d0
                    coefs(1) = 4d0 * dot_product(dv, da) + 4d0 * dot_product(dx, dda) / 3d0
                    coefs(2) = 4d0 * (dot_product(dv,dv) + dot_product(dx, da))
                    coefs(3) = 8d0 * dot_product(dx, dv)

                    call cardano_formula(coefs, roots)
                    
                    ! if a positive real root < dt exists
                    do i = 1, 3
                        rero = realpart(roots(i))
                        imro = imagpart(roots(i))
                        dd = dx + dv * rero + da / 2d0 * rero**2
                        rtmp = dot_product(dd, dd) - dot_product(dx, dx) &
                        - (dot_product(dv, dv) + dot_product(dx, da)) * rero**2 &
                        - 2d0 * dot_product(dv, dx) * rero &
                        - dot_product(dv, da) * rero**3 &
                        - dot_product(da, da) / 4d0 * rero**4
                        
                        collision = (abs(imro/ rero) < 1e-16 .and. &
                                    rero < dt .and. rero > 0.0)
                        if(collision) then
                            exit
                        endif
                    enddo
                    if(collision) then
                        dd = dx + dv * roots(ind) + da / 2d0 * roots(ind)**2
                        ctmp = dot_product(dd, dd) - dot_product(dx, dx) &
                        - (dot_product(dv, dv) + dot_product(dx, da)) * roots(ind)**2 &
                        - 2d0 * dot_product(dv, dx) * roots(ind) &
                        - dot_product(dv, da) * roots(ind)**3 &
                        - dot_product(da, da) / 4d0 * roots(ind)**4
                        if(abs(ctmp) > 1d-8) then
                            collision = .False.
                        endif
                    endif
                endif
                
            end subroutine collision_check




            ! compute \Delta\phi from \theta, velocity and spacecraft position
            function deltaphi(theta, r, r0, v, muR, pericenter)
                use const
                implicit none
                real(8), intent(in) :: theta
                real(8), intent(in) :: r, r0, v, muR
                real(8) deltaphi, e, pp, cosphi, cosphim
                real(8) v2, Ekep, h, u, sinpsi
                logical pericenter
                
                v2 = v * v
                h = r * v * sin(theta)
                Ekep = v2 / 2d0 - muR / r
                pp = h * h / muR
                e = sqrt(1d0 + 2d0 * pp * Ekep / muR)
                u = sqrt2d0 * sqrt(Ekep + muR / r0)
                sinpsi = h / r0 / u
                
                if(abs(muR) > 1d-9) then
                    if(muR > 0) then
                        cosphim = (pp / r0 - 1d0) / e
                        cosphi = (pp / r - 1d0) / e
                    else
                        cosphim = (1d0 - pp / r0) / e
                        cosphi = (1d0 - pp / r) / e
                    endif
                else
                    cosphim = sign(1d0, muR) * sinpsi
                    cosphi = sign(1d0, muR) * sin(theta)
                endif
                if(cosphim > 1d0) then
                    cosphim = sign(0.999995d0, cosphim)
                endif
                if(cosphi > 1d0) then
                    cosphi = sign(0.999995d0, cosphi)
                endif
                if(pericenter) then
                    deltaphi = acos(cosphi) + acos(cosphim)
                else
                    deltaphi = acos(cosphi) - acos(cosphim)
                endif

            end function deltaphi



            ! derivative of \Delta \phi by velocity v
            function derdphidv(rr, vv, theta, mu, rrm, pericenter)
                use const
                implicit none
                real(8), intent(in) :: rr, vv, theta, mu, rrm
                real(8) r, v, sintheta2, bracket, bracket2
                real(8) rv, rv2, r2v2, cos2theta, mu2, v2, rv2sin, r2v2sin
                real(8) derdphidv
                real(8) dphi1, dphi2, dphi3, dphi4, delta, numder
                real(8) Ekep
                logical pericenter

                sintheta2 = sin(theta)**2
                
                if(mu < 0.0) then
                    
                    rv = rr * vv
                    rv2 = rv * vv
                    r2v2 = rr * rv2
                    v2 = vv * vv
                    Ekep = v2 / 2d0 - mu / rr
                    mu2 = mu * mu
                    rv2sin = rv2 * sintheta2
                    r2v2sin = rr * rv2sin
                    bracket = 2d0 * r2v2 * Ekep * sintheta2
                    bracket2 = mu + rrm * Ekep
                    
                    derdphidv = (rv2sin * (2d0 * (mu2 + rr * Ekep * mu) &
                    + bracket + rv2 * (mu - rv2sin)) &
                    / (sqrt((rv2sin * (2d0 * (mu + rr * Ekep) - rv2sin)) &
                    / (mu2 + bracket)) * mu * (mu2 + bracket)) &
                    - (rrm * sqrt((2d0 * r2v2 * rrm * bracket2 * sintheta2 &
                    - r2v2sin**2) / (rrm**2 * (mu2 + bracket))) &
                    * (2d0 * (mu2 + Ekep * mu * rrm) &
                    + bracket + (mu * rrm - r2v2sin) * v2) &
                    / (mu * (2d0 * rrm * bracket2 - r2v2sin)))) &
                    / (vv * sqrt(1d0 + bracket / mu2))
                else
                    r = rr / rrm
                    v = vv / sqrt(2.0 * abs(mu) / rrm)
                    rv = r * v 
                    r2v2 = rv**2
                    rv2 = rv * v
                    
                    bracket = (-1d0 + rv2)
                    bracket2 = (1d0 + 4d0 * rv2 * bracket * sintheta2)
                    derdphidv = (2d0 * rv * sintheta2 * (-1d0 + r &
                              + 2d0 * rv2 - 2d0 * r2v2 * sintheta2)) &
                              / (sqrt((rv2 * sintheta2 * (-1d0 + r &
                              + rv2 - r * r2v2 * sintheta2)) /bracket2) &
                              * (1d0 + 4d0 * rv2 * bracket * sintheta2)**1.5) &
                              - (2d0 * sqrt((r2v2 * v**2 * sin(2d0 * theta)**2) &
                              / bracket2)) / (v * sqrt(bracket2))
                    derdphidv = derdphidv / sqrt2d0 / sqrt(abs(mu) / rrm)
                    

                endif
                    
                if(derdphidv /= derdphidv) then
                    delta = 1d-4
                    dphi1 = deltaphi(theta, rr, rrm, vv+2d0*delta, mu, pericenter)
                    dphi2 = deltaphi(theta, rr, rrm, vv+delta, mu, pericenter)
                    dphi3 = deltaphi(theta, rr, rrm, vv-delta, mu, pericenter)
                    dphi4 = deltaphi(theta, rr, rrm, vv-2d0*delta, mu, pericenter)
                    numder = (-dphi1 + 8d0 * dphi2 - 8d0 * dphi3 + dphi4) / 12d0 / delta 
                    write(666,*) 'the derivative d\Delta\phi/d\theta was obtained numerically &
                                because the analytical expression contains numerically difficult parts'
                    N_of_warnings = N_of_warnings + 1
                    if(N_of_warnings > maxNofWarnings) then
                        write(*,*) 'too many warnings have been printed to fort.666'
                        stop
                    endif
                    derdphidv = numder
                endif
                
            end function derdphidv
                    
            
            ! derivative of the angle \Delta\phi by the angle \theta
            function derdphidtheta(rr, vv, theta, mu, rrm, &
                                    pericenter, close2pericenter)
                use const
                implicit none
                real(8), intent(in) :: rr, vv, theta, mu, rrm
                real(8) r, v, sintheta2, bracket, bracket2, costhx2
                real(8) rv, rv2, r2v2, sin2theta, sintheta, Ekep, mu2
                real(8) derdphidtheta, r2v2sin
                real(8) dphi1, dphi2, dphi3, dphi4, delta, numder
                logical close2pericenter, pericenter
                
                if(.not. close2pericenter) then
                    r = rr / rrm
                    v = vv / sqrt(2.0 * abs(mu) / rrm)
                    sintheta = sin(theta)
                    sintheta2 = sintheta**2
                    sin2theta = sin(2d0 * theta)
                    
                    if(mu > 0.0) then
                        rv = r * v 
                        r2v2 = rv**2
                        rv2 = rv * v
                        bracket = -1.0 + rv2
                        bracket2 = 1.0 + 4.0 * rv2 * bracket * sintheta2

                        derdphidtheta = (rv2 * (-1.0 + r + rv2 + 2.0 * r * rv2 &
                                        * bracket * sintheta2) * 2d0 * cos(theta)) / &
                                        (sqrt((rv2 * (-1.0 + r + rv2 - r**3 * v**2 * sintheta2)) / &
                                        bracket2) * bracket2**1.5) &
                                        - (sign(1d0, sin2theta) * 2.0 * (1.0 + 2.0 * bracket * sintheta2) &
                                        * sqrt((rv2**2)/bracket2)) / &
                                        sqrt(bracket2)
                    else
                        Ekep = vv * vv / 2d0 - mu / rr
                        rv = rr * vv
                        rv2 = rv * vv
                        r2v2 = rr * rv2
                        mu2 = mu * mu
                        costhx2 = 2d0 * cos(theta)
                        bracket = mu2 + 2d0 * r2v2 * Ekep * sintheta2
                        bracket2 = rrm * (mu + rrm * Ekep)
                        r2v2sin = r2v2 * sintheta2
                        
                        derdphidtheta = ((rv2 * sintheta2 * (mu2 * costhx2 &
                        + rr * Ekep * (mu * costhx2 + rv2 * sintheta2 * costhx2))) &
                        / (sqrt((rv2 * sintheta2 * (2d0 * (mu + rr * Ekep)  &
                        - rv2 * sintheta2)) / bracket) * mu * bracket) &
                        - (rrm * sqrt((2d0 * bracket2 * r2v2sin - r2v2sin**2) &
                        / (rrm**2 * bracket)) &
                        * (costhx2 * (mu2 + Ekep * (mu * rrm + r2v2sin)))) &
                        / mu / (2d0 * bracket2 - r2v2 * sintheta2)) &
                        / sqrt(1d0 + (2d0 * Ekep * r2v2sin) / mu2) / sin(theta)
                                      

                    endif
                endif
                if(derdphidtheta /= derdphidtheta .or. close2pericenter .or. pericenter) then
                    if(close2pericenter) then
                        delta = 3d-2
                    else
                        delta = 1d-3
                    endif
                    dphi1 = deltaphi(theta+2d0*delta, rr, rrm, vv, mu, pericenter)
                    dphi2 = deltaphi(theta+delta, rr, rrm, vv, mu, pericenter)
                    dphi3 = deltaphi(theta-delta, rr, rrm, vv, mu, pericenter)
                    dphi4 = deltaphi(theta-2d0*delta, rr, rrm, vv, mu, pericenter)
                    numder = (-dphi1 + 8d0 * dphi2 - 8d0 * dphi3 + dphi4) / 12d0 / delta 
                    !write(*,*) dphi1, dphi2, dphi3, dphi4
                    !write(*,*) 'numder1', numder
                    !numder = (dphi2 - dphi3) / 2d0 / delta 
                    !write(*,*) 'numder2', numder
                    if((.not. close2pericenter) .and. (.not. pericenter)) then
                        write(666,*) 'the derivative d\Delta\phi/d\theta was obtained numerically &
                                because the analytical expression contains numerically difficult parts'
                        N_of_warnings = N_of_warnings + 1
                        if(N_of_warnings > maxNofWarnings) then
                            write(*,*) 'too many warnings have been printed to fort.666'
                            stop
                        endif
                    endif
                    derdphidtheta = numder
                endif
            end function derdphidtheta
                
            
            



            ! Evaluates the integrand in the case of \delta-ejection method
            ! (only one value of the integrand is computed for each pair 
            ! 'a source + a point of interest'
            subroutine Integrand_delta_ejection(Integrand, uu, velocity, &
                                          ee, psi, theta, point, &
                                         dphi, dbeta, s, muR, dt, &
                                         comet, Rast_AU)
                use const
                use define_types
                use help
                use distributions_fun
                implicit none
                integer i
                real(8), intent(in) :: uu, dphi, dbeta, psi
                real(8), intent(in) :: velocity, theta, ee
                real(8), intent(in) :: muR, dt, Rast_AU
                type(ephemeris), intent(in) :: comet
                type(source_properties), intent(in) :: s
                type(position_in_space), intent(in) :: point
                real(8) lambdaM, wpsi, jetdir(3)
                real, intent(out) :: Integrand
                real(8) ddphidtheta, ddphidv
                real(8) fac1, fac2, fac3, Jpsi
                real(8) sindphi, hcJpsi, Vastvechorz(3)
                real(8) ueject, uejectvec(3), psieject, lambdaMeject, uvec(3)
                real(8) ddeltatdv, ddeltatdtheta, tmp, vvec(3)
                real(8) zvec(3), xvec(3), yvec(3)
                logical collision
                
                Integrand = 0d0
                collision = .FALSE.
                wpsi = halfpi 
                ! azimuth of the particle velocity vector at the moment of ejection
                lambdaM = azimuth(point%rvector - s%rrM, s%rrM)
                
                ! particle velocity vector at the moment of ejection
                uvec = uu * (/sin(psi) * cos(lambdaM), sin(psi) * sin(lambdaM), cos(psi)/)
                ! particle velocity at the point of interest expressed in the horizontal CS defined in the point of ejection
                vvec = velocity * (/sin(theta) * cos(lambdaM), sin(theta) * sin(lambdaM), cos(theta)/)
                
                ! ejection velocity vector in the local horizontal CS
                if(comet%Vast /= 0d0) then
                    call find_astapex(s%rrM, s%r, comet%Vastvec, comet%Vast, Vastvechorz)
                    uejectvec = uvec - Vastvechorz
                else
                    Vastvechorz = 0d0
                    uejectvec = uvec
                endif
                ueject = norma3d(uejectvec)
                
                if(Rast_AU > 0d0) then
                    call collision_check(muR, dt, comet%coords, s, &
                            point, uejectvec, collision, Rast_AU)
                else
                    collision = .FALSE.
                endif
                       
                if((.not. collision) .and. &
                    ueject > s%ud%umin .and. ueject < s%ud%umax) then
                    fac1 = ejection_speed_distribution(s%ud, ueject)
                    fac1 = velocity * fac1 / uu / uu
                
                    psieject = acos(uejectvec(3) / ueject)

                    lambdaMeject = atan(uejectvec(2), uejectvec(1))
                    if(lambdaMeject /= lambdaMeject) lambdaMeject = 0d0
                    if(lambdaMeject < 0d0) lambdaMeject = lambdaMeject + twopi
                    
                    ! ejection symmetry axis in the local horizontal CS
                    jetdir = (/sin(s%zeta) * cos(s%eta), sin(s%zeta) * sin(s%eta), cos(s%zeta)/)
                    ! angle between the ejection velocity angle and the ejection axis of symmetry
                    tmp = dot_product(jetdir, uejectvec / ueject)
                    if(abs(tmp) > 1d0) tmp = sign(1d0,tmp)
                    wpsi = acos(tmp)
                    
                    ! the distribution of ejection angle is defined in coordinates (wpsi, wlambdaM) 
                    ! where wpsi is an angle between the jet main axis of symmetry and the direction of ejection
                    ! wlambdaM is a longitude in the plane perpendicular to the jet's main axis
                    ! However, the factor 1/cos(psi) comes from the Jacobian
                    ! of transformation (alphaM, betaM, u, psi, lambdaM) -> (alpha, beta, v, theta, lambda)
                    ! and here psi is the angle between the direction of ejection and the normal to surface
                    fac2 = ejection_direction_distribution(s%ejection_angle_distr, &
                                    wpsi, psieject, lambdaMeject, s%zeta, s%eta)
                                                    
                    if(fac2 == 0d0) then
                        Integrand = 0d0
                    else
                        hcJpsi = hc_jacobian_ueject(ueject, &
                             sqrt(uejectvec(1)**2 + uejectvec(2)**2), &
                             uu, sqrt(uvec(1)**2 + uvec(2)**2))
                        
                        fac2 = fac2 / abs(cos(psi))
                        fac2 = fac2 * abs(hcJpsi)
                        
                        ddphidtheta = derdphidtheta(point%r, velocity, theta, &
                                             muR, s%r, .FALSE., tmp < 1d-3)
                        
                        ddphidv = derdphidv(point%r, velocity, theta, muR, s%r, .FALSE.)
                        
                        call deltat_num_derivatives(muR, velocity, theta, &
                               point%r, s%r, uu, psi, .FALSE. , ddeltatdv, ddeltatdtheta)
                        
                               
                        fac3 = abs(ddeltatdv * ddphidtheta - ddeltatdtheta * ddphidv)
                        
                        Integrand = fac1 * fac2 / fac3 * s%Nparticles
                    endif
                else
                    Integrand = 0d0
                endif
                if((Integrand /= Integrand .or. Integrand < 0.0d0 .or. Integrand > check_inf)) then
                    write(666,*) ' '
                    write(666,*) 'a bad value is obtained for the integrand &
                    in case of delta-ejection:', Integrand
                    write(666,*) 'factor of ejection speed distribution', fac1
                    write(666,*) 'factor of ejection direction distribution', fac2
                    write(666,*) 'Jacobian', fac3
                    write(666,*) 'rM r dphi', s%r, point%r, dphi
                    write(666,*) 'ueject, u, v [m/s]', ueject * AUdays2SI, &
                                uu * AUdays2SI, velocity * AUdays2SI
                    write(666,*) 'uejectvec [AU/day]', uejectvec
                    write(666,*) 'jetdir', jetdir
                    write(666,*) 'eject dir distr', s%ejection_angle_distr
                    write(666,*) 'wpsi psi theta [deg]', wpsi*rad2deg, psi * rad2deg, theta * rad2deg
                    write(666,*) dot_product(jetdir, uejectvec / ueject), s%zeta * rad2deg, s%eta * rad2deg
                    write(666,*) 'lambdaM [deg]', lambdaM * rad2deg
                    write(666,*) 'Jpsi', Jacobian_tilt(psi, lambdaM, s%zeta, s%eta)
                    write(666,*) 'hc_Jpsi', hcJpsi
                    write(666,*) 'ddphidtheta ddphidv', &
                                ddphidtheta, ddphidv
                    write(666,*) 'ddeltatdv, ddeltatdtheta', ddeltatdv, ddeltatdtheta
                    write(666,*) 'umin umax ueject', s%ud%umin, s%ud%umax, ueject
                    write(666,*) 'psi_ast =', acos(dot_product(comet%Vastvec, comet%coords) &
                               /norma3d(comet%coords) / comet%Vast) * rad2deg
                    write(666,*) 'collision', collision
                    N_of_warnings = N_of_warnings + 1
                    if(N_of_warnings > maxNofWarnings) then
                        write(*,*) 'too many warnings have been printed to fort.666'
                        stop
                    endif
                endif

                
            end subroutine Integrand_delta_ejection



            
            ! computes the time interval `dt´ needed for a particle to travel
            ! from `r0´ to `r´ under the influence of gravity defined
            ! by the gravitational parameter `mu´.
            ! v is the speed the particle assumes at `r´
            subroutine deltat(r, r0, v, theta, mu, dt, pericenter)
                implicit none
                real(8), intent(in) :: r, r0, v, theta, mu
                logical, intent(in) :: pericenter
                real(8), intent(out) :: dt
                real(8) a, meanmotion, phi, phi0, ea, ea0, e
                real(8) Ekep, h, e2, u, psi
                real(8) one_minus_e, one_plus_e, pp, thesqrt
                real(8) tra0, tmpsin
                logical ascend
                
                ascend = r > r0 .and. .not. pericenter
                
                a = 1d0 / (2d0 / r - v * v / mu)
                meanmotion = sqrt(abs(mu / a**3))
                Ekep = v * v / 2d0 - mu / r
                u = sqrt(2d0 * (Ekep + mu / r0))
                h = r * v * sin(theta)
                tmpsin = h / r0 / u
                if(tmpsin > 1d0 .or. tmpsin < -1d0) tmpsin = sign(1d0, tmpsin)
                psi = asin(tmpsin)
                e = sqrt(1d0 + 2d0 * Ekep * (h / mu) * (h / mu))
                one_minus_e = 1d0 - e
                one_plus_e = 1d0 + e
                e2 = e * e
                pp = a * (1d0 - e2)
                ea = compute_ea(r, v, theta, mu, ascend)
                ea0 = compute_ea(r0, u, psi, mu, ascend)
                
                tra0 = compute_phi(r0, u, psi, mu, ascend)
                
                if(mu < 0d0) then
                    dt = (e * (sinh(ea) - sinh(ea0)) + ea - ea0) / meanmotion
                else
                    if(a < 0d0) then
                        dt = (e * (sinh(ea) - sinh(ea0)) - ea + ea0) / meanmotion
                        dt = abs(dt)
                    else
                        dt = (ea - ea0 - e * (sin(ea) - sin(ea0))) / meanmotion
                    endif
                endif
                
            end subroutine deltat
            
            
            ! computes derivatives of the time interval \Delta t
            ! for different trajectories: ellipse, hyperbola, hyperbola in case
            ! of a repulsive force
            subroutine deltatderivatives(mu, v, theta, r, r0, u, psi, ddeltatdv, ddeltatdtheta)
                use const
                implicit none
                real(8), intent(in) :: mu, v, theta, r, r0, u, psi
                real(8), intent(out) :: ddeltatdtheta, ddeltatdv
                real(8) a, meanmotion, phi, phi0, ea, ea0, e
                real(8) Ekep, h, dt1, dt0, dt2, dt3
                real(8) one_minus_e, one_plus_e, pp, tg05phi, tg05phi0, thesqrt
                real(8) dadv, dedv, de0dv, dmeanmotiondv, dphidv, dphi0dv, deadv, dea0dv
                real(8) dedtheta, de0dtheta, dphidtheta, dphi0dtheta, deadtheta, dea0dtheta
                real(8) therel, cos05phi2, cos05phi02, e2, sqr1, sqr0
                real(8) factor, factor0, bracket1
                real(8) eps1, eps2, signdr
                
                signdr = sign(1d0, r - r0)
                
                a = 1d0 / (2d0 / r - v * v / mu)
                meanmotion = sqrt(abs(mu / a**3))
                Ekep = v * v / 2d0 - mu / r
                h = r * v * sin(theta)
                e = sqrt(1d0 + 2d0 * Ekep * (h / mu) * (h / mu))
                one_minus_e = 1d0 - e
                one_plus_e = 1d0 + e
                therel = one_minus_e / one_plus_e
                e2 = e * e
                pp = a * (1d0 - e2)
                
                dedv = (2d0 * h * h * v + 4d0 * Ekep * h * r * sin(theta)) &
                       / 2d0 / e / mu**2
                dedtheta = 2d0 * h * Ekep * r * v * cos(theta) / e / mu**2
                
                de0dv = (2d0 * h * h * u + 4d0 * Ekep * h * r0 * sin(psi)) &
                       / 2d0 / e / mu**2
                de0dtheta = 2d0 * h * Ekep * r0 * u * cos(psi) / e / mu**2
                
                dadv = 2d0 / (2d0 / r - v * v / mu)**2 * v / mu
                dmeanmotiondv = -sign(1d0,a) * 1.5d0 * sqrt(abs(mu)) * abs(a)**(-2.5) * dadv
                
                
                if(mu < 0d0) then
                    ea = acosh((r / a - 1d0) / e)
                    ea0 = acosh((r0 / a - 1d0) / e)
                    
                    deadv = -(r * dadv / a**2 / e + (r / a - 1d0) * dedv / e2) &
                            / sqrt((r / a - 1d0)**2 / e2 - 1d0)
                    dea0dv = -(r0 * dadv / a**2 / e + (r0 / a - 1d0) * dedv / e2) &
                            / sqrt((r0 / a - 1d0)**2 / e2 - 1d0)
                    
                    deadtheta = -(r / a - 1d0) * dedtheta / e2 / sqrt((r / a - 1d0)**2 / e2 - 1d0)
                    dea0dtheta = -(r0 / a - 1d0) * dedtheta / e2 / sqrt((r0 / a - 1d0)**2 / e2 - 1d0)
                    
                    ddeltatdv = (sinh(ea) * dedv + deadv + cosh(ea) * e * deadv) / meanmotion &
                    - (ea + e * sinh(ea)) * dmeanmotiondv / meanmotion**2 &
                    -(sinh(ea0) * dedv + dea0dv + cosh(ea0) * e * dea0dv) / meanmotion &
                    + (ea0 + e * sinh(ea0)) * dmeanmotiondv / meanmotion**2
                    
                    ddeltatdtheta = (sinh(ea) * dedtheta + deadtheta + cosh(ea) * e * deadtheta &
                    -sinh(ea0) * dedtheta - dea0dtheta - cosh(ea0) * e * dea0dtheta) / meanmotion
                else
                    phi = acos((pp / r - 1d0) / e)
                    phi0 = acos((pp / r0 - 1d0) / e)
                    
                    if(r0 > r) then
                        phi0 = twopi - phi0
                        phi = twopi - phi
                    endif
                    
                    tg05phi = tan(phi / 2d0)
                    tg05phi0 = tan(phi0 / 2d0)
                    
                    cos05phi2 = cos(phi / 2d0)**2
                    cos05phi02 = cos(phi0 / 2d0)**2
                
                    sqr1 = sqrt(1d0 - (1d0 - pp / r)**2 / e2)
                    sqr0 = sqrt(1d0 - (1d0 - pp / r0)**2 / e2)
                    
                    dphidv = -signdr * ((dadv * (1d0 - e2) - 2d0 * a * e * dedv) / e &
                              - dedv * (a * (1d0 - e2) - r) / e2) / r / sqr1
                    dphi0dv = -signdr * ((dadv * (1d0 - e2) - 2d0 * a * e * de0dv) / e &
                              - de0dv * (a * (1d0 - e2) - r0) / e2) / r0 / sqr0
                    
                    dphidtheta = signdr * (a * e2 - r + a) / r / e2 * dedtheta &
                               / sqr1
                    dphi0dtheta = signdr * (a * e2 - r0 + a) / r0 / e2 * de0dtheta &
                               / sqr0
                    
                    factor = 2d0 / (1d0 + therel * tg05phi**2)
                    factor0 = 2d0 / (1d0 + therel * tg05phi0**2)
                    if(a < 0d0) then
                        thesqrt = sqrt(-therel)
                        bracket1 = one_plus_e**2 * thesqrt
                        
                        ea = 2d0 * atanh(tg05phi * thesqrt)
                        ea0 = 2d0 * atanh(tg05phi0 * thesqrt)
                        if(ea < 0d0) then
                            ea = twopi + ea
                            ea0 = twopi + ea0
                        endif
                        
                        
                        deadv = factor * (tg05phi * dedv / bracket1 &
                              + thesqrt / 2d0 / cos05phi2 * dphidv)
                        dea0dv = factor0 * (tg05phi0 * de0dv / bracket1 &
                              + thesqrt / 2d0 / cos05phi02 * dphi0dv)
                              
                        deadtheta = factor * (tg05phi * dedtheta / bracket1 &
                                  + thesqrt / 2d0 / cos05phi2 * dphidtheta)
                        dea0dtheta = factor0 * (tg05phi0 * de0dtheta / bracket1 &
                                  + thesqrt / 2d0 / cos05phi02 * dphidtheta)
                                  
                        ddeltatdv = -(e * (sinh(ea) - sinh(ea0)) - ea + ea0) &
                                  * dmeanmotiondv / meanmotion**2 &
                                  + ((sinh(ea) - sinh(ea0)) * dedv - ea + ea0 &
                                  + (cosh(ea) * deadv - cosh(ea0) * dea0dv) * e) &
                                  / meanmotion
                        
                        ddeltatdtheta = (dedtheta * (sinh(ea) - sinh(ea0)) &
                                      - deadtheta + dea0dtheta &
                                      + e * (cosh(ea) * deadtheta &
                                      - cosh(ea0) * dea0dtheta)) / meanmotion
                    else
                        thesqrt = sqrt(therel)
                        bracket1 = one_plus_e**2 * thesqrt
                        
                        ea = 2d0 * atan(tg05phi * thesqrt)
                        ea0 = 2d0 * atan(tg05phi0 * thesqrt)
                        
                        if(ea < 0d0) then
                            ea = ea + twopi
                            ea0 = ea0 + twopi
                        endif
                        
                        deadv = -factor * (tg05phi * dedv / bracket1 &
                              - dphidv / 2d0 * thesqrt / cos05phi2)
                        dea0dv = -factor0 * (tg05phi0 * de0dv / bracket1 &
                              - dphi0dv / 2d0 * thesqrt / cos05phi02)
                        
                        deadtheta = -factor * (tg05phi * dedtheta / bracket1  &
                              - dphidtheta / 2d0 * thesqrt / cos05phi2)
                        dea0dtheta = -factor0 * (tg05phi0 * de0dtheta / bracket1  &
                              - dphi0dtheta / 2d0 * thesqrt / cos05phi02)
                        
                        ddeltatdv = (deadv - dea0dv - dedv * (sin(ea) - sin(ea0)) &
                        - e * (cos(ea) * deadv - cos(ea0) * dea0dv)) / meanmotion &
                        - dmeanmotiondv * (ea - ea0 - e * (sin(ea) - sin(ea0))) / meanmotion**2
                        
                        ddeltatdtheta = (deadtheta * (1d0 - e * cos(ea)) - dea0dtheta * (1d0 - e * cos(ea0)) &
                                         - dedtheta * (sin(ea) - sin(ea0))) / meanmotion
                    endif
                endif
                
            
            end subroutine deltatderivatives
            
            
            
            ! derivatives of the time interval \Delta t computed numerically
            subroutine deltat_num_derivatives(mu, v, theta, r, r0, u, psi, &
                              pericenter, ddeltatdv, ddeltatdtheta)
                use const
                implicit none
                real(8), intent(in) :: mu, v, theta, r, r0, u, psi
                real(8), intent(out) :: ddeltatdtheta, ddeltatdv
                logical, intent(in) :: pericenter
                real(8) dt1, dt0, dt2, dt3
                real(8), parameter :: eps1 = 1d-3, eps2 = 1d-4
                
                
                call deltat(r, r0, v+eps2, theta, mu, dt1, pericenter)
                call deltat(r, r0, v-eps2, theta, mu, dt0, pericenter)
                call deltat(r, r0, v+2d0*eps2, theta, mu, dt2, pericenter)
                call deltat(r, r0, v-2d0*eps2, theta, mu, dt3, pericenter)
                ddeltatdv = (-dt2 + 8d0 * dt1 - 8d0 * dt0 + dt3) / 12d0 / eps2
                
                call deltat(r, r0, v, theta+eps1, mu, dt1, pericenter)
                call deltat(r, r0, v, theta-eps1, mu, dt0, pericenter)
                call deltat(r, r0, v, theta+2d0*eps1, mu, dt2, pericenter)
                call deltat(r, r0, v, theta-2d0*eps1, mu, dt3, pericenter)
                ddeltatdtheta = (-dt2 + 8d0 * dt1 - 8d0 * dt0 + dt3) / 12d0 / eps1
                
            
            end subroutine deltat_num_derivatives
            
            
            ! true anomaly for positive mu
            function compute_phi(r, v, theta, mu, ascend) result(phi)
                use const
                implicit none
                real(8) r, v, theta, mu, phi
                real(8) a, e, Ekep, h, tmpcos
                logical ascend
                
                a = 1d0 / (2d0 / r - v * v / mu)
                Ekep = v * v / 2d0 - mu / r
                h = r * v * sin(theta)
                e = sqrt(1d0 + 2d0 * Ekep * (h / mu) * (h / mu))
                tmpcos = (a * (1d0 - e**2) / r - 1d0) / e
                if(tmpcos < -1d0 .or. tmpcos > 1d0) tmpcos = sign(1d0, tmpcos)
                phi = acos(tmpcos)
                if(.not. ascend) then
                    phi = twopi - phi
                endif
            
            end function compute_phi
            
            
            ! eccentric anomaly
            function compute_ea(r, v, theta, mu, ascend) result(ea)
                use const
                implicit none
                real(8) r, v, theta, mu, phi, ea
                real(8) a, e, Ekep, h
                logical ascend
                
                
                if(mu > 0d0) then
                    a = 1d0 / (2d0 / r - v * v / mu)
                    Ekep = v * v / 2d0 - mu / r
                    h = r * v * sin(theta)
                    e = sqrt(1d0 + 2d0 * Ekep * (h / mu) * (h / mu))
                    phi = compute_phi(r, v, theta, mu, ascend)
                    if(a > 0d0) then
                        ea = 2d0 * atan(tan(phi / 2d0) * sqrt((1d0 - e) / (1d0 + e)))
                        if(ea < 0d0) ea = twopi + ea
                    else
                        ea = 2d0 * atanh(tan(phi / 2d0) * sqrt(-(1d0 - e) / (1d0 + e)))
                    endif
                else
                    a = 1d0 / (2d0 / r - v * v / mu)
                    Ekep = v * v / 2d0 - mu / r
                    h = r * v * sin(theta)
                    e = sqrt(1d0 + 2d0 * Ekep * (h / mu) * (h / mu))
                    ea = acosh((r / abs(a) - 1d0) / e)
                endif
                
            end function compute_ea


            
    
            
            ! function control tests if the obtained value of theta is correct
            ! the criterion is: using the formulae from Sremcevic, 2003 paper and the obtained value of theta
            ! one gets the same value of dphi which was used to calculate the theta
            ! also the eccentricity values are compared
            ! input parameters: r0 and rm0 - lengths of two vectors - positions on the orbit,
            ! vv - speed at position r0, phi - angle between r0 and rm0 used in subroutine theta_geometry_...to calculate theta
            ! theta angle between r0 and velosity at the position r0
            ! ee is eccentricity obtained in theta_geometry...
            ! the function uses formulae for energy and angular momentum to obtain values of ee1 and dphi
            ! they must coincide with ee and phi up to 7th digit correspondingly, in this case the function returns TRUE
            ! otherwise it returns FALSE
            ! it also returns the estimated accuracy stored in the variable discr
            subroutine control(muR, theta, ee, phi, vv, r0, rm0, &
                               solved, discr, pericenter)
                use const
                implicit none
                real(8), parameter :: eps = 5d-3
                real(8) hh, hh2, ee1, cosp, cospm, phi1, phi1m, dphi, Ekep
                real(8), intent(in) :: muR, vv, phi, ee, r0, rm0, theta
                logical, intent(in) :: pericenter
                logical, intent(out) :: solved
                real(8), intent(out) :: discr
                
                Ekep = vv * vv / 2d0 - muR / r0
                hh = r0 * vv * sin(theta)
                hh2 = hh * hh
            !  eccentricity (eq 31)
                ee1 = sqrt(1d0 + 2d0 * Ekep * (hh / muR) * (hh / muR))

            !   delta phi (equation 32)
                cosp = (hh2 / r0 / muR - 1d0) / ee1
                cospm = (hh2 / rm0 / muR - 1d0) / ee1
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if(cosp > 1d0) then
                    cosp = 1d0
                    write(666,*) 'cos(phi) > 1 obtained, corrections applied'
                    N_of_warnings = N_of_warnings + 1
                    if(N_of_warnings > maxNofWarnings) then
                        write(*,*) 'too many warnings have been printed to fort.666'
                        stop
                    endif
                endif
                if(cosp < -1d0) then
                    cosp = -1d0
                    write(666,*) 'cos(phi) < -1 obtained, corrections applied'
                    N_of_warnings = N_of_warnings + 1
                    if(N_of_warnings > maxNofWarnings) then
                        write(*,*) 'too many warnings have been printed to fort.666'
                        stop
                    endif
                endif
                
                if(cospm > 1d0) then
                    cospm = 1d0
                    write(666,*) 'control << cos(phiM) > 1 obtained, corrections applied'
                    N_of_warnings = N_of_warnings + 1
                    if(N_of_warnings > maxNofWarnings) then
                        write(*,*) 'too many warnings have been printed to fort.666'
                        stop
                    endif
                endif
                if(cospm < -1d0) then
                    cospm = -1d0
                    write(666,*) 'cos(phiM) < -1 obtained, corrections applied'
                    N_of_warnings = N_of_warnings + 1
                    if(N_of_warnings > maxNofWarnings) then
                        write(*,*) 'too many warnings have been printed to fort.666'
                        stop
                    endif
                endif
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                phi1 = acos(cosp)
                phi1m = acos(cospm)
                if(pericenter) then
                    dphi = twopi - (phi1 + phi1m)
                else
                    dphi = abs(phi1 - phi1m)
                endif
                if(dphi > pi) dphi = twopi - dphi
                discr = abs(dphi - phi) / abs(phi)
                solved = discr < eps
                if(pericenter) solved = solved .and. theta > 0.1
                
            end subroutine control
            


        ! subroutine theta_geometry_hyperbola recieves vectors' absolute
        ! values and the angle between these vectors
        ! it's assumed that the point (0, 0, 0) is a focus of an hyperbola
        ! value of a semi major axis of the hyperbola, the points lay on is
        ! also an input parameter of the subroutine
        ! the subroutine returns theta the angle between radius-vector r
        ! and a tangent to the hyperbola in the point r
        ! the choice between two possible values of this angle is made
        ! in the way that movement from point rm to the point r along
        ! the hyperbola is possible
        ! (means: r and rm lay on the same quadrant in the CS in which
        ! the equation of the hyperbola is in canonical form)
        ! In general case there are two possible hyperbolae
        ! => two possible values of theta are to be investigated
        ! Returns array of 2 vallues of theta.
        ! If theta = large negative number,
        ! it means "no physically plausible solution can be found"
        subroutine theta_geometry_hyperbola(muR, r0, rm0, vv, phi, a0, &
                                        ee, theta_res, deltat_res, pericenter)
            use help
            use const
            implicit none
            real(8), intent(in) :: r0, rm0, a0, phi,  vv, muR
            real(8), intent(out) :: ee, theta_res, deltat_res
            logical, intent(in) :: pericenter
            real(8) r, rmoon, a
            real(8) theta, deltat
            real(8) c, x(2), y(2), xm, ym, b
            real(8) shift(2), r2d(2), rm2d(2), angle, tangent(2), tmp(2)
            real(8) cosf1, eanm, ean, psi, f1, f2, discr
            real(8) one_plus_e, one_minus_e, one_minus_e2, aux, tmp1, tmp2
            integer i
            logical solved
            
            solved = .FALSE.
            theta = -555d0
                                                            
            r = r0 / rm0 ; rmoon = rm0 / rm0; a = a0 / rm0                            
            ! define x-axis in the same direction as r-vector
            r2d(1) = r ; r2d(2) = 0d0                    
            
            ! we don't have enough information to define the sign of r
            ! vector in the right-handed coordinate system
            ! but the value of theta that we are looking for are the same
            ! in both cases
            
            rm2d(1) = rmoon * cos(phi) ; rm2d(2) = rmoon * sin(phi)
            
            ! (x(1),y(1)) and (x(2),y(2)) are coordinates
            ! of 2 possible position of the hyperbola's second focus
            call circle_intersection(rm2d(1), rm2d(2), 2d0 * a + rmoon, &
                                        r2d(1), 2d0 * a + r, x, y)
            do i = 1, 2
!~             write(*,*) x(i), y(i), 0.5d0 * sqrt(x(i)**2 + y(i)**2) / a
                ! shift is coordinates of the hyperbola's center
                ! in the CS centered at the focus
                shift(1) = x(i) / 2d0 ; shift(2) = y(i) / 2d0        
                
                ! coords of vector r in CS centered at center of the hyperbola
                r2d = r2d - shift    
                ! coords of vector rm in CS centered at center of the hyperbola                                
                rm2d = rm2d - shift                                
                
                ! angle between major axis of the ellipse and the current x-axis
                angle = atan(y(i) / x(i))
                ! vector r in the CS with its center at the center of the ellipse
                ! and the x-axis along major axis of the ellipse                                    
                r2d = rot2d(r2d, -angle)    
                ! vector rm in the CS with its center at the center of the ellipse
                ! and the x-axis along major axis of the ellipse                                
                rm2d = rot2d(rm2d, -angle)    
                ! vector shift in the CS with its center at the center of the ellipse
                ! and the x-axis along major axis of the ellipse                            
                shift = rot2d(shift, -angle)
                ! the hyperbola solves our problem only if r and rm lay
                ! on the same branch and the trajectory doesn't intersect
                ! the moon's surface (the particle doesn't pass the pericenter
                if(r2d(1) / rm2d(1) > 0d0) then                    
                    c = 0.5d0 * sqrt(x(i)**2 + y(i)**2)                                    
                    ee = c / a
                    one_plus_e = 1d0 + ee
                    one_minus_e = 1d0 - ee
                    one_minus_e2 = one_plus_e * one_minus_e
                    aux = -a * one_minus_e2
                    cosf1 = (aux - 1d0) / ee
                    f1 = acos(cosf1)
                    if(r0 < rm0 .or. pericenter) f1 = -f1
                    f2 = f1 + phi
                    
                    theta = halfpi - atan((ee * sin(f2)) / (1d0 + ee * cos(f2)))
                                            
                    call control(muR, theta, ee, phi, vv, r0, rm0, solved, discr, pericenter)
                    ean = 2d0 * atanh(tan(f2/2d0) * sqrt(-one_minus_e / one_plus_e))
                    eanm = 2d0 * atanh(tan(f1/2d0) * sqrt(-one_minus_e / one_plus_e))
                    
                    deltat = (ee * (sinh(ean) - sinh(eanm)) - ean + eanm) / sqrt(muR / a0**3)
                    if(solved) then
                        theta_res = theta 
                        deltat_res = deltat 
                    endif
                        
                else
                    theta = -444d0
                    ee = -444d0
                endif
                ! we have changed the vectors r and rmoon we started from
                ! so we need to go back to the beginning
                ! to find the second value of theta
                r = r0 / rm0 ; rmoon = rm0 / rm0; a = a0 / rm0
                rm2d(1) = rmoon * cos(phi) ; rm2d(2) = rmoon * sin(phi)
                r2d(1) = r ; r2d(2) = 0d0
                
            enddo

        end subroutine theta_geometry_hyperbola


        ! subroutine theta_geometry_hyperbola recieves vectors' absolute
        ! values and the angle between these vectors
        ! it's assumed that the point (0, 0, 0) is a focus of an hyperbola
        ! value of a semi major axis of the hyperbola, the points lay on is
        ! also an input parameter of the subroutine
        ! the subroutine returns theta the angle between radius-vector r
        ! and a tangent to the hyperbola in the point r
        ! the choice between two possible values of this angle is made
        ! in the way that movement from point rm to the point r along
        ! the hyperbola is possible
        ! (means: r and rm lay on the same quadrant in the CS in which
        ! the equation of the hyperbola is in canonical form)
        ! In general case there are two possible hyperbolae
        ! => two possible values of theta are to be investigated
        ! Returns array of 2 vallues of theta.
        ! If theta = large negative number,
        ! it means "no physically plausible solution can be found"
        subroutine theta_geometry_other_branch(muR, r0, rm0, vv, phi, a0, &
                                        ee, theta_res, deltat_res, pericenter)
            use help
            use const
            implicit none
            real(8), intent(in) :: r0, rm0, a0, phi,  vv, muR
            real(8), intent(out) :: ee, theta_res, deltat_res
            logical, intent(in) :: pericenter
            real(8) r, rmoon, a
            real(8) theta, deltat
            real(8) c, x(2), y(2), xm, ym, b, vec1(2), vec2(2)
            real(8) shift(2), r2d(2), rm2d(2), angle, tangent(2), tmp(2)
            real(8) eanm, ean, psi, discr
            real(8) one_plus_e, one_minus_e, one_minus_e2, aux, tmp1, tmp2
            integer i
            logical solved
            real(16) x1, y1, x2, r1, r2, s1(2), s2(2)
            
            solved = .FALSE.
            theta_res = -505d0
            deltat_res = -222d0
                                                            
            r = r0 / rm0 ; rmoon = rm0 / rm0; a = a0 / rm0                            
            ! define x-axis in the same direction as r-vector
            r2d(1) = r ; r2d(2) = 0d0                    
            
            ! we don't have enough information to define the sign of r
            ! vector in the right-handed coordinate system
            ! but the value of theta that we are looking for are the same
            ! in both cases
            
            rm2d(1) = rmoon * cos(phi) ; rm2d(2) = rmoon * sin(phi)
            
            ! (x(1),y(1)) and (x(2),y(2)) are coordinates
            ! of 2 possible position of the hyperbola's second focus
            call circle_intersection(rm2d(1), rm2d(2), rmoon - 2d0 * a, &
                                        r2d(1), r - 2d0 * a, x, y)
                
            do i = 1, 2
                ! shift is coordinates of the hyperbola's center
                ! in the CS centered at the focus
                shift(1) = x(i) / 2d0 ; shift(2) = y(i) / 2d0        
                
                ! coords of vector r in CS centered at center of the hyperbola
                r2d = r2d - shift    
                ! coords of vector rm in CS centered at center of the hyperbola                                
                rm2d = rm2d - shift                                
                
                ! angle between major axis of the hyperbola and the current x-axis
                angle = atan(y(i) / x(i))
                ! vector r in the CS with its center at the center of the hyperbola
                ! and the x-axis along major axis of the ellipse                                    
                r2d = rot2d(r2d, -angle)    
                ! vector rm in the CS with its center at the center of the hyperbola
                ! and the x-axis along major axis of the ellipse                                
                rm2d = rot2d(rm2d, -angle)    
                ! vector shift in the CS with its center at the center of the hyperbola
                ! and the x-axis along major axis of the hyperbola                            
                shift = rot2d(shift, -angle)
                ! the hyperbola solves our problem only if r and rm lay
                ! on the same branch and the trajectory doesn't intersect
                ! the moon's surface (the particle doesn't pass the pericenter
                if(r2d(1) / rm2d(1) > 0d0) then
                    c = 0.5d0 * sqrt(x(i)**2 + y(i)**2)
                    b = sqrt(c**2 - a**2)
                    ee = c / a
                    one_plus_e = 1d0 + ee
                    one_minus_e = 1d0 - ee
                    one_minus_e2 = one_plus_e * one_minus_e
                    aux = b**2 / a
                    
                    vec1 = (/r, 0d0/)
                    vec2 = vec1 - (/x(i), y(i)/)
                    theta = 0.5d0 * acos(dot_product(vec1, vec2) / r / norma2d(vec2))
                    if(r < rmoon .and. .not. pericenter) theta = pi - theta
                    ean = acosh((r / abs(a) - 1d0) / ee)
                    eanm = acosh((rmoon / abs(a) - 1d0) / ee)
                    if(pericenter) eanm = - eanm
                    deltat = abs(ee * (sinh(ean) - sinh(eanm)) + ean - eanm) / sqrt(-muR / a0**3)
                    call control(muR, theta, ee, phi, vv, r0, rm0, solved, discr, pericenter)
                    
                    if(solved) then
                        theta_res = theta 
                        deltat_res = deltat 
                    endif
                else
                    theta = -440d0
                    ee = -440d0
                endif
                ! we have changed the vectors r and rmoon we started from
                ! so we need to go back to the beginning
                ! to find the second value of theta
                r = r0 / rm0 ; rmoon = rm0 / rm0; a = a0 / rm0
                rm2d(1) = rmoon * cos(phi) ; rm2d(2) = rmoon * sin(phi)
                r2d(1) = r ; r2d(2) = 0d0
                
            enddo
            
            
        end subroutine theta_geometry_other_branch






        ! subroutine theta_geometry_ellipse recieves absolute values
        ! of two vectors and the angle between these vectors
        ! it's assumed that the point (0, 0, 0) is a focus of an ellipse
        ! value of a semi major axis of an ellipse, the points lay on
        ! is also an input parameter of the subroutine
        ! the subroutine returns theta that is the angle between
        ! radius-vector r and a tangent to the ellipse in the point r
        ! the choice between two possible values of this angle is made
        ! in the way that movement happens from point rm to the point r 
        ! In general case there are two possible ellipses
        ! => two possible values of theta are to be found
        ! Returns array of 2 vallues of theta.
        ! If theta = large negative number
        ! it means "no physically plausible solution can be found"
        subroutine theta_geometry_ellipse(muR, r0, rm0, vv, phi, a0, &
                                ee_res, theta_res, deltat, pericenter)
            use help
            use const
            implicit none
            real(8), intent(in) :: r0, rm0, a0, phi,  vv, muR
            real(8), intent(out) :: ee_res, theta_res, deltat
            logical, intent(in) :: pericenter
            real(8) r, rmoon, a, discr, theta(2)
            real(8) cc(2), x(2), y(2), ee(2)
            real(8) r2d(2), rm2d(2), tmp(2), tmp1, tmp2
            real(8) E1(1), f1, f2, sinf1, cosf1, ean, rtest
            real(8) one_plus_e, one_minus_e, one_minus_e2, eanm
            real(8) aux
            integer i
            logical solved
            
            solved = .FALSE.
            theta = -888d0
            deltat = -333d0
            
            r = r0 / rm0  ; rmoon = rm0 / rm0; a = a0 / rm0
            ! define x-axis in the same direction as rmoonvector
            r2d(1) = r; r2d(2) = 0d0                    
            ! we don't have enough information to define the sign of r
            ! vector in the right-handed coordinate system
            ! but the value of theta that we are looking for are the same
            ! in both cases
            rm2d(1) = rmoon * cos(phi) ; rm2d(2) = rmoon * sin(phi)
            
            ! (x(1),y(1)) and (x(2),y(2)) are coordinates of 2 possible
            ! position of the ellsipse's second focus
            call circle_intersection(rm2d(1), rm2d(2), 2d0 * a - rmoon, &
                                    r2d(1), 2d0 * a - r, x, y) 
            ! distance between the foci of the ellipse
            cc = sqrt(x**2 + y**2)
            
            do i = 1, 2
                if(cc(i) == cc(i)) then
                    ee(i) = cc(i) / 2d0 / a

                    one_plus_e = 1d0 + ee(i)
                    one_minus_e = 1d0 - ee(i)
                    one_minus_e2 = one_plus_e * one_minus_e
                    cosf1 = (-x(i) * rm2d(1) - y(i) * rm2d(2)) / rmoon / cc(i)
                    f1 = acos(cosf1)
                    if((r0 < rm0 .or. pericenter)) f1 = -f1
                    f2 = f1 + phi
                    
                    theta(i) = halfpi - atan((ee(i) * sin(f2)) &
                                            / (1d0 + ee(i) * cos(f2)))
                        
                    call control(muR, theta(i), ee(i), phi, vv, r0, rm0, &
                                solved, discr, pericenter)
                    if(solved) then
                        ean = 2d0 * atan(tan(f2/2d0) &
                                * sqrt(one_minus_e / one_plus_e))
                        if(ean < 0d0) ean = twopi + ean
                        eanm = 2d0 * atan(tan(f1/2d0) &
                                * sqrt(one_minus_e / one_plus_e))
                        if(eanm < 0d0 .and. .not. pericenter) eanm = twopi + eanm
                        tmp1 = a0 * sqrt(a0 / muR) * (eanm - ee(i) * sin(eanm))
                        tmp2 = a0 * sqrt(a0 / muR) * (ean - ee(i) * sin(ean))
                        deltat = tmp2 - tmp1
                        theta_res = theta(i)
                        ee_res = ee(i)
                    endif
                else
                    theta(i) = -888d0
                    ee(i) = -888d0
                endif
                                    
            enddo
        
        
        end subroutine theta_geometry_ellipse




        ! Integrand_v_integration performs integration over theta and lambda,
        ! returns the expression standing under the integral over v
        ! this subroutine is called `order_v´ times for each pair
        ! 'a source + a point'
        subroutine Integrand_v_integration(Integrand, velocity,  dt, &
                                            point, dphi, dbeta, s, &
                                            tnow, muR, comet, &
                                            Rast_AU, pericenter)
            use const
            use define_types
            use help
            use distributions_fun
            implicit none
            integer i
            real(8), intent(in) :: velocity, dphi, dbeta, tnow
            real(8), intent(out) :: dt
            real(8), intent(in) :: muR, Rast_AU
            logical, intent(in) :: pericenter
            type(ephemeris), intent(in) :: comet
            type(source_properties), intent(in) :: s
            type(position_in_space), intent(in) :: point
            real(8) lambdaM, sinlambdaM, coslambdaM, lambda
            real(8) hcJpsi, Vastvechorz(3), tmp, jetdir(3)
            real(8) ueject, uejectvec(3), psieject, lambdaMeject, uvec(3)
            real(8) uu, Ekep, psi, wpsi, hh
            real(8) ee, semi_major_axis, moment_of_ejection
            real(8), intent(out) :: Integrand
            real(8) theta, ddphidtheta
            real(8) fac1, fac2, rate
            logical collision

            semi_major_axis = (2d0 / point%r - velocity**2 / muR)**(-1)
            theta = -999d0
            dt = -111d0
            collision = .False.
            ! find solutions for theta
            if(semi_major_axis > 0d0 .and. muR > 0d0) then
                call theta_geometry_ellipse(muR, point%r, s%r, velocity, &
                                  dphi, semi_major_axis, ee, theta, dt, &
                                  pericenter)
            endif
            if(semi_major_axis < 0d0 .and. muR > 0d0) then

                call theta_geometry_hyperbola(muR, point%r, s%r, velocity, &
                               dphi, abs(semi_major_axis), ee, theta, dt, &
                                  pericenter)

            endif
            if(muR < 0d0) then
                call theta_geometry_other_branch(muR, point%r, s%r, velocity, &
                               dphi, abs(semi_major_axis), ee, theta, dt, &
                                  pericenter)
            endif
            Integrand = 0d0
            wpsi = halfpi
            moment_of_ejection = tnow - dt
            if(s%Tj-s%dtau/2d0 <= moment_of_ejection &
                           .and. moment_of_ejection <= s%Tj + s%dtau/2d0) then
                Ekep = velocity**2 / 2d0 - muR / point%r
                uu = sqrt(2d0 * (Ekep + muR / s%r))
                ! angular momentum (eq 27)
                hh = point%r * velocity * sin(theta)
                !  psi
                tmp = hh / s%r / uu
                if(tmp <= -0.99999d0) tmp = -0.99999d0
                if(tmp >= 0.999999d0) tmp = 0.999999d0
                psi = asin(tmp)
                if((theta > halfpi .and. .not. pericenter) .or. pericenter) psi = pi - psi
                ! azimuth of the particle velocity vector at the moment of ejection
                lambdaM = azimuth(point%rvector - s%rrM, s%rrM)
                
                ! particle velocity vector at the moment of ejection
                uvec = uu * (/sin(psi) * cos(lambdaM), sin(psi) * sin(lambdaM), cos(psi)/)
                
                ! ejection velocity vector in the local horizontal CS
                if(comet%Vast /= 0d0) then
                    call find_astapex(s%rrM, s%r, comet%Vastvec, comet%Vast, Vastvechorz)
                    uejectvec = uvec - Vastvechorz
                else
                    Vastvechorz = 0d0
                    uejectvec = uvec
                endif
                ueject = norma3d(uejectvec)
                
                if(Rast_AU > 0d0) then
                    call collision_check(muR, dt, comet%coords, s, point, &
                           uejectvec, collision, Rast_AU)
                else
                    collision = .False.
                endif
                
                fac1 = 0d0 ; fac1 = 0d0
                if((.not. collision) .and. &
                    ueject > s%ud%umin .and. ueject < s%ud%umax) then
                    fac1 = ejection_speed_distribution(s%ud, ueject)
                    fac1 = velocity * fac1 / uu / uu

                    psieject = acos(uejectvec(3) / ueject)
                    if(psieject /= psieject) then
                        write(*,*) uejectvec*AUdays2SI, ueject*AUdays2SI
                    endif
                    lambdaMeject = atan(uejectvec(2), uejectvec(1))
                    if(lambdaMeject /= lambdaMeject) lambdaMeject = 0d0
                    if(lambdaMeject < 0d0) lambdaMeject = lambdaMeject + twopi
                    
                    ! ejection symmetry axis in the local horizontal CS
                    jetdir = (/sin(s%zeta) * cos(s%eta), sin(s%zeta) * sin(s%eta), cos(s%zeta)/)
                    ! angle between the ejection velocity angle and the ejection axis of symmetry
                    tmp = dot_product(jetdir, uejectvec / ueject)
                    if(abs(tmp) > 1d0) tmp = sign(1d0,tmp)
                    wpsi = acos(tmp)
                    
                    ! the distribution of ejection angle is defined in coordinates (wpsi, wlambdaM) 
                    ! where wpsi is an angle between the jet main axis of symmetry and the direction of ejection
                    ! wlambdaM is a longitude in the plane perpendicular to the jet's main axis
                    ! However, the factor 1/cos(psi) comes from the Jacobian
                    ! of transformation (alphaM, betaM, u, psi, lambdaM) -> (alpha, beta, v, theta, lambda)
                    ! and here psi is the angle between the direction of ejection and the normal to surface
                    fac2 = ejection_direction_distribution(s%ejection_angle_distr, &
                                        wpsi, psieject, lambdaMeject, s%zeta, s%eta)
                    if(fac2 == 0d0 .or. fac1 == 0d0) then
                        Integrand = 0d0
                    else
                        
                        hcJpsi = hc_jacobian_ueject(ueject, &
                             sqrt(uejectvec(1)**2 + uejectvec(2)**2), &
                             uu, sqrt(uvec(1)**2 + uvec(2)**2))
                        
                        tmp = abs(cos(psi))
                        fac2 = fac2 / tmp
                        fac2 = fac2 * abs(hcJpsi)
                        
                        ddphidtheta = derdphidtheta(point%r, velocity, theta, &
                        muR, s%r, pericenter, tmp < 1d-3)
                        
                        Integrand = fac1 * fac2 / abs(ddphidtheta) * s%Nparticles / s%dtau
                    endif
                else
                    Integrand = 0d0
                endif
            else
                Integrand = 0d0
            endif
            if(Integrand /= Integrand .or. Integrand < 0.0d0 &
                    .or. Integrand > check_inf) then
                write(666,*) ' '
                write(666,*) 'a bad value is obtained for the integrand &
                    in case of delta-ejection:', Integrand
                write(666,*) 'factor of ejection speed distribution', fac1
                write(666,*) 'factor of ejection direction distribution', fac2
                write(666,*) 'Jacobian', 1d0 / abs(ddphidtheta)
                write(666,*) 'rM r dphi', s%r, point%r, dphi
                write(666,*) 'semi major axis', semi_major_axis, 'muR', muR
                write(666,*) 'ueject, u, v [m/s]', uejectvec * AUdays2SI, &
                uu * AUdays2SI, velocity * AUdays2SI
                write(666,*) 'dt', s%Tj-s%dtau, moment_of_ejection, s%Tj+s%dtau
                write(666,*) 'wpsi psi theta [deg]', wpsi*rad2deg, psi * rad2deg, &
                           theta * rad2deg
                write(666,*) 'lambdaM [deg]', lambdaM * rad2deg
                write(666,*) 'Jpsi', Jacobian_tilt(psi, lambdaM, s%zeta, s%eta)
                write(666,*) 'hc_Jpsi', hcJpsi
                write(666,*) 'ddphidtheta', &
                            ddphidtheta
                write(666,*) 'psi_ast =', acos(dot_product(comet%Vastvec, &
                comet%coords)/norma3d(comet%coords) / comet%Vast) * rad2deg
                write(666,*) 'collision', collision
                write(666,*) '  '
                N_of_warnings = N_of_warnings + 1
                if(N_of_warnings > maxNofWarnings) then
                    write(*,*) 'too many warnings have been printed to fort.666'
                    stop
                endif
            endif
            
        end subroutine Integrand_v_integration


end module twobody_fun
