! This file is a part of DUDI-heliocentric, the Fortran-95 implementation 
! of the two-body model for the dynamics of dust ejected from an atmosphereless
! body moving around the Sun
! Version 1.0.1
! This is free software. You can use and redistribute it 
! under the terms of the GNU General Public License (http://www.gnu.org/licenses/)
! If you do, please cite the following paper
! Anastasiia Ershova and Jürgen Schmidt, 
! Two-body model for the spatial distribution of dust ejected from
! an atmosphereless body, 2021, A&A, 650, A186 

! File: DUDIhc.f90
! Description: Contains ssubroutines that implement three solution methods
! on a high level and auxiliary subroutines that support these methods.


! Author: Anastasiia Ershova
! E-mail: vveyzaa@gmail.com

module DUDIhc
    implicit none
        contains
        
        ! The subroutine finds the speed and direction of a dust grain
        ! ejection. Using it, the subroutine calculates the dust number
        ! density at the given point.
        subroutine hc_DUDI_simple_expansion(density, source, dt, &
                                             cloudcentr, point)
            use const
            use define_types
            use distributions_fun
            use help
            implicit none
            real, intent(out) :: density
            real(8), intent(in) :: dt, cloudcentr(3)
            type(position_in_space), intent(in) :: point 
            type(source_properties), intent(in) :: source 
            real(8) ueject, psi, lambdaM, dist, fac1, fac2
            real(8) uejectvec(3), tmp
            ! guard cloud center norm
            real(8) :: cnorm, sinpsi, dens8
            
            dist = norma3d(point%rvector - cloudcentr)
            if (dist <= tiny(1.0d0)) then
                density = 0.0
                return
            end if

            if (dt <= 0.0d0) then
                density = 0.0
                return
            end if

            ueject = dist / dt
            uejectvec = (point%rvector - cloudcentr) / dist        ! unit vector

            cnorm = norma3d(cloudcentr)
            if (cnorm <= tiny(1.0d0)) then
                density = 0.0
                return
            end if

            tmp = dot_product(uejectvec, cloudcentr) / cnorm
            if (tmp >= 1d0 .or. tmp <= -1d0) tmp = sign(9.9999d-1, tmp)
            psi = acos(tmp)
            lambdaM = azimuth(point%rvector - cloudcentr, cloudcentr)

            fac1 = 0d0 ; fac2 = 0d0
            if (source%ud%umin <= ueject .and. ueject <= source%ud%umax) then
                fac1 = ejection_speed_distribution(source%ud, ueject)
                fac2 = ejection_direction_distribution(source%ejection_angle_distr, &
                       psi, psi, lambdaM, 0d0, 0d0)

                ! compute in double, guard sin(psi)
                sinpsi = sin(psi)
                if (abs(sinpsi) <= 1.0d-15) then
                    density = 0.0
                    return
                end if

                dens8 = source%Nparticles * fac1 * fac2 / dt
                dens8 = dens8 / (AU**3 * dist**2 * sinpsi)
                density = real(dens8, kind=kind(0.0))
            else
                density = 0.0
            end if

        
        end subroutine hc_DUDI_simple_expansion
        
        

        ! The function finds the only possible values of speed u
        ! and zenith angle theta of velocity vector of the particle ejected 
        ! from the source and detected at the point
        ! these values are supplied to the function calculating Integrand
        ! Thus, the integration itself is performed analytically
        subroutine hc_DUDI_delta_ejection(density, point, source, &
                    muR, dt, comet, Rast_AU)
            use const
            use define_types
            use nan_utils
            use help
            use twobody_fun
            implicit none
            real(8), parameter :: eps = 1d-8
            integer i, ii
            real, intent(out) :: density
            real(8), intent(in) :: dt, muR, Rast_AU
            type(ephemeris), intent(in) :: comet
            real(8) u, psi,  v, tmp2(2)
            real(8) dr(3),  midval
            type(position_in_space), intent(in) :: point
            type(source_properties), intent(in) :: source
            real(8) dphi, dbeta
            real(8) r2d(2), rm2d(2), ax, ux, uy, theta
            real(8)  Vastx, Vasty, discr
            real(8) rm2, rm3, rm4, rm5, rm6, rm7, rm8, rm9, rm10, rm11
            real(8) ux2, ux3, uy2, uy3, ux4, uy4, ux5, uy5
            real(8) dt2, dt3, dt4, dt5, dt6, dt7, dt8
            real(8) Ekep, hh, ee
            real(8) coefsx(0:4), coefsy(0:3), tmpux, tmpuy
            real(8) xcoefs(0:6), ycoefs(0:5), muR2, muR3, muR4
            
            call ApuTrajectoryLight(point, dphi, dbeta, source)
            
            ! 2d coordinates of r and rM in the orbital plane
            r2d(1) = point%r * cos(dphi) ; r2d(2) = point%r * sin(dphi)
            rm2d(1) = source%r ; rm2d(2) = 0d0
            
            ! auxiliary variables
            rm2 = source%r**2
            rm3 = source%r * rm2
            rm4 = rm2 * rm2
            rm5 = rm3 * rm2
            rm6 = rm4 * rm2
            rm7 = rm4 * rm3
            rm8 = rm4 * rm4
            rm9 = rm6 * rm3
            rm10 = rm5 * rm5
            rm11 = rm9 * rm2
            
            dt2 = dt * dt
            dt3 = dt2 * dt
            dt4 = dt3 * dt
            dt5 = dt4 * dt
            dt6 = dt5 * dt
            dt7 = dt6 * dt
            dt8 = dt7 * dt

            muR2 = muR * muR
            muR3 = muR2 * muR
            muR4 = muR3 * muR

!~             ! x-coordinate of the comet velosity
            Vastx = dot_product(comet%Vastvec, source%rrM) / source%r
            Vasty = sqrt(comet%Vast**2 - Vastx**2)
            if(is_zero_r8(Vastx)) Vastx = 2d5 / AUdays2SI
            ! initial approximate solution
            tmpux = Vastx
            tmpuy = Vasty 
            ux = 0d0
            uy = 0d0
            ii = 0
            discr = 1d+13
            do while(discr > eps)
                ux = tmpux
                uy = tmpuy

                ux2 = ux * ux ; uy2 = uy * uy
                ux3 = ux2 * ux ; uy3 = uy2 * uy
                ux4 = ux2 * ux2 ; uy4 = uy2 * uy2
                ux5 = ux3 * ux2 ; uy5 = uy2 * uy3
                
                xcoefs(0) = (-73d0 * dt8 * muR4) / (5040d0 * rm11) &
                          - (11d0 * dt6 * muR3) / (360d0 * rm8) &
                          - (dt4 * muR2) / (12d0 * rm5) - (dt2 * muR) / (2d0 * rm2) &
                          + rm2d(1) + (73d0 * dt8 * muR3 * uy2) / (1120d0 * rm10) &
                          + (11d0 * dt6 * muR2 * uy2) / (120d0 * rm7) &
                          + (dt4 * muR * uy2) / (8d0 * rm4) & 
                          -(201d0 * dt8 * muR2 * uy4) / (2240d0 * rm9) &
                          - (dt6 * muR * uy4) / (16d0 * rm6) &
                           + (5d0 * dt8 * muR * uy4 * uy2) / (128d0 * rm8) - r2d(1)
                           
                xcoefs(1) = dt + (73d0 * dt7 * muR3) / (630d0 * rm9) &
                          + (11d0 * dt5 * muR2) / (60d0 * rm6) &
                          + (dt3 * muR) / (3d0 * rm3) - (13d0 * dt7 * muR2 * uy2) / (35d0 * rm8) &
                          - (3d0 * dt5 * muR * uy2) / (10d0 * rm5) &
                          + (15d0 * dt7 * muR * uy4) / (56d0 * rm7)
                xcoefs(2) = (-61d0 * dt8 * muR3) / (224d0 * rm10) &
                          - (17 * dt6 * muR2) / (60d0 * rm7) &
                          - (dt4 * muR) / (4d0 * rm4) &
                          + (519d0 * dt8 * muR2 * uy2) / (560d0 * rm9) &
                          + (dt6 * muR * uy2) / (2d0 * rm6) &
                          - (45d0 * dt8 * muR * uy4) / (64d0 * rm8)
                xcoefs(3) = (53d0 * dt7 * muR2) / (140d0 * rm8) &
                          + (dt5 * muR) / (5d0 * rm5) - (5d0 * dt7 * muR * uy2) &
                          / (7d0 * rm7)
                xcoefs(4) = (-131d0 * dt8 * muR2) / (280d0 * rm9) &
                          - (dt6 * muR) / (6d0 * rm6) + (15d0 * dt8 * muR * uy2) &
                          / (16d0 * rm8)
                xcoefs(5) = (dt7 * muR) / (7d0 * rm7)
                xcoefs(6) = -(dt8 * muR) / (8d0 * rm8)
                
                ycoefs(0) = -r2d(2)
                ycoefs(1) = dt - (43d0 * dt7 * muR3) / (1260d0 * rm9) &
                          - (dt5 * muR2) / (15d0 * rm6) - (dt3 * muR) &
                          / (6d0 * rm3) + (13d0 * dt8 * muR3 * ux) &
                          / (80d0 * rm10) + (5d0 * dt6 * muR2 * ux) &
                          / (24d0 * rm7) + (dt4 * muR * ux) / (4d0 * rm4) &
                          - (59d0 * dt7 * muR2 * ux2) / (140d0 * rm8) &
                          - (3d0 * dt5 * muR * ux2) / (10d0 * rm5) &
                          + (7d0 * dt8 * muR2 * ux3) / (10d0 * rm9) &
                          + (dt6 * muR * ux3) / (3d0 * rm6) - (5d0 * dt7 * muR * ux4) &
                          / (14d0 * rm7) + (3d0 * dt8 * muR * ux5) / (8d0 * rm8)
                ycoefs(2) = 0d0
                ycoefs(3) = (11d0 * dt7 * muR2) / (140d0 * rm8) &
                          + (3d0 * dt5 * muR) / (40d0 * rm5) &
                          - (63d0 * dt8 * muR2 * ux) / (160d0 * rm9) &
                          - (dt6 * muR * ux) / (4d0 * rm6) &
                          + (15d0 * dt7 * muR * ux2) / (28d0 * rm7) &
                          - (15d0 * dt8 * muR * ux3) / (16d0 * rm8)
                ycoefs(4) = 0d0
                ycoefs(5) = (-5d0 * dt7 * muR) / (112d0 * rm7) &
                          + (15d0 * dt8 * muR * ux) / (64d0 * rm8)
                
                xcoefs = xcoefs / xcoefs(1)
                ycoefs = ycoefs / ycoefs(1)
                tmpux = -xcoefs(0) - xcoefs(2) * ux2 - xcoefs(3) * ux3 &
                       - xcoefs(4) * ux2**2 - xcoefs(5) * ux5 - xcoefs(6) * ux * ux5
                tmpuy = -ycoefs(0) - ycoefs(3) * uy3 - ycoefs(5) * uy5
                ii = ii + 1
                discr = sqrt((ux - tmpux)**2 + (uy - tmpuy)**2)
            enddo
            
            u = sqrt(ux**2 + uy**2)
            psi = dble(atan(uy / ux))
            if(ux < 0.0) psi = pi + psi
                        
            Ekep = u * u / 2d0 - muR / source%r
            hh = source%r * u * sin(psi) 
            ee = sqrt(1d0 + 2d0 * Ekep * (hh / muR) * (hh / muR))
            v = sqrt2d0 * sqrt(Ekep + muR / point%r)
            theta = asin(hh / point%r / v)
            if(psi > halfpi) theta = pi - theta
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            if(is_nan_r8(theta) .or. abs(theta-halfpi) < 1d-8) then
                if(is_nan_r8(theta)) write(666,*) 'delta_eject << sin(theta) =', &
                                        hh / point%r / v, 'corrections applied'
                if(abs(theta-halfpi) < 1d-8) write(666,*) &
                                  'theta is close to pi/2, corrections applied'
                theta = halfpi - 1d-5
                N_of_warnings = N_of_warnings + 1
                if(N_of_warnings > maxNofWarnings) then
                    write(*,*) 'too many warnings have been printed to fort.666'
                    stop
                endif
            endif
!~                     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            call Integrand_delta_ejection(density, u, v, psi, theta, &
                                        point, dphi, source, &
                                        muR, dt, comet, Rast_AU)
            density = density / real(sin(dphi) * point%r * source%r * AU**3)

        end subroutine hc_DUDI_delta_ejection



        ! finds dt -- the time necessary to reach from point r0 to point r
        ! if at point r the particle has speed of v
        ! r and r0 are lengths of position-vectors r and r0
        ! phi is the angle between vectors r and r0
        ! muR is a pseudo gravitational parameter (may be negative)
        subroutine find_delta_t(muR, r, r0, v, dphi, dt, pericenter)
            use twobody_fun
            implicit none
            real(8), intent(in) :: muR, r, r0, v, dphi
            real(8), intent(out) :: dt
            logical, intent(in) :: pericenter
            real(8) ee, theta, semi_major_axis
            
            semi_major_axis = (2d0 / r - v**2 / muR)**(-1)
            dt  = 100d0
            
            if(semi_major_axis > 0d0 .and. muR > 0d0) then
                call theta_geometry_ellipse(muR, r, r0, v, &
                                  dphi, semi_major_axis, ee, &
                                  theta, dt, pericenter)
            endif
            if(semi_major_axis < 0d0 .and. muR > 0d0) then 
                call theta_geometry_hyperbola(muR, r, r0, v, &
                               dphi, abs(semi_major_axis), ee, &
                               theta, dt, pericenter)
            endif
            if(muR < 0d0) then
                call theta_geometry_other_branch(muR, r, r0, v, &
                               dphi, abs(semi_major_axis), ee, &
                               theta, dt, pericenter)
            endif
        
        end subroutine find_delta_t
        
        
        
        
        
        subroutine Qinterpol(vmin, vmax, dt1, dt2, v1, v2, muR, r, r0, dphi, pericenter)
            implicit none
            real(8), intent(in) :: vmin, vmax, dt1, dt2, muR, r, r0, dphi
            logical, intent(in) :: pericenter
            real(8) vmid, a, b, c, t1, t2, t3, tt1, tt2, tt3
            real(8), intent(out) :: v1, v2
            
            vmid = (vmax + vmin) / 2d0
            call find_delta_t(muR, r, r0, vmin, dphi, t1, pericenter)
            call find_delta_t(muR, r, r0, vmax, dphi, t3, pericenter)
            call find_delta_t(muR, r, r0, vmid, dphi, t2, pericenter)
            ! interpolating on the three points with a quadratic polinomial to find
            ! the velocity-boundaries corresponding to the desired interval of dt
            tt3 = t3 * t3
            tt1 = t1 * t1
            tt2 = t2 * t2
            b = (vmin - vmid - (tt1 - tt2) * (vmid - vmax) / (tt2 - tt3)) &
               / ((t1 - t2) - (tt1 - tt2) / (t2 + t3))
            a = ((vmid - vmax) - (t2 - t3) * b) / (tt2 - tt3)
            c = vmin - a * tt1 - b * t1
            v1 = a * dt1**2 + b * dt1 + c
            v2 = a * dt2**2 + b * dt2 + c
        end subroutine Qinterpol
        
        
        
        ! Finds value v corresponding to the given value of dt
        ! using bisection method
        subroutine bisection(dt, v, vmin, vmax, muR, r, r0, dphi, pericenter)
            implicit none
            real(8), parameter :: eps = 1e-12
            real(8), intent(in) :: dt, vmin, vmax, muR, r, r0, dphi
            logical, intent(in) :: pericenter
            real(8), intent(out) :: v
            real(8) v1, dttmp, v2
            integer i
            
            dttmp = dt * 1d3
            v1 = vmin ; v2 = vmax
            v = v1
            i = 0
            do while(abs(dttmp - dt) > eps .and. i < 30)
                i = i + 1
                v = (v1 + v2) / 2d0
                call find_delta_t(muR, r, r0, v, dphi, dttmp, pericenter)
                if(dttmp < dt) then
                    v2 = v
                else
                    v1 = v
                endif
            enddo
        
        end subroutine bisection
        
        
        

        ! Finds the limits of integration (v1, v2)
        ! based on the allowed values of dt
        ! Then performs integration over v
        subroutine hc_DUDI_v_integration(density, point, source, &
                    muR, tnow, comet, Rast_AU, pericenter)
            use const
            use define_types
            use help
            use twobody_fun
            implicit none
            integer i
            real, intent(out) :: density
            real(8), intent(in) :: tnow, muR, Rast_AU
            logical, intent(in) :: pericenter
            type(ephemeris), intent(in) :: comet
            real(8) vmax, vmin, v1, v2, vmid
            type(position_in_space), intent(in) :: point
            type(source_properties), intent(in) :: source
            real(8) term, tmp
            real(8) dt1, dt2
            real(8) dphi, dbeta, vtmp
            real(8) ldif, lsum
            real xi(order_v), wi(order_v)
            
            call ApuTrajectoryLight(point, dphi, dbeta, source)
            tmp = 2d0 * muR * (1d0 / point%r - 1d0 / source%r)
            
            vmin = sqrt((comet%Vast - source%ud%umax)**2 + tmp)
            ! v_max is then the maximal *possible* speed at radius r,
            ! assuming that the ejection velocity is limited by gas velocity
            vmax = sqrt((comet%Vast + source%ud%umax)**2 + tmp)
            vmid = (vmax + vmin) / 2d0
            density = 0.0
            ! if the maximal possible velocity is enough to get from rm to rr
            if(vmax > vmin) then
                dt1 = tnow - source%Tj + source%dtau / 2d0
                dt2 = tnow - source%Tj - source%dtau / 2d0
                call bisection(dt1, v1, vmin, vmax, muR, &
                               point%r, source%r, dphi, pericenter)
                call bisection(dt2, v2, vmin, vmax, muR, &
                               point%r, source%r, dphi, pericenter)
                call GaussLegendreQuadra(xi, wi, order_v)
                ldif = v2 - v1
                if(ldif > 1d-12) then
                    ldif = ldif * 0.5d0
                    lsum = v2 + v1
                    lsum = lsum * 0.5d0
                    do i = 1, order_v
                        vtmp = ldif * xi(i) + lsum
                        call Integrand_v_integration(term, vtmp, &
                                                dt1, point, dphi, &
                                                source, tnow, muR, &
                                                comet, Rast_AU, pericenter)
                        density = density + real(ldif * wi(i) * term)
                    enddo
                    ! factor independent on velocity and conversion to m^-3
                    density = density / real(point%r * source%r * sin(dphi) * AU**3)
                else
                    density = 0.0
                endif
            endif

        end subroutine hc_DUDI_v_integration
        
        
 
 end module DUDIhc
 
 
 
 
 

