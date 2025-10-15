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

! File: help.f90
! Description: Contains functions and subroutines that perform general
! mathematical tasks.

! Author: Anastasiia Ershova
! E-mail: vveyzaa@gmail.com

module help
  implicit none
    contains

      ! rotation of a vector in the plain
      pure function rot2d(r, phi) result(x)
        implicit none
        real(8), intent(in) :: r(2), phi
        real(8) x(2)
        
        x(1) = r(1) * cos(phi) - r(2) * sin(phi)
        x(2) = r(1) * sin(phi) + r(2) * cos(phi)
        
      end function rot2d
      
      
    ! returns the position vector r computed for a point mass that had
    ! a position vector `r0´ and a velocity vector `v0´ `time´ ago
    ! the trajectory is integrated by the 4th order  Runge-Kutta method
    ! within two-body problem with the gravitational parameter `mu´
        subroutine runge_kutta_point_position(r0, v0, mu, time, r)
            implicit none
            real(8), intent(in) :: r0(3), v0(3), mu, time
            real(8), intent(out) :: r(3)
            integer i
            integer Nstep
            real(8) dt, v(3)
            real(8) k1(3), k2(3), k3(3), k4(3), l1(3), l2(3), l3(3), l4(3)
            
            r = r0 ; v = v0
            
            Nstep = 200
            dt = time / dble(Nstep)
            do while(dt > 3d-4)
                Nstep = int(Nstep * 1.2)
                dt = time / dble(Nstep)
            enddo
            do i = 1, Nstep
                k1 = -mu / sqrt(sum(r**2))**3 * r
                l1 = v
                k2 = -mu / sqrt(sum((r + l1 * dt/2d0)**2))**3 * (r + l1 * dt/2d0)
                l2 = v + k1 * dt/2d0
                k3 = -mu / sqrt(sum((r + l2 * dt/2d0)**2))**3 * (r + l2 * dt/2d0)
                l3 = v + k2 * dt/2d0
                k4 = -mu / sqrt(sum((r + l3 * dt)**2))**3 * (r + l3 * dt)
                l4 = v + k3 * dt
                v = v + dt / 6d0 * (k1 + 2d0 * k2 + 2d0 * k3 + k4)
                r = r + dt / 6d0 * (l1 + 2d0 * l2 + 2d0 * l3 + l4)
            enddo
            
        end subroutine runge_kutta_point_position      
      
      
      
      ! solves a system of 2 quadratic equations
      ! (x - x0)**2 + (y - y0)**2 = R0**2, x0**2 + y0**2 = 1
      !(x - x1)**2 + y**2 = R1**2
      ! returns 2-d vectors x and y where (x(1), y(1))
      ! and (x(2), y(2)) are solutions of the SoE
      pure subroutine circle_intersection(x0, y0, R0, x1, R1, x, y)
        implicit none
        real(8), intent(in) :: x0, y0, R0, x1, R1
        real(8), intent(out) :: x(2), y(2)
        real(8) sR, dx, sx, dx2, dy2, R02, R12, x1x0
        real(8) sumdifs2, sqrtshort, ybracket
        real(8) tmp(2), eps, eps1, eps0
        
        sR = R0 + R1
        dx = x0 - x1 ; sx = x0 + x1
        x1x0 = x0 * x1
        dx2 = dx * dx ; dy2 = y0 * y0
        R02 = R0 * R0 ; R12 = R1 * R1
        
        sumdifs2 = dx2 + dy2
        
        sqrtshort = sqrt(-2d0 * (x1x0 - abs(x1)) * (sR**2 - sumdifs2) * dy2)
        
        ybracket = (R12 + sumdifs2) * dy2

        x(1) = ((R12 - R02) * dx + sx * sumdifs2 &
          - sqrtshort) / (2d0 * sumdifs2)


        y(1) = (-R02 * dy2 + (x0 - x1) * sqrtshort + ybracket) &
            / (2d0 * sumdifs2 * y0)

        x(2) = ((R12 - R02) * dx + sx * sumdifs2 &
          + sqrtshort) / (2d0 * sumdifs2)


        y(2) = (-R02 * dy2 + (x1 - x0) * sqrtshort + ybracket) &
            / (2d0 * sumdifs2 * y0)
        
      end subroutine circle_intersection
      
      
      ! transforms Cartesian coordinates `x´, `y´, and `z´
      ! to cylindrical radius `r´, longitude `phi´, and height `h´
      pure subroutine cartesian2cylindrical(x, y, z, r, phi, h)
        implicit none
        real(8), intent(in) :: x, y, z
        real(8), intent(out) :: r, phi, h
        
        r = sqrt(x * x + y * y)
        phi = atan(y, x)
        h = z
        
      end subroutine cartesian2cylindrical
      
    
      ! cross product of two 3d vectors `x´ and `y´
      pure function vector_product(x,y) result(z)
        implicit none
        real(8), intent(in) :: x(3), y(3)
        real(8) z(3)
        
        z(1) = x(2) * y(3) - x(3) * y(2)
        z(2) = x(3) * y(1) - x(1) * y(3)
        z(3) = x(1) * y(2) - x(2) * y(1)
      
      end function vector_product
    
  
      
    
      ! Linear interpolation. Returns a value of `y´ at `xout´
      ! x must be in ascending order
      ! if `xout´ < x(1) or xout > x(N) then 
      ! `yout´ = y(1) or `yout´ = y(N) respectively
      pure function LiNTERPOL(N, y, x, xout) result(yout)
        implicit none
        real(8), intent(in) ::  xout, y(N), x(N)
        real(8) yout
        real(8) :: x1 = 1d0
        real(8) :: x2 = 1d0
        integer i, i1
        integer :: i2 = -1
        integer, intent(in) :: N
        
        if(xout < x(1)) then
          yout = y(1)
          return
        endif
        
        if(xout > x(N)) then
          yout = y(N)
          return
        endif
        
        i = 1
        i1 = 0
        do while(i1 == 0 .and. i .lt. N) 
          if(xout .lt. x(i+1) .and. xout .ge. x(i)) then
            i1 = i;       i2 = i+1;
            x1 = x(i1);   x2 = x(i2);
          endif
          i = i + 1
        enddo

        yout = y(i1) + (y(i2) - y(i1)) / (x2 - x1) * (xout - x1)
        return
        
      end function LiNTERPOL  
  
  
  
      ! Euclidean norm of a 3d vector
      pure function norma3d(v)
        implicit none
        real(8), intent(in) :: v(3)
        real(8) norma3d
        
          norma3d = sqrt(dot_product(v,v))
          
      end function norma3d
      
      
      ! Euclidean norm of a 2d vector
      pure function norma2d(v)
        implicit none
        real(8), intent(in) :: v(2)
        real(8) norma2d
        
          norma2d = sqrt(dot_product(v,v))
          
      end function norma2d
        
  
      
      
    ! modified inverse tanget function 
    ! returns values in [0, 2pi]      
      pure function myatan(N, re0, im)
        use const
        implicit none
        integer i
        integer, intent(in):: N
        real(8), intent(in) :: re0(N), im(N)
        real(8) myatan(N), re(N)
        
        re = re0

      ! avoid problems with zero re
        do i = 1, N
          if(re(i) < 1d-12 .and. re(i) > 0d0) re(i) = 1d-12 
          if(re(i) > 1d-12 .and. re(i) < 0d0) re(i) = -1d-12 
        enddo

        myatan = atan(im / re)

        do i = 1, N
          if(re(i) < 0d0) myatan(i) = myatan(i) + pi
          if(re(i) > 0d0 .and. im(i) < 0d0) myatan(i) = myatan(i) + twopi
        enddo

      end function myatan
      

      
      ! Euler's rotation
      pure subroutine eulrot(phiE, thetaE, psiE, xin, yin, zin, xout, yout, zout, inverse)
        implicit none
        integer i
        logical, intent(in) :: inverse
        real(8), intent(in) :: phiE, thetaE, psiE, xin, yin, zin
        real(8), intent(out) :: xout, yout, zout
        real(8) cp, sp, ct, st, cps, sps
        
        cp = cos(phiE)
        sp = sin(phiE)
        ct = cos(thetaE)
        st = sin(thetaE)
        cps = cos(psiE)
        sps = sin(psiE)

        if (.not. inverse) then      ! if NOT INVERSE case
        
          xout= (cps*cp-ct*sp*sps)*xin + (cps*sp+ct*cp*sps)*yin + sps*st*zin
          yout=(-sps*cp-ct*sp*cps)*xin + (-sps*sp+ct*cp*cps)*yin + cps*st*zin
          zout=st*sp*xin - st*cp*yin + ct*zin  

        else 

          xout=(cps*cp-ct*sp*sps)*xin + (-sps*cp-ct*sp*cps)*yin + st*sp*zin
          yout=(cps*sp+ct*cp*sps)*xin + (-sps*sp+ct*cp*cps)*yin - st*cp*zin
          zout=st*sps*xin + st*cps*yin + ct*zin

        endif

      end subroutine eulrot


      ! for given coordinates of the source and the spacecraft compute:
      ! dphi (\Delta\phi) the angle between
      ! the moon-centered radius-vectors of the source and the spacecraft
      ! dbeta (\Delta\beta) the angle between the projections of these vectors
      ! to the moon equatorial plane
      pure subroutine ApuTrajectoryLight(point, dphi, dbeta, source)
        use const
        use define_types
        use nan_utils
        implicit none
        integer i
        type(position_in_space), intent(in) :: point
        type(source_properties), intent(in) :: source
        real(8), intent(out) :: dphi, dbeta
        real(8) Rsource(3)
        
        Rsource = source%rrM / source%r
        
        ! dphi is an angle between a vector pointing to the source
        ! and a vector pointing to the spacecraft
        dphi = acos(dot_product(Rsource, point%rvector) / point%r)
        if(is_nan_r8(dphi) .or. dphi < 2d-8) &
            dphi = sign(1d0, dot_product(Rsource, point%rvector)) * 1d-8
        
        ! dbeta is an angle between proections of the same vectors
        ! in the longitudinal plane
        dbeta = acos((Rsource(1) * point%rvector(1) &
            + Rsource(2) * point%rvector(2)) &
            / (norma2d(Rsource(1:2)) * norma2d(point%rvector(1:2))))
        if(is_nan_r8(dbeta) .or. dbeta < 1d-8) &
        dbeta = sign(1d0, dot_product(Rsource, point%rvector)) * 1d-8
        
      end subroutine ApuTrajectoryLight
      
      
      ! returns azimuth of the vector direction at the point location
      ! the azimuth is counted CLOCKWISE from the local north
      function azimuth(direction, location)
        use const
        real(8), intent(in) :: direction(3), location(3)
        real(8), parameter, dimension(3) :: zvec = (/0d0, 0d0, 1d0/)
        real(8) azimuth
        real(8) tmpr(3), north(3), dirplane(3), testhir(3), tmpd(3)
        real(8) dotpr
        
        ! unit vector normal to surface of a sphere
        ! passing through the point location
        tmpr = location / norma3d(location)
        
        ! the component of the vector along the Z-axis
        ! which is tangential to the sphere
        north = zvec - tmpr * dot_product(zvec, tmpr)
        
        ! normalized
        north = north / norma3d(north)
        
        ! the component of the vector direction
        ! which is tangential to the sphere
        tmpd = direction / norma3d(direction)
        dirplane = tmpd - tmpr * dot_product(tmpd, tmpr)
        ! normalized
        dirplane = dirplane / norma3d(dirplane)
        
        ! angle between the vector pointing to the local north
        ! and the horizontal component of direction
        dotpr = dot_product(dirplane, north)
        
        if(abs(dotpr) > 1d0) dotpr = sign(1d0, dotpr)
        azimuth = acos(dotpr)
        
        ! check on which side of the vector north
        ! lays the vector dirplane by comparing
        ! their vector product to the local normal to the surface
        ! (pointing outwards)
        testhir = vector_product(north, dirplane)
        testhir = testhir / norma3d(testhir)
        if(norma3d(testhir - tmpr) < 1d-3) azimuth = twopi - azimuth
        
      end function azimuth
      
      
            ! solves a cubic equation a(0)*x^3 + a(1) * x^2 + a(2)*x + a(3)
      subroutine cardano_formula(a, x)
        use const
        implicit none
        real(8), intent(in) :: a(0:3)
        integer, parameter :: dk = kind(0.0d0)
        complex(kind=dk), intent(out) :: x(3)
        complex(kind=dk) :: omega1, omega2, y(3), roots1(3), roots2(3)
        integer k, i, ii
        real(8) B1, B2, B3, pp, qq, phi1, phi2
        
        B1 = a(1) / a(0)
        B2 = a(2) / a(0)
        B3 = a(3) / a(0)
        
        pp = -B1**2 / 3d0 + B2
        qq = 2d0 * B1**3 / 27d0 - B1 * B2 / 3d0 + B3
        
        omega1 = -qq / 2d0 + sqrt( cmplx(qq**2 / 4d0 + pp**3 / 27d0, kind=dk) )
        omega2 = -qq / 2d0 - sqrt( cmplx(qq**2 / 4d0 + pp**3 / 27d0, kind=dk) )
        
        phi1 = atan(imag(omega1), real(omega1))
        phi2 = atan(imag(omega2), real(omega2))
        do k = 1, 3
          roots1(k) = abs(omega1)**(1.0d0/3.0d0) * cmplx( cos((phi1 + twopi * (k-1)) / 3d0), &
                       sin((phi1 + twopi * (k-1)) / 3d0), kind=dk )

          roots2(k) = abs(omega2)**(1.0d0/3.0d0) * cmplx( cos((phi2 + twopi * (k-1)) / 3d0), &
                      sin((phi2 + twopi * (k-1)) / 3d0), kind=dk )

        enddo
        
        k = 1
        do i = 1, 3
          do ii = 1, 3
            if(k <= 3 .and. abs(roots1(i) * roots2(ii) + pp / 3d0) < 1d-8) then
              y(k) = roots1(i) + roots2(ii)
              x(k) = y(k) - B1 / 3d0
              k = k + 1
            endif
          enddo
        enddo
      
      end subroutine cardano_formula

            
            ! computes the invert matrix 3x3
      pure subroutine invert_matrix3(A, B)
        implicit none
        real(8), intent(in) :: A(3,3)
        real(8), intent(out) :: B(3,3)
        real(8) detinv

        detinv = (A(1,1) * A(2,2) * A(3,3) - A(1,1) * A(2,3)  * A(3,2) &
            - A(1,2) * A(2,1) * A(3,3) + A(1,2) * A(2,3) * A(3,1) &
            + A(1,3) * A(2,1) * A(3,2) - A(1,3) * A(2,2) * A(3,1))
        detinv = 1d0 / detinv

        B(1,1) = detinv * (A(2,2) * A(3,3) - A(2,3) * A(3,2))
        B(2,1) = - detinv * (A(2,1) * A(3,3) - A(2,3) * A(3,1))
        B(3,1) = detinv * (A(2,1) * A(3,2) - A(2,2) * A(3,1))
        B(1,2) = - detinv * (A(1,2) * A(3,3) - A(1,3) * A(3,2))
        B(2,2) = detinv * (A(1,1) * A(3,3) - A(1,3) * A(3,1))
        B(3,2) = - detinv * (A(1,1) * A(3,2) - A(1,2) * A(3,1))
        B(1,3) = detinv * (A(1,2) * A(2,3) - A(1,3) * A(2,2))
        B(2,3) = - detinv * (A(1,1) * A(2,3) - A(1,3) * A(2,1))
        B(3,3) = detinv * (A(1,1) * A(2,2) - A(1,2) * A(2,1))
        
      end subroutine invert_matrix3


            ! computes the invert matrix 4x4
      pure function invert_matrix4(A) result(B)
        real(8), intent(in) :: A(4,4)   
        real(8) B(4,4)   
        real(8) detinv

        detinv = &
          1/(A(1,1)*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))&
           - A(1,2)*(A(2,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))&
           + A(1,3)*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))&
           - A(1,4)*(A(2,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(2,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(2,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1))))

        B(1,1) = detinv*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))
        B(2,1) = detinv*(A(2,1)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(2,3)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(2,4)*(A(3,3)*A(4,1)-A(3,1)*A(4,3)))
        B(3,1) = detinv*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
        B(4,1) = detinv*(A(2,1)*(A(3,3)*A(4,2)-A(3,2)*A(4,3))+A(2,2)*(A(3,1)*A(4,3)-A(3,3)*A(4,1))+A(2,3)*(A(3,2)*A(4,1)-A(3,1)*A(4,2)))
        B(1,2) = detinv*(A(1,2)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(1,3)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(1,4)*(A(3,3)*A(4,2)-A(3,2)*A(4,3)))
        B(2,2) = detinv*(A(1,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(1,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(1,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))
        B(3,2) = detinv*(A(1,1)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(1,2)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(1,4)*(A(3,2)*A(4,1)-A(3,1)*A(4,2)))
        B(4,2) = detinv*(A(1,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(1,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(1,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
        B(1,3) = detinv*(A(1,2)*(A(2,3)*A(4,4)-A(2,4)*A(4,3))+A(1,3)*(A(2,4)*A(4,2)-A(2,2)*A(4,4))+A(1,4)*(A(2,2)*A(4,3)-A(2,3)*A(4,2)))
        B(2,3) = detinv*(A(1,1)*(A(2,4)*A(4,3)-A(2,3)*A(4,4))+A(1,3)*(A(2,1)*A(4,4)-A(2,4)*A(4,1))+A(1,4)*(A(2,3)*A(4,1)-A(2,1)*A(4,3)))
        B(3,3) = detinv*(A(1,1)*(A(2,2)*A(4,4)-A(2,4)*A(4,2))+A(1,2)*(A(2,4)*A(4,1)-A(2,1)*A(4,4))+A(1,4)*(A(2,1)*A(4,2)-A(2,2)*A(4,1)))
        B(4,3) = detinv*(A(1,1)*(A(2,3)*A(4,2)-A(2,2)*A(4,3))+A(1,2)*(A(2,1)*A(4,3)-A(2,3)*A(4,1))+A(1,3)*(A(2,2)*A(4,1)-A(2,1)*A(4,2)))
        B(1,4) = detinv*(A(1,2)*(A(2,4)*A(3,3)-A(2,3)*A(3,4))+A(1,3)*(A(2,2)*A(3,4)-A(2,4)*A(3,2))+A(1,4)*(A(2,3)*A(3,2)-A(2,2)*A(3,3)))
        B(2,4) = detinv*(A(1,1)*(A(2,3)*A(3,4)-A(2,4)*A(3,3))+A(1,3)*(A(2,4)*A(3,1)-A(2,1)*A(3,4))+A(1,4)*(A(2,1)*A(3,3)-A(2,3)*A(3,1)))
        B(3,4) = detinv*(A(1,1)*(A(2,4)*A(3,2)-A(2,2)*A(3,4))+A(1,2)*(A(2,1)*A(3,4)-A(2,4)*A(3,1))+A(1,4)*(A(2,2)*A(3,1)-A(2,1)*A(3,2)))
        B(4,4) = detinv*(A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))+A(1,2)*(A(2,3)*A(3,1)-A(2,1)*A(3,3))+A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1)))
        end function



      
      
      ! returns nodes and weights of the Gauss-Legendre quadrature formula
      subroutine GaussLegendreQuadra(xi, wi, order)
        use const
        implicit none
        integer, intent(in) :: order 
        real, intent(out) :: xi(order), wi(order)
        
        if(order /= 5 .and. order /= 10 &
              .and. order /= 20 .and. order /= 30 .and. order /= 40 &
              .and. order /= 50 .and. order /= 64) then
          write(*,*) 'order of integration must be equal to 5, 10, 20, 30, 40, 50 or 64'
          write(*,*) 'check variables order_v and order_R'
          stop
        endif
        
        select case(order)
          case(64)
             wi(1) = 0.0486909570091397;    xi(1) = -0.0243502926634244
             wi(2) = 0.0486909570091397;    xi(2) = 0.0243502926634244
             wi(3) = 0.0485754674415034;    xi(3) = -0.0729931217877990
             wi(4) = 0.0485754674415034;    xi(4) = 0.0729931217877990
             wi(5) = 0.0483447622348030;    xi(5) = -0.1214628192961206
             wi(6) = 0.0483447622348030;    xi(6) = 0.1214628192961206
             wi(7) = 0.0479993885964583;    xi(7) = -0.1696444204239928
             wi(8) = 0.0479993885964583;    xi(8) = 0.1696444204239928
             wi(9) = 0.0475401657148303;    xi(9) = -0.2174236437400071
             wi(10) = 0.0475401657148303;    xi(10) = 0.2174236437400071
             wi(11) = 0.0469681828162100;    xi(11) = -0.2646871622087674
             wi(12) = 0.0469681828162100;    xi(12) = 0.2646871622087674
             wi(13) = 0.0462847965813144;    xi(13) = -0.3113228719902110
             wi(14) = 0.0462847965813144;    xi(14) = 0.3113228719902110
             wi(15) = 0.0454916279274181;    xi(15) = -0.3572201583376681
             wi(16) = 0.0454916279274181;    xi(16) = 0.3572201583376681
             wi(17) = 0.0445905581637566;    xi(17) = -0.4022701579639916
             wi(18) = 0.0445905581637566;    xi(18) = 0.4022701579639916
             wi(19) = 0.0435837245293235;    xi(19) = -0.4463660172534641
             wi(20) = 0.0435837245293235;    xi(20) = 0.4463660172534641
             wi(21) = 0.0424735151236536;    xi(21) = -0.4894031457070530
             wi(22) = 0.0424735151236536;    xi(22) = 0.4894031457070530
             wi(23) = 0.0412625632426235;    xi(23) = -0.5312794640198946
             wi(24) = 0.0412625632426235;    xi(24) = 0.5312794640198946
             wi(25) = 0.0399537411327203;    xi(25) = -0.5718956462026340
             wi(26) = 0.0399537411327203;    xi(26) = 0.5718956462026340
             wi(27) = 0.0385501531786156;    xi(27) = -0.6111553551723933
             wi(28) = 0.0385501531786156;    xi(28) = 0.6111553551723933
             wi(29) = 0.0370551285402400;    xi(29) = -0.6489654712546573
             wi(30) = 0.0370551285402400;    xi(30) = 0.6489654712546573
             wi(31) = 0.0354722132568824;    xi(31) = -0.6852363130542333
             wi(32) = 0.0354722132568824;    xi(32) = 0.6852363130542333
             wi(33) = 0.0338051618371416;    xi(33) = -0.7198818501716109
             wi(34) = 0.0338051618371416;    xi(34) = 0.7198818501716109
             wi(35) = 0.0320579283548516;    xi(35) = -0.7528199072605319
             wi(36) = 0.0320579283548516;    xi(36) = 0.7528199072605319
             wi(37) = 0.0302346570724025;    xi(37) = -0.7839723589433414
             wi(38) = 0.0302346570724025;    xi(38) = 0.7839723589433414
             wi(39) = 0.0283396726142595;    xi(39) = -0.8132653151227975
             wi(40) = 0.0283396726142595;    xi(40) = 0.8132653151227975
             wi(41) = 0.0263774697150547;    xi(41) = -0.8406292962525803
             wi(42) = 0.0263774697150547;    xi(42) = 0.8406292962525803
             wi(43) = 0.0243527025687109;    xi(43) = -0.8659993981540928
             wi(44) = 0.0243527025687109;    xi(44) = 0.8659993981540928
             wi(45) = 0.0222701738083833;    xi(45) = -0.8893154459951141
             wi(46) = 0.0222701738083833;    xi(46) = 0.8893154459951141
             wi(47) = 0.0201348231535302;    xi(47) = -0.9105221370785028
             wi(48) = 0.0201348231535302;    xi(48) = 0.9105221370785028
             wi(49) = 0.0179517157756973;    xi(49) = -0.9295691721319396
             wi(50) = 0.0179517157756973;    xi(50) = 0.9295691721319396
             wi(51) = 0.0157260304760247;    xi(51) = -0.9464113748584028
             wi(52) = 0.0157260304760247;    xi(52) = 0.9464113748584028
             wi(53) = 0.0134630478967186;    xi(53) = -0.9610087996520538
             wi(54) = 0.0134630478967186;    xi(54) = 0.9610087996520538
             wi(55) = 0.0111681394601311;    xi(55) = -0.9733268277899110
             wi(56) = 0.0111681394601311;    xi(56) = 0.9733268277899110
             wi(57) = 0.0088467598263639;    xi(57) = -0.9833362538846260
             wi(58) = 0.0088467598263639;    xi(58) = 0.9833362538846260
             wi(59) = 0.0065044579689784;    xi(59) = -0.9910133714767443
             wi(60) = 0.0065044579689784;    xi(60) = 0.9910133714767443
             wi(61) = 0.0041470332605625;    xi(61) = -0.9963401167719553
             wi(62) = 0.0041470332605625;    xi(62) = 0.9963401167719553
             wi(63) = 0.0017832807216964;    xi(63) = -0.9993050417357722
             wi(64) = 0.0017832807216964;    xi(64) = 0.9993050417357722
          case(50)
             xi(1) = -0.9988664044;    wi(1) = 0.0029086226
             xi(2) = -0.9940319694;    wi(2) = 0.0067597992
             xi(3) = -0.9853540841;    wi(3) = 0.0105905484
             xi(4) = -0.9728643851;    wi(4) = 0.0143808228
             xi(5) = -0.9566109552;    wi(5) = 0.0181155607
             xi(6) = -0.9366566189;    wi(6) = 0.0217802432
             xi(7) = -0.9130785567;    wi(7) = 0.0253606736
             xi(8) = -0.8859679795;    wi(8) = 0.0288429936
             xi(9) = -0.8554297694;    wi(9) = 0.032213728
             xi(10) = -0.8215820709;   wi(10) = 0.035459836
             xi(11) = -0.7845558329;   wi(11) = 0.0385687566
             xi(12) = -0.7444943022;   wi(12) = 0.0415284631
             xi(13) = -0.7015524687;   wi(13) = 0.0443275043
             xi(14) = -0.6558964657;   wi(14) = 0.0469550513
             xi(15) = -0.6077029272;   wi(15) = 0.0494009384
             xi(16) = -0.5571583045;   wi(16) = 0.0516557031
             xi(17) = -0.5044581449;   wi(17) = 0.0537106219
             xi(18) = -0.449806335;    wi(18) = 0.0555577448
             xi(19) = -0.393414312;    wi(19) = 0.05718992565
             xi(20) = -0.335500245;    wi(20) = 0.05860084981
             xi(21) = -0.2762881938;   wi(21) = 0.0597850587
             xi(22) = -0.2160072369;   wi(22) = 0.0607379708
             xi(23) = -0.15489059;     wi(23) = 0.0614558996
             xi(24) = -0.0931747016;   wi(24) = 0.0619360674
             xi(25) = -0.0310983383;   wi(25) = 0.0621766167
             xi(26) = 0.0310983383;    wi(26) = 0.0621766167
             xi(27) = 0.0931747016;    wi(27) = 0.0619360674
             xi(28) = 0.15489059 ;     wi(28) = 0.0614558996
             xi(29) = 0.2160072369;    wi(29) = 0.0607379708
             xi(30) = 0.2762881938;    wi(30) = 0.0597850587
             xi(31) = 0.335500245;     wi(31) = 0.0586008498
             xi(32) = 0.3934143119;    wi(32) = 0.057189926
             xi(33) = 0.449806335;     wi(33) = 0.0555577448
             xi(34) = 0.504458145;     wi(34) = 0.0537106219
             xi(35) = 0.5571583045;    wi(35) = 0.051655703
             xi(36) = 0.6077029272;    wi(36) = 0.04940093845
             xi(37) = 0.6558964657;    wi(37) = 0.0469550513
             xi(38) = 0.7015524687;    wi(38) = 0.044327504
             xi(39) = 0.7444943022;    wi(39) = 0.0415284631
             xi(40) = 0.7845558329;    wi(40) = 0.03856875661
             xi(41) = 0.8215820709;    wi(41) = 0.0354598356
             xi(42) = 0.8554297694;    wi(42) = 0.0322137282
             xi(43) = 0.8859679795;    wi(43) = 0.0288429936
             xi(44) = 0.9130785567;    wi(44) = 0.025360674
             xi(45) = 0.9366566189;    wi(45) = 0.021780243
             xi(46) = 0.9566109552;    wi(46) = 0.018115561
             xi(47) = 0.9728643851;    wi(47) = 0.0143808228
             xi(48) = 0.9853540841;    wi(48) = 0.0105905484
             xi(49) = 0.9940319694;    wi(49) = 0.006759799
             xi(50) = 0.9988664044;    wi(50) = 0.00290862255
          case(40)
             xi(1) = -0.998237710;     wi(1) = 0.0045212771
             xi(2) = -0.990726239;     wi(2) = 0.0104982845
             xi(3) = -0.977259950;     wi(3) = 0.0164210584
             xi(4) = -0.957916819;     wi(4) = 0.0222458492
             xi(5) = -0.932812808;     wi(5) = 0.0279370070
             xi(6) = -0.902098807;     wi(6) = 0.0334601953
             xi(7) = -0.865959503;     wi(7) = 0.0387821680
             xi(8) = -0.824612231;     wi(8) = 0.0438709082
             xi(9) = -0.778305651;     wi(9) = 0.0486958076
             xi(10) = -0.727318255;     wi(10) = 0.053227847
             xi(11) = -0.671956685;     wi(11) = 0.057439769
             xi(12) = -0.612553890;     wi(12) = 0.061306242
             xi(13) = -0.549467125;     wi(13) = 0.064804013
             xi(14) = -0.483075802;     wi(14) = 0.067912046
             xi(15) = -0.413779204;     wi(15) = 0.070611647
             xi(16) = -0.341994091;     wi(16) = 0.072886582
             xi(17) = -0.268152185;     wi(17) = 0.074723169
             xi(18) = -0.192697581;     wi(18) = 0.076110362
             xi(19) = -0.116084071;     wi(19) = 0.077039818
             xi(20) = -0.038772418;     wi(20) = 0.077505948
             xi(21) = 0.0387724175;     wi(21) = 0.077505948
             xi(22) = 0.1160840707;     wi(22) = 0.077039818
             xi(23) = 0.1926975807;     wi(23) = 0.076110362
             xi(24) = 0.2681521850;     wi(24) = 0.074723169
             xi(25) = 0.3419940908;     wi(25) = 0.072886582
             xi(26) = 0.4137792044;     wi(26) = 0.070611647
             xi(27) = 0.4830758017;     wi(27) = 0.067912046
             xi(28) = 0.5494671251;     wi(28) = 0.064804013
             xi(29) = 0.6125538897;     wi(29) = 0.061306242
             xi(30) = 0.6719566846;     wi(30) = 0.057439769
             xi(31) = 0.7273182552;     wi(31) = 0.053227847
             xi(32) = 0.7783056514;     wi(32) = 0.048695808
             xi(33) = 0.8246122308;     wi(33) = 0.043870908
             xi(34) = 0.8659595032;     wi(34) = 0.038782168
             xi(35) = 0.9020988070;     wi(35) = 0.033460195
             xi(36) = 0.9328128083;     wi(36) = 0.027937007
             xi(37) = 0.9579168192;     wi(37) = 0.022245849
             xi(38) = 0.9772599500;     wi(38) = 0.016421058
             xi(39) = 0.9907262387;     wi(39) = 0.010498285
             xi(40) = 0.9982377097;     wi(40) = 0.004521277
          case(30)
             xi(1) = -0.996893466;     wi(1) = 7.96819292E-03
             xi(2) = -0.983668149;     wi(2) = 1.84664689E-02
             xi(3) = -0.960021853;     wi(3) = 2.87847072E-02
             xi(4) = -0.926200032;     wi(4) = 3.87991928E-02
             xi(5) = -0.882560551;     wi(5) = 4.84026745E-02
             xi(6) = -0.829565763;     wi(6) = 5.74931577E-02
             xi(7) = -0.767777443;     wi(7) = 6.59742281E-02
             xi(8) = -0.697850466;     wi(8) = 7.37559721E-02
             xi(9) = -0.620526195;     wi(9) = 8.07558969E-02
             xi(10) = -0.536624134;    wi(10) = 8.68997872E-02
             xi(11) = -0.447033763;    wi(11) = 9.21225250E-02
             xi(12) = -0.352704734;    wi(12) = 9.63687375E-02
             xi(13) = -0.254636914;    wi(13) = 9.95934233E-02
             xi(14) = -0.153869912;    wi(14) = 0.101762392    
             xi(15) = -5.14718443E-02; wi(15) = 0.102852650    
             xi(16) = 5.14718443E-02;  wi(16) = 0.102852650    
             xi(17) = 0.153869912;     wi(17) = 0.101762392    
             xi(18) = 0.254636914;     wi(18) = 9.95934233E-02
             xi(19) = 0.352704734;     wi(19) = 9.63687375E-02
             xi(20) = 0.447033763;     wi(20) = 9.21225250E-02
             xi(21) = 0.536624134;     wi(21) = 8.68997872E-02
             xi(22) = 0.620526195;     wi(22) = 8.07558969E-02
             xi(23) = 0.697850466;     wi(23) = 7.37559721E-02
             xi(24) = 0.767777443;     wi(24) = 6.59742281E-02
             xi(25) = 0.829565763;     wi(25) = 5.74931577E-02
             xi(26) = 0.882560551;     wi(26) = 4.84026745E-02
             xi(27) = 0.926200032;     wi(27) = 3.87991928E-02
             xi(28) = 0.960021853;     wi(28) = 2.87847072E-02
             xi(29) = 0.983668149;     wi(29) = 1.84664689E-02
             xi(30) = 0.996893466;     wi(30) = 7.96819292E-03
          case(20)
            xi(1) = -0.993128598;     wi(1) = 1.76140070E-02
            xi(2) = -0.963971913;     wi(2) = 4.06014286E-02
            xi(3) = -0.912234426;     wi(3) = 6.26720488E-02
            xi(4) = -0.839116991;     wi(4) = 8.32767412E-02
            xi(5) = -0.746331930;     wi(5) = 0.101930119    
            xi(6) = -0.636053681;     wi(6) = 0.118194535    
            xi(7) = -0.510867000;     wi(7) = 0.131688640    
            xi(8) = -0.373706102;     wi(8) = 0.142096102    
            xi(9) = -0.227785856;     wi(9) = 0.149172992    
            xi(10) = -7.65265226E-02; wi(10) = 0.152753383    
            xi(11) = 7.65265226E-02;  wi(11) = 0.152753383    
            xi(12) = 0.227785856;     wi(12) = 0.149172992    
            xi(13) = 0.373706102;     wi(13) = 0.142096102    
            xi(14) = 0.510867000;     wi(14) = 0.131688640    
            xi(15) = 0.636053681;     wi(15) = 0.118194535    
            xi(16) = 0.746331930;     wi(16) = 0.101930119    
            xi(17) = 0.839116991;     wi(17) = 8.32767412E-02
            xi(18) = 0.912234426;     wi(18) = 6.26720488E-02
            xi(19) = 0.963971913;     wi(19) = 4.06014286E-02
            xi(20) = 0.993128598;     wi(20) = 1.76140070E-02
          case(10)
            xi(1) = -0.973906517;     wi(1) = 6.66713417E-02
            xi(2) = -0.865063369;     wi(2) = 0.149451345    
            xi(3) = -0.679409564;     wi(3) = 0.219086364    
            xi(4) = -0.433395386;     wi(4) = 0.269266725    
            xi(5) = -0.148874342;     wi(5) = 0.295524240    
            xi(6) = 0.148874342;      wi(6) = 0.295524240    
            xi(7) = 0.433395386;      wi(7) = 0.269266725    
            xi(8) = 0.679409564;      wi(8) = 0.219086364    
             xi(9) = 0.865063369;      wi(9) = 0.149451345    
            xi(10) = 0.973906517;     wi(10) = 6.66713417E-02
          case(5)
            xi(1) = -0.906179845;    wi(1) = 0.236926883    
            xi(2) = -0.538469315;    wi(2) = 0.478628665    
            xi(3) = 0.00000000;      wi(3) = 0.568888903    
            xi(4) = 0.538469315;     wi(4) = 0.478628665    
            xi(5) = 0.906179845;     wi(5) = 0.236926883 
        
        endselect
      
      
      end subroutine GaussLegendreQuadra
      
      

end module help
