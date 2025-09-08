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

! File: distributions_fun.f90
! Description: Contains unctions that describe the ejection process and
! an auxiliary function supporting these calculations.

! Author: Anastasiia Ershova
! E-mail: vveyzaa@gmail.com

module distributions_fun
	use const
	implicit none
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! parameters and variables for modeling the distribution of 
	! the ejection direction using the impact-ejecta density maps
	! from Szalay et al., 2019
	real(8) lonmax, lonmin
	integer, parameter :: nlats = 90
	integer, parameter :: nlons = 180
	real lats(nlats), lons(nlons)
	real(8) rmap1(nlons, nlats), rmap2(nlons, nlats), ratemap(nlons, nlats)
	real(8) rMtmp(3)
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        contains
            
            subroutine ratematr_interpolate(rhel, rhel1, rhel2)
                implicit none
                real(8), intent(in) :: rhel, rhel1, rhel2
                
                ratemap = rmap1 + (rhel - rhel1) / (rhel2 - rhel1) * (rmap2 - rmap1)
                
            
            end subroutine ratematr_interpolate
            
            
            ! reads a matrix from Jamey Szalay's file fname
            ! rhel is the heliocentric distance of the asteroid
            ! corresponding to the moments the map is computed for
            ! lats and lons are the latitudes and longitudes to which
            ! the columns and the rows of the matrix correspond
            ! they are not exactly uniform
            ! the subroutine will return the lats and lons in radians
            ! and rhel in AU
            subroutine read_ratemap(fname, rhel)
                character(*), intent(in) :: fname
                real(8), intent(out) :: rhel
                integer i
                character(len = 18) strtmp
                character(len = 14) strtmp2
                real(8) x, y, z
                
                open(100, file = fname, status = 'old')
                    do i = 1, 5
                        read(100,*) strtmp
                    enddo
                    read(100,*) strtmp, x, y, z
                    rhel = sqrt(sum((/x, y, z/)**2)) * 1d3        ! in meters
                    rhel = rhel / AU                        ! in AU
                    read(100,*) strtmp
                    read(100,*) strtmp2, lons
                    read(100,*) strtmp2, lats
                    do i = 1, 3
                        read(100,*) strtmp
                    enddo
                    do i = 1, nlats
                        read(100,*) rmap2(:,i)
                    enddo
                close(100)
                lats = lats * deg2rad
                lons = lons * deg2rad
                
                lonmin = minval(lons)
                lonmax = maxval(lons)
            
            end subroutine read_ratemap
            
            ! reads a matrix from Jamey Szalay's file fname
            ! rhel is the heliocentric distance of the asteroid
            ! corresponding to the moments the map is computed for
            ! lats and lons are the latitudes and longitudes to which
            ! the columns and the rows of the matrix correspond
            ! they are not exactly uniform
            ! the subroutine will return the lats and lons in radians
            ! and rhel in AU
            subroutine read_first_ratemap(fname, rhel)
                character(*), intent(in) :: fname
                real(8), intent(out) :: rhel
                integer i
                character(len = 18) strtmp
                character(len = 14) strtmp2
                real(8) x, y, z
                
                open(100, file = fname, status = 'old')
                    do i = 1, 5
                        read(100,*) strtmp
                    enddo
                    read(100,*) strtmp, x, y, z
                    rhel = sqrt(sum((/x, y, z/)**2)) * 1d3        ! in meters
                    rhel = rhel / AU                        ! in AU
                    read(100,*) strtmp
                    read(100,*) strtmp2, lons
                    read(100,*) strtmp2, lats
                    do i = 1, 3
                        read(100,*) strtmp
                    enddo
                    do i = 1, nlats
                        read(100,*) rmap1(:,i)
                    enddo
                close(100)
                lats = lats * deg2rad
                lons = lons * deg2rad
                
                lonmin = minval(lons)
                lonmax = maxval(lons)
            
            end subroutine read_first_ratemap
		
		
			! The axisymmetric distribution of ejection direction
			! distribution_shape is the parameter used to select
			! the expression for the PDF
			! wpsi is the polar angle in the coordinate system where
			! the distribution is axisymmetrical
			! psi is the polar angle in the horizontal coordinate system
			! lambdaM is the azimuth in the horizontal CS
			! zeta and eta are respectively zenith angle and azimuth
			! of the distribution symmetry axis in the horizontal CS
			function ejection_direction_distribution(distribution_shape, &
			                    wpsi, psi, lambdaM, zeta, eta) result(fpsi)
				use const
				use help
				implicit none
				integer N, i
				real(8), parameter :: normconst1 = 0.7723999d0
				real(8), parameter :: omega60 = 1.047198d0
				integer ii1, ii2, ind1, ind2
				real lat, lon, klat, klon, prerate1, prerate2
				real(8) utmp(3), uvec(3), xvec(3), yvec(3)
				real(8) xvec1(3), yvec1(3), zvec1(3)
				real(8) :: zvec(3) = (/0d0, 0d0, 1d0/)
				real(8) fpsi, Jpsi
				integer, intent(in) :: distribution_shape
				real(8), intent(in) :: psi, wpsi, lambdaM, zeta, eta
				
				select case(distribution_shape)
					case(0)
						fpsi = 0.25d0 / pi
					case(1)
					! pseudo Gaussian distribution of polar angle, uniform distribution of azimuth
						if(wpsi < pi) then
							fpsi = Exp(-(psi)**2 / 2d0 / omega60 / omega60)
							! this factor is normalization due to the fact that fpsi
							! domain is from 0 to pi and not from -infinity to +infinity
							fpsi = fpsi / normconst1
							fpsi = fpsi / twopi
						else
							fpsi = 0d0
						endif
					case(2)
					! Uniform distribution of polar angle inside a cone,
					! uniform distribution of azimuth
						if(wpsi <= omega60) then
							fpsi = 1d0 / (1d0 - cos(omega60)) / twopi
						else
							fpsi = 0d0
						endif
					case(3)
                    ! orts of the CS where we can easily know coords
                    ! of the ejection velocity vector
                        ! z-axis is along heliocentric radius
                        zvec1 = rMtmp / norma3d(rMtmp)
                        ! x-axis is towards local north
                        xvec1 = zvec - zvec1 * dot_product(zvec, zvec1)
                        xvec1 = xvec1 / norma3d(xvec1)
                        yvec1 = vector_product(zvec1, xvec1)
                        ! ejection velocity vector in this CS
                        utmp(1) = sin(psi) * cos(lambdaM) ! lambdaM is azimuth
                        utmp(2) = sin(psi) * sin(-lambdaM) ! counted from the local north clockwise
                        utmp(3) = cos(psi) 
                        ! ejection velocity vector in the same CS as rMtmp
                        uvec = utmp(1) * xvec1 + utmp(2) * yvec1 &
                             + utmp(3) * zvec1
                        uvec = uvec / norma3d(uvec)
                    ! the CS where the maps are defined: z-axis = (0, 0, 1)
                        ! projection of rMtmp to xy-plane
                        xvec = rMtmp - dot_product(zvec, rMtmp) * zvec
                        xvec = -xvec / norma3d(xvec)
                        yvec = vector_product(zvec, xvec)
                        ! ejection velocity vector in this CS
                        utmp(1) = dot_product(uvec, xvec)
                        utmp(2) = dot_product(uvec, yvec)
                        utmp(3) = uvec(3)
                        
                        lat = acos(utmp(3))            ! polar angle
                        lat = halfpi - lat                ! latitude
                        lon = atan(utmp(2), utmp(1))
                        if(lon < 0d0) lon = twopi + lon
                        lon = twopi - lon
                        
                        if(lat > lats(1) .and. lat < lats(nlats)) then
                            ii1 = 1
                            do while(lat < lats(ii1) .or. lat >= lats(ii1+1))
                                ii1 = ii1 + 1
                            enddo
                            ind1 = ii1
                            ii1 = ii1 + 1
                            klat = (lat - lats(ind1)) / (lats(ii1) - lats(ind1))
                        else
                        ! if the given latitude is beyond the possible
                        ! interpolation limits, we use the marginal values
                            if(lat < lats(1)) then
                                ii1 = 1 ; ind1 = 1
                                klat = 0d0
                            else
                                ii1 = nlats ; ind1 = nlats
                                klat = 0d0
                            endif
                        endif
                        if(lon < lonmax .and. lon > lonmin) then
                            ii2 = 1
                            do while(lon < lons(ii2) .or. lon >= lons(ii2+1))
                                ii2 = ii2 + 1
                            enddo
                            ind2 = ii2
                            ii2 = ii2 + 1
                            klon = (lon - lons(ind2)) / (lons(ii2) - lons(ind2))
                        else
                        ! if lon is beyond the values available for interpolation
                        ! we round the circle and use the values corresponding 
                        ! to the minimal and maximum longitudes
                            ii2 = 1
                            ind2 = nlons
                            if(lon >= lonmax) then
                            ! the number is 2 degrees between lon(180) = 359 deg and lon(1) = 1deg
                                klon = (lon - lonmax) / 0.03490659d0
                            else
                            ! 1 deg / 2 deg + lon / 2 deg
                                klon = 0.5d0 + lon / 0.03490659d0
                            endif
                        endif
                        ! interpolation over latitudes
                        prerate1 = ratemap(ind2, ii1) + klat * (ratemap(ind2, ind1) - ratemap(ind2, ii1))
                        prerate2 = ratemap(ii2, ii1) + klat * (ratemap(ii2, ind1) - ratemap(ii2, ii1))
                        ! interpolation over longitudess
                        fpsi = prerate1 + klon * (prerate2 - prerate1)
					case(4)
                        ! Here is the place to write the PDF 
                        ! of the distribution tailored for your needs
                        fpsi = 0d0

				endselect
				Jpsi = 1d0
				if(zeta /= 0d0) then
					Jpsi = Jacobian_tilt(psi, lambdaM, zeta, eta)
					fpsi = fpsi * abs(Jpsi)
				endif
				fpsi = fpsi * sin(wpsi)
					
			end function ejection_direction_distribution

			
						
			! Jacobian of coordinate transformation from vertical CS
			! to the CS with z-axis coinciding with the jet axis of symmetry
			! zeta and A are respectively zenith angle ang azimuth
			! of the distribution symmetry axis in the horizontal CS
			! psi and lambdaM are respectively polar angle
			! and azimuth of ejection in the horizontal CS
			function Jacobian_tilt(psi, lambdaM, zeta, A) result(J)
				use const
				implicit none
				real(8), intent(in) :: psi, lambdaM, zeta, A
				real(8) sinpsi, cospsi, sinzeta, coszeta
				real(8) J
				
				if(psi < 0.12d0 .and. zeta < 0.12d0) then
					J = psi / sqrt(psi*psi + zeta*zeta &
					               - 2d0 * psi * zeta * cos(lambdaM - A))
				else
					sinpsi = sin(psi) ; cospsi = cos(psi)
					sinzeta = sin(zeta) ; coszeta = cos(zeta)

					J = 4d0 * sinpsi / sqrt(10d0 - 2d0 * (cospsi - sinpsi) * (cospsi + sinpsi) &
						- 3d0 * cos(2d0 * (psi - zeta)) - 2d0 * (coszeta - sinzeta) * (coszeta + sinzeta) &
						- 3d0 * cos(2d0 * (psi + zeta)) &
						- 8d0 * cos(2d0 * (lambdaM - A)) * sinpsi * sinpsi * sinzeta * sinzeta &
						- 32d0 * cos(lambdaM - A) * sinpsi * cospsi * sinzeta * coszeta)
				endif
				
			end function Jacobian_tilt
			
			
			! This function represents the ejection speed distribution
			! (possibly time-dependent)
			! ud is  the parameter used to select the expression for the distribution
			! u is the ejection speed
			! R is the particle size
			function ejection_speed_distribution(ud, u) result(fu)
				use const
				use define_types
				implicit none
                real(8), parameter :: Rad = 1737d0 ! km
                real(8), parameter :: lambda = 200d0 ! km
                real(8), parameter :: hrel = Rad / lambda
                real(8), parameter :: lambda1 = 1d0 ! km
                real(8), parameter :: hrel1 = Rad / lambda1
                real(8), parameter :: uesc = 2.4d3 ! m/s
				real(8) Rc, nu
				real(8) urel, Rrel, fu
				type(ejection_speed_properties) ud
				real(8) u, q
								
				select case(ud%ud_shape)
					case(0)
						fu = 1d0 / (ud%umax - ud%umin) 
					
                    case(1)
                    ! from Szalay & Horányi, 2016, Lunar meteoritic gardening rate
                        nu = (u * AUdays2SI) / uesc
                        fu = 2d0 * hrel * nu / uesc / (1d0 - nu**2)**2 &
                             * exp(- hrel / (nu**(-2) - 1d0))
                        fu = fu * AUdays2SI
					case(2)
                        ! Custom distribution specification
                        fu = 0d0
				endselect
			
			end function ejection_speed_distribution

			

end module distributions_fun
