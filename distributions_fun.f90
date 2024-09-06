! This file is a part of DUDI-heliocentric, the Fortran-95 implementation 
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
	implicit none
		contains
		
		
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
				real(8) Rc, nu
				real(8) urel, Rrel, fu
				type(ejection_speed_properties) ud
				real(8) u, q
								
				select case(ud%ud_shape)
					case(0)
						fu = 1d0 / (ud%umax - ud%umin) 
					
					case(1)
                        ! Custom distribution specification
                        fu = 0d0
				endselect
			
			end function ejection_speed_distribution

			

end module distributions_fun
