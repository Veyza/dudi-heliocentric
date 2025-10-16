! This file is a part of DUDI-heliocentric, the Fortran-90 implementation 
! of the two-body model for the dynamics of dust ejected from an atmosphereless
! body moving around the Sun
! Version 1.0.2
! This is free software. You can use and redistribute it 
! under the terms of the GNU General Public License (http://www.gnu.org/licenses/)
! If you do, please cite the following paper
! Anastasiia Ershova and Jürgen Schmidt, 
! Two-body model for the spatial distribution of dust ejected from
! an atmosphereless body, 2021, A&A, 650, A186 

! File: const.f90
! Description: Fundamental constants and numerical parameters utilized across
! various subroutines in the package.

! Author: Anastasiia Ershova
! E-mail: vveyzaa@gmail.com

module const
    implicit none
        
        ! fundamental constants and auxilary numbers
        real(8), parameter :: pi = 3.141592653589793d0
        real(8), parameter :: halfpi = pi / 2d0
        real(8), parameter :: sqrtpi = sqrt(pi)
        real(8), parameter :: twopi = 2d0 * pi
        real(8), parameter :: sqrt2d0 = sqrt(2d0)
        real(8), parameter :: deg2rad = pi / 1.8d+2
        real(8), parameter :: rad2deg = 1.8d+2 / pi
        real(8), parameter :: gravity_constant = 6.674d-11      ! m^3 / kg / s^2
        real(8), parameter :: speed_of_light = 2.99792458d+8   ! m/s
        
        ! parameters defining the motion around the Sun
        real(8), parameter :: Msun = 1.9891d+30        ! kg
        real(8), parameter :: AU = 1.495978707d+11     ! Astronomical Unit in SI
        real(8), parameter :: GMsun = 1.327124400419393e+20 / AU**3 * 864d2**2
        real(8), parameter :: Lsun = 3.828d+26         ! Solar luminosity in SI    
        real(8), parameter :: s_in_day = 86400d0
        real(8), parameter :: AUdays2SI = AU / s_in_day

        real(8), parameter :: rho = 2.5d3        ! kg/m^3 

        ! number of steps in Gauss-Legendre integration
        ! for the v-integration method
        ! must be 5, 10, 20, 30, 40, 50, or 64
        integer, parameter :: order_v = 20
        
        real(8), parameter :: check_inf = 9d+38  
        integer :: N_of_warnings = 0        ! warnings counter
    integer, parameter :: maxNofWarnings = 100  ! maximum number of
                                                ! warnings printed 
                                                ! out to the text file
                                                ! fort.666
            
end module const
