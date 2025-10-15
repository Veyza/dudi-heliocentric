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

! File: define_types.f90
! Description: Contains definitions for structures used throughout
! the routines in the package.

! Author: Anastasiia Ershova
! E-mail: vveyzaa@gmail.com

module define_types
  use const
  implicit none
  
  type ejection_speed_properties
    integer ud_shape                ! parameter defining which distribution is used for ejection speed
    real(8) umax                    ! gas velocity
    real(8) umin                    ! a parameter in velocity distribution
  end type ejection_speed_properties
  
  type source_properties           ! parameters of dust ejection
    real(8) rrM(3)               ! Cartesian coordinates of a point
                                 ! source in the heliocentric coordinate 
                                 ! system
    real(8) r                    ! heliocentric distance of the point source
    real(8) alphaM               ! polar angle of the point source
    real(8) betaM                ! eastern longitude of the point source
    real(8) zeta                 ! zenith angle of the axis around which 
                                 ! ejection is symmetrical
    real(8) eta                  ! azimuth of this axis (counted from 
                                 ! the local North, clockwise)
    real(8) symmetry_axis(3)             ! unit vector in heliocentric 
                                         ! coordinate system pointing to the
                                         ! direction of the axis around
                                         ! which ejection is symmetrical
    type(ejection_speed_properties) ud   ! parameters of ejection speed
                                         ! distribution
    integer ejection_angle_distr         ! parameter defining which ejection
                                         ! angle distribution is used
    real(8) Nparticles                   ! number of particles ejected by
                                         ! the source in total
    real(8) Tj                    ! moment of dust ejection (for simple 
                                  ! expansion and delta-ejection methods, 
                                  ! central moment of the activity interval 
                                  ! for v-integration method
    real(8) dtau                  ! Length of the time interval centered at
                                      ! Tj over which the given source ejected
                                      ! dust
  end type source_properties

  type ephemeris
    real(8) coords(3)                ! heliocentric coordinates of the comet
    real(8) Vastvec(3)               ! heliocentric velocity of the comet
    real(8) Vast                     ! comet's heliocentric speed
  end type ephemeris
  
  type position_in_space       ! where the dust density is to be calculated
    real(8) r                        ! distance from the Sun's center
    real(8) alpha                    ! polar angle
    real(8) beta                     ! eastern longitude
    real(8) rvector(3)               ! Cartesian coordinates of the point 
  end type position_in_space


end module define_types
