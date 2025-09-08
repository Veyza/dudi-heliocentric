= README for DUDI-heliocentric

This repository contains DUDI-heliocentric, a Fortran-95 code that models dust
ejection from an atmosphereless body orbiting the Sun. This software builds upon
the research described by 
Ershova, A. & Schmidt, J. (2021). Two-body model for the spatial distribution 
of dust ejected from anatmosphereless body, 2021, A&A, 650, A186, 
and extends the initial DUDI model applicable to moons of giant planets.
Find the original DUDI at https://github.com/Veyza/dudi.

Derivations of the formulae for DUDI-heliocentric are detailed in a forthcoming
publication 
Ershova, A., Schmidt, J., Liu, X., Szalay, J., Kimura, H., Hirai, T., Arai,
T., and Kobayashi, M., A computationally efficient semi-analytical model for the dust
environment of comets and asteroids, A&A 693, A80 (2025).

== Prerequisites 

To compile and run DUDI-heliocentric, the following is required: 
* gfortran - Fortran compiler 

For generating plots from the example applications described in this README, 
the following additional software is needed: 
* Python3, NumPy, and Matplotlib 
Note: If these are not installed, DUDI-heliocentric will still produce results, 
but the plot will not be generated. 

== License

DUDI-heliocentric is distributed under GNU GENERAL PUBLIC LICENSE Version 3.

== Table of Contents

1. Specifying Functions Describing the Dust Ejection
2. Defining Pseudo Gravitational Parameter `muR`
3. Supplying Input Data
4. Choosing Calculation Method
5. Compilation and Running an Example
6. Controlling the Accuracy of Calculations
7. Modeling the dust environment of Phaethon
8. Possible Issues
9. Latest updates

== 1. Specifying Functions Describing the Dust Ejection

Customize dust ejection characteristics in distributions_fun.f90. This file 
includes functions that define the speed and direction distribution of dust 
ejection. Each function utilizes Fortran's `select case` operator, where 
a selector specifies a different distribution to compute. Several distributions 
are already implemented and are used in the example applications.

*`ejection_speed_distribution(ud, u)` calculates the probability
density function (PDF) of ejection speeds based on predefined parameters:

- `ud`(type(ejection_speed_properties)): Structure containing parameters of 
the ejection speed distribution, defined as follows: 
  - `ud%ud_shape`(integer): Selector for the ejection speed distribution.
  - `ud%umin`(real(8)): Minimum ejection speed (m/s).
  - `ud%umax`(real(8)): Maximum ejection speed (m/s).
- `u`(real(8)): Ejection speed (m/s)

The output of this function is stored in the variable `fu`. To customize the
distribution to fit specific needs, locate the `case(1)` block in the function
implementation:

case(1)
! Custom distribution specification
fu = 0d0

Replace `fu = 0d0` with the desired expression for the PDF.

*`ejection_direction_distribution` defines the angular distribution
of ejected particles:
- `distribution_shape` (integer): Selector for the size distribution type.
- `wpsi`(real(8)): Polar angle in the coordinate system in which dust ejection 
is axisymmetric (radians).
- `psi`(real(8)):  Polar angle in the local horizontal coordinate system (radians).
- `lambdaM`(real(8)):  Azimuth in the local horizontal coordinate system (radians)
counted clockwise from the local north.
- `zeta`(real(8)):  Zenith angle of the distribution symmetry axis in 
the local horizontal coordinate system (radians).
- `eta`(real(8)):  Azimuth of the distribution symmetry axis in the local 
horizontal coordinate system (radians).

The value of this function is stored in the variable `fpsi`. The PDF is primarily
defined for the variable `wpsi` in a coordinate system where the ejection is
axisymmetric, assuming azimuth of ejection is distributed uniformly (1/2π).

To implement a custom distribution, modify the function within the `case(3)` block:

case(3)
! Custom PDF specification
fpsi = 0d0

Replace `fpsi = 0d0` with the appropriate PDF expression for your simulation.


== 2. Defining Pseudo Gravitational Parameter `muR`

The pseudo gravitational parameter, `muR`, is essential for simulating the
dynamics of dust particles under solar gravity and radiation pressure. It is
defined as:

----
muR = GMsun * (1 - beta(R, Qpr))
----

Where:
- `GMsun` is the solar gravitational parameter (AU^3/day^2).
- `R` is the radius of the dust particle.
- `beta(R)` is the ratio of radiation pressure force to gravitational force.

### Calculation of muR

The subroutine `reduced_gravitational_parameter(Rg, Qpr, muR)` in `data_in.f90`
calculates `muR` using the expression for `beta` from Burns et al., 1979
(Equation 19).

#### Subroutine Parameters:
- `Rg` (real(8), Input): Radius of the dust particle (meters).
- `Qpr` (real(8), Input): Radiation pressure efficiency (dimensionless).
- `muR` (real(8), Output): Computed pseudo gravitational parameter (AU^3/day^2).


== 3. Supplying Input Data

All the calculations are performed by the subroutines 
`hc_DUDI_v_integration(density, point, source, &
                        muR, tnow, comet, Rast_AU, pericenter)`
                     or
`hc_DUDI_delta_ejection(density, point, source, muR, dt, comet, Rast_AU)`
                     or
`hc_DUDI_simple_expansion(density, source, dt, cloudcentr, point)`
                     
which will return the result in the variable density.

- `density` (real(8)): Stores the computed number density of dust particles.

The input variables of the subroutines: 

- `tnow` (real(8)): Specifies the time moment (days) at which number density 
                    is calculated. (Input of the v-integration method)
- `dt` (real(8)): Specifies the time (days) elapsed since the moment of dust 
                  ejection. (Input of the delta-ejection and simple expansion 
                  methods)
- `muR` (real(8)): The pseudo gravitational parameter adjusted 
                   for radiation pressure.
- `Rast_AU` (real(8)): Radius of the dust-emitting body (AU),
                       assuming a spherical shape. This can be set to zero if 
                       body size and dust re-collisions are neglected.
- `pericenter` (logical): Indicates if dust particles have passed perihelion in
                          their orbit from the source to the point of interest, 
                          relevant only for `hc_DUDI_v_integration`.
- `cloudcentr` (real(8)) : A 3d-vector of heliocentric position of the prime 
                           cloud center (AU), relevant only for 
                           `hc_DUDI_simple_expansion`

* `source` structure details the properties of the dust source:

- `source%r` (real(8)):      Heliocentric radial distance of the source (AU).
- `source%alphaM` (real(8)): Polar angle of the source (radians)
- `source%betaM` (real(8)):  Eastern longitude of the source
- `source%rrM` (real(8)):    3d Cartesian coordinates of the source:

`source%rrM(1) = source%r * sin(source%alphaM) * cos(source%betaM)`
`source%rrM(2) = source%r * sin(source%alphaM) * sin(source%betaM)`
`source%rrM(3) = source%r * cos(source%alphaM)`

- `source%zeta` (real(8)): Zenith angle of the ejection symmetry axis counted 
                           from the anti-solar direction (radians).
-`source%eta` (real(8)): Azimuth of the ejection symmetry axis counted clockwise
                         from the local north (radians).
- `source%symmetry_axis` (real(8) 3d-vector): Unit vector pointing along the
                                              ejection symmetry axis, derived 
                                              from the orientation parameters.
- `source%ejection_angle_distr` (integer): Selector for the ejection direction 
                                           distribution
- `source%ud` (type(ejection_speed_properties)): Structure containing parameters
                                  of the ejection speed distribution, including:
   -- `ud_shape` (integer): Selector for the ejection speed distribution PDF.
   -- `umin` (real(8)): Minimum possible ejection speed (AU/day).
   -- `umax` (real(8)): Maximum possible ejection speed (AU/day).

*`point` structure specifies the point in space where the dust number density
is to be calculated (referred in the paper as the "point of interest" or 
"spacecraft position")

- `point%r` (real(8)): Heliocentric radial distance of the point (AU)
- `point%alpha` (real(8)): Polar angle of the point (radians)
- `point%beta` (real(8)): Estern longitude of the point (radians).
- `point%rvector` (real(8)): 3d-vector of the point's Cartesian coordinates:

`point%rvector(1) = point%r * sin(point%alpha) * cos(point%beta)
point%rvector(2) = point%r * sin(point%alpha) * sin(point%beta)
point%rvector(3) = point%r * cos(point%alpha)`

`comet` structure describes the dust-emitting body at the time of ejection:

- `coords`(real(8)): Heliocentric coordinates of the asteroid (AU)
- `Vastvec` (real(8)): Heliocentric velocity vector of the asteroid (AU/day)
- `Vast` (real(8)): Magnitude of the body's velocity (AU/day).


== 4. Choosing Calculation Method

DUDI-heliocentric offers three algorithms to calculate the dust number density,
each with unique applicability limits and CPU time requirements. Detailed
mathematical descriptions and applicability analyses are available in the
associated paper.

=== Overview of Methods

- **Simple Expansion Method:**
  This is the fastest algorithm and is capable of modeling dust dynamics under
  arbitrary forces. However, its accuracy decreases fastest with time elapsed
  since the dust ejection.

- **Delta-Ejection Method:**
  Focuses on dust dynamics influenced solely by solar gravity and radiation
  pressure. It requires less CPU time compared to v-integration but has stricter
  applicability limits regarding time since dust ejection and orbital parameters
  of the dust-emitting body.

- **V-Integration Method:**
  Also handles dust dynamics under solar solar gravity and radiation
  pressure, but with broader applicability in terms of time and orbital 
  parameters than the delta-ejection method.

=== Method Selection Routine

A specific routine in this GitHub repository helps determine the most suitable
method for given scenarios. This routine utilizes parameters defined in
`orbit_and_time_test.dat` located in the `input_data_files` directory:

- **Line 1**: Cartesian coordinates of the comet at the moment of dust ejection 
             (AU).
- **Line 2**: Velocity vector of the comet at the moment of dust ejection 
             (AU/day).
- **Line 3**: Time interval from ejection to density computation (∆t, days).
- **Line 4**: Interval of the source activity for v-integration method
              (∆τ, days).

Executing `make select` compiles and runs this routine, outputting recommendations
to the terminal based on the input parameters. The routine calculates the dust 
number density by the three methods and compares the results. A non-uniform 
distribution of the ejection direction is employed to ensure a proper test for 
the simple expansion method.

=== Results

The routine `select_method_dudihc` also outputes the obtained results to 
the following  files in  the folder "./results":
* "test_simple_exp_meth.dat": the number density obtained by the simple expansion 
method ("S-matrix" 200x200)
* "test_delta-eject_meth.dat": the number density obtained by the delta-ejection
method ("D-matrix" 200x200)
* "test_v-integr_meth.dat": the number density obtained by the v-integration
method ("V-matrix" 200x200)
* "test_simp_exp_vs_delta-eject.dat": discrepancies between the result of the 
simple expansion and delta-ejection methods computed as (S/D - 1)
* "test_delta-eject_vs_v-integr.dat": discrepancies between the result of the 
delta-ejection and v-integration methods computed as (D/V - 1)

=== Practical Recommendations

The routine evaluates the spatial distribution of dust in a single prime cloud 
using all three methods. It suggests the simple expansion or delta-ejection 
methods if their results closely match those from the v-integration method.
The desired accuracy can be set directly in the file `select_method.f90` by
defining the real paremeter `accuracy`. The accuracy is set in per cent.

For specific orbital shapes and phases:
- A small ∆t may favor delta-ejection.
- A suitable ∆τ interval, is typically between 0.001 and 0.2 seconds, though 
adjustments may be necessary for eccentric orbits or ejections near perihelion.

It's recommended to start with ∆τ = 0.05 and ∆t = 0.1 and to run the command
`make select` adjusting these values to find the value of ∆τ and ∆t for
which the delta-ejection method is recommended. In this way one secures the 
optimal ∆τ for future simulations. Then start increasing ∆t to identify the 
limits of delta-ejection and simple expansion methods applicability.


== 5. Compilation and Running an Example

To compile DUDI-heliocentric and create the executable:

----
make example
----

This command compiles the source files using the makefile provided in the
repository, resulting in an executable named `dudihc`.

=== Running the Supplied Example

== Example Routine Overview in DUDI-heliocentric

The "examples/example.f90" program demonstrates the package's basic functionality and is
set up as follows:

- **Asteroid Model:**
  - Radius: 3d3 meters (`Rast`), also expressed in astronomical units (`Rast_AU`).
  - Trajectory and velocity: Loaded from "ephemeridae.dat" in "input_data_files".
    Each line records time (days), 3D position (AU), and velocity (AU/day).

- **Dust Ejection Simulation:**
  - Number of positions (`Np`) and sources per position (`Ns`) are configurable.
  - Sources eject dust uniformly into 60-degree cones, with ejection speeds 
    uniformly distributed between 5 and 100 m/s.
  - The number of particles ejected by each source is defined as
    `N = 1d5 + 1d5 * cos(a)`, where `a` is the angle between asteroid velocity 
    source position normal, enhancing ejection toward the asteroid's apex.

- **Gravitational Parameter (`muR`):**
  - Calculated for a dust grain radius (Rg) and radiation pressure efficiency (Qpr),
    and set in such a way that beta = 0.4, resulting in muR = 0.6 * GMsun.

- **Density Calculation:**
  - Conducted at the grid of points in the asteroid orbital frame centered on 
    the asteroid last position loaded from "ephemeridae.dat".
  - Number density at grid points is the cumulative dust from (Np-1) * Ns
    sources, computed using the delta-ejection method. Sources at the final
    position (Np) eject dust at the moment for which the number density is 
    calculated, thus their contribution is excluded.

This configuration ensures efficient processing, with a typical run time of
about 40 seconds on a 4-cored PC for `Np = 41` and `Ns = 50`.
`Np = 41` corresponds to dust ejection for 20 min, as the time intervals between
the ephemeridae in the file "ephemeridae.dat" is 30 sec.


=== Running the Supplied Example with Visualization

To visualize the results of the dust ejection simulation:

----
make example_image
----

Executing this command will:

- Compile the `bins/dudihc` executable, if it is not already compiled.
- Run the executable to simulate dust dynamics, generating the necessary
  output data.
- Utilize the Python script `show_image.py` to visualize the simulation results.

-- Numerical results are saved to `./results/result.dat`. 
   This will be a matrix 200x200.
-- A graphical representation is saved as `./results/result.png`.


== 6. Controlling the Accuracy of Calculations

The program monitors numerical accuracy. If a poorly conditioned case occurs
where accuracy falls below standard levels, it generates a file named `fort.666`.
This file logs details about the issues, aiding in pinpointing the problems.

Even with several indications of poor accuracy in `fort.666`, this does not
completely discredit the overall solution. Numerical difficulties may appear
at some integration steps or for some "source and spacecraft position" pairs, 
but not necessarily at all steps.

Furthermore, if the number of warnings in `fort.666` exceeds the limit set by
`maxNofWarnings` in the `const.f90` module, the program will stop and print a
warning message to the command line. 


== 7. Modeling dust environment of Phaethon

This repository includes the routines and input files to reproduce 
Figure 10 from Ershova et al. (2025), Astronomy & Astrophysics. 
The input data are the ephemeridae of the asteroid saved in the file 

"input_data_files/Phaethon_2025-02-22_last_int=10min_ECLIPJ2000.dat,"

a table of beta values corresponding to the given particle radii saved in 

"input_data_files/beta_vs_forsterite_Rg.dat," 

and four maps of impact ejecta in the directory 

"input_data_files/impact_ejecta_maps."

The command:

----
make phaethon
----

compiles the program, runs it, and, if Python 3 is installed, generates 
the plot. The results of the calculations and the plot are saved in the "results" 
directory. The results are text files containing matrices of number densities at 
points around Phaethon. These number densities are not yet weighted by the size 
distribution. Integration over particle radii is performed by the Python script.
The result of the integration is written to the file "integrated_Rmin=0.10.dat". 
On a 4-core PC, the program’s typical run time is under 15 minutes. 

To estimate dust number density at different points, assume different ejection 
parameters, or use a custom beta(R) relation, modify the subroutines in the 
"phaethon_input.f90" module.


== 8. Possible Issues

If you encounter a "Segmentation fault" when running a program compiled with 
the `-openmp` key (as specified in our Makefile), it may be due to the stack 
size limit.
To address this:

* In `bash`, increase the stack size by running:
  `ulimit -s unlimited`
* In `tcshell`, use:
  `limit stacksize unlimited`
  
  
== 9. Latest updates
    Version 1.0.1
  - Restructured repository with clear src/, examples/, scripts/, bin/, build/, and results/ folders.
  - Updated Makefile and .gitignore accordingly
