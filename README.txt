MINORBIT2
========================
A simple minor planet orbit simulator with n-body gravitational perturbations from major planets.
Uses SPICE kernels for planets, JPL Horizons for initial conditions and validation, and an 8th order symplectic propagator for minor planet orbits.
Results can be visualized in an animated 3D orbit plot by using the visualizer utility.

You can run minorbit2.py with '--help' argument to print out the help text.

The repository includes de421.bsp as the only SPICE kernel due to size limitations, but you can also download and use other kernels such as DE440.
SPICE kernels are retrieved from NAIF. See REFERENCES.

PREPARING INPUT FILES
========================
example.txt is an example input file.

You should define the following inputs for each input file:

  T0 YYYY-MM-DD
  TF YYYY-MM-DD
  DT XX
  RF example_result.txt

These lines define the simulation start date (T0), simulation ending date (TF), fixed time-step size (DT) and output filename (RF).
After the headers, you can choose which minor planets to include in the simulation as follows:

  MP 2017 BX232
  MP 2017 AC64

LICENSE
========================
Minorbit2 is licensed under the terms of MIT License.
https://github.com/arda-guler/minorbit2/blob/master/LICENSE

CITATION
========================
If you use Minorbit2 or derived software for your research, I'd like to hear about it!

Example APA citation:

Guler, H. A., (2024). Minorbit2 (Version 2.2.0) [Source code].
Retrieved from https://github.com/arda-guler/minorbit2

Please also credit NAIF (https://naif.jpl.nasa.gov/naif/credit.html) for SPICE data, and NASA JPL SSD (https://ssd.jpl.nasa.gov/about/) for JPL Horizons System API data.

REFERENCES
========================
    Acton, C.H.; "Ancillary Data Services of NASA's Navigation and Ancillary Information Facility;" Planetary and Space Science, Vol. 44, No. 1, pp. 65-70, 1996.
    DOI 10.1016/0032-0633(95)00107-7

    Charles Acton, Nathaniel Bachman, Boris Semenov, Edward Wright; A look toward the future in the handling of space science mission geometry; Planetary and Space Science (2017);
    DOI 10.1016/j.pss.2017.02.013

