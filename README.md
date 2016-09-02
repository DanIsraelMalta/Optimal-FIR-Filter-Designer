# optimal filter design

This project is a part of my private assement of two frameworks:
* the Ultimate++ (c++) framework, assesed as a rapid GUI prototyper.
* Eigen linear algebra framework, assesed as a Matlab replacement.
 
Anyway, given cutoff frequency (bandwidth; 3db frequency), sampling frequency and filter order,
OptimalFilter will return an all-pole, near-linear phase low pass filter with optimized
magnitude response in the passband region.

remarks:
* since the filter is "near linear" in phase, it should be forward/backward filtered
 on implementation if non linear phase shift is not acceptable.
* left click with the mouse on the graphic output to open a label displaying 
 point information.
