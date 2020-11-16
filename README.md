# Correlation_integral_SoS
A programme to compute the dimension of the trace of an orbit on the Surface of Section

Surface of section is useful for finding the presence of a third integral of motion in an axis-symmetric potential, so this program analizes a .txt file where the first 7 colomns are supposed to be

$ t   R   z   phi   vR    vz    vphi$

and compute the dimention of the trace of the orbit.

The .txt file you want to analize should be in the same folder where you're launching the programme

At the moment the programs works just with a .txt file already containing just surface of section data, maybe future implementation will make the programme computer the SoS on its own.

This little programme has been written using Python 3.8.3, but probably works even with earlier version of Python, you just have to try :)
