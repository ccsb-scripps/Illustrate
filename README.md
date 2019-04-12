--------------------------------------------------------------------
                 ******************************
                     I L L U S T R A T O R
                   Biomolecular Illustration
		  ******************************
		        David S Goodsell
--------------------------------------------------------------------

		no anti-aliasing of edges
	July 19 2007--conical shadows, smooth outlines with kernels,
                     fixed bond intersections, ppm output

--------------------------------------------------------------------
Compilation
>gfortran illustrator-2016.f -o illustrator-2016
>chmod +x  illustrator-2016
--------------------------------------------------------------------
Examples
>cd example
>../illustrator-2016 < example/2hhb.inp
