# Bezier extraction of TSplines

These few code files demonstrate the Bezier extraction for TSplines based 
on the work "Isogeometric finite element data structures based on Bezier 
extraction of T-Splines", by M. A. Scott, M. J. Borden, C. V.
Verhoosel, T. W. Sederberg,  and T. J. R. Hughes, in International
Journal for Numerical Methods in Engineering 2011, DOI: 10.1002/nme.3167

The algorithm from this paper had some minor errors which are corrected
in this application. 

To get a demo of how Bezier extraction for TSplines in 1D (aka B-splines)
simply call visualize_bezier_extraction(xi, knot_in) with a local knot 
vector xi of a BSpline together with an additional argument knot_in to 
specify whether to insert additional knots inside the knot vector xi 
for the bezier extraction. Specify knot_in = [] for no inputs. 

There are some examples given in visualize_bezier_extraction, however,
we encourage you to independently try this on your own. 

The results of the referenced paper are verified to rounding errors at the 
end of visualize_bezier_extraction. 
