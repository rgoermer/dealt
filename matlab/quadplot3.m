function hh = quadplot3(quad,varargin)
%TRIPLOT Plots a 2D triangulation
%   QUADPLOT(QUAD,X,Y,Z) displays the quadrilaterals defined in the
%   M-by-8 matrix QUAD.  A row of QUAD contains indices into X,Y,Z that
%   define a single quadrilateal. The default line color is blue.
%
%   QUADPLOT3(...,COLOR) uses the string COLOR as the line color.
%
%   H = QUADPLOT3(...) returns a vector of handles to the displayed 
%   quadrilaterals
%
%   QUADPLOT3(...,'param','value','param','value'...) allows additional
%   line param/value pairs to be used when creating the plot.
%
%   See also TRISURF, TRIMESH, DELAUNAY, TriRep, DelaunayTri.
%
%   Script code based on copyrighted code from mathworks for TRIPLOT.
%   Allan P. Engsig-Karup, apek@imm.dtu.dk.

error(nargchk(1,inf,nargin,'struct'));
start = 1;
x = varargin{1};
y = varargin{2};
z = varargin{3};
quads = quad;
if (nargin == 4) || (mod(nargin-4,2) == 0)
    c = 'black';
    start = 4;
else
    c = varargin{4};
    start = 5;
end
  
d = quads(:,[1 2 4 3 1 5 6 8 7 5 1 3 7 8 4 2 6])';
h = plot3(x(d), y(d), z(d) ,c, varargin{start:end});
if nargout == 1, hh = h; end
