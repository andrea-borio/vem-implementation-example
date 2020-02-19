function [PiNabla_VemDofs, cellMatrixD] = computePiNablaVemDofs(VemData, ...
    VanderBoundary, area, cellMatrixH, PiNabla)
%COMPUTEPINABLAVEMDOFS Pi^{nabla,dof}_k projection matrix on a polygon.
%   PINABLA_VEMDOFS = COMPUTEPINABLAVEMDOFS(VEMDATA, VANDERBOUNDARY, AREA,
%   CELLMATRIXH, PINABLA) returns the matrix representation of the operator
%   Pi^nabla_k defined as an operator onto the VEM space defined by
%   VEMDATA. AREA is the area of the polygon, CELLMATRIXH is the mass
%   matrix of scaled monoimals on the polygon, PINABLA the matrix
%   representation of the same operator defined onto the space of
%   polynomials.
%
%   [PINABLA_VEMDOFS, CELLMATRIXD] = COMPUTEPINABLAVEMDOFS(...) also
%   returns the auxiliary matrix D.
%
%   See also INITVEMDATA, COMPUTEPINABLA

cellMatrixD = [ VanderBoundary ; (1/area)*cellMatrixH(1:VemData.N_km2,:)];

PiNabla_VemDofs = cellMatrixD*PiNabla;

end

