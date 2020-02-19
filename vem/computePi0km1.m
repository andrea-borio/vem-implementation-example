function Pi0km1 = computePi0km1(VemData, area, cellMatrixH, cellPiNabla, ...
    cellMatrixH_Cholkm1)
%COMPUTEPI0KM1 Pi^0_{k-1} projection matrix
%   PI0KM1 = COMPUTEPI0KM1(VEMDATA, AREA, CELLMATRIXH, CELLPINABLA,
%   CELLMATRIXH_CHOLKM1) returns the Pi^0_{k-1} projection matrix for VEM
%   whose parameters are in VEMDATA. AREA is the area of the polygon,
%   CELLMATRIXH is the local mass matrix of scaled monomials of order k,
%   CELLPINABLA is the Pi^{nabla}_k projection matrix of the polygon,
%   CELLMATRIXH_CHOLKM1 is the Cholesky decomposition of the mass matrix of
%   scaled monomials of order k-1
%
%   See also INITVEMDATA, COMPUTEPINABLA

N_boundaryDofs = size(cellPiNabla,2)-VemData.N_km2;

cellMatrixC = [ ...
    [zeros(VemData.N_km2,N_boundaryDofs), area*eye(VemData.N_km2,VemData.N_km2)];
    cellMatrixH((VemData.N_km2+1:VemData.N_km1),:)*cellPiNabla ];

Pi0km1 = cellMatrixH_Cholkm1\(cellMatrixH_Cholkm1'\cellMatrixC);

end
