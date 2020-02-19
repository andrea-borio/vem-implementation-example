function Pi0km1Grad = computePi0km1Grad(VemData, VanderBoundary, ...
    BoundaryWeightsNormal, area, cellDerivationMatrices, cellMatrixH_Cholkm1)
%COMPUTEPI0KM1GRAD Pi^0_{k-1}nabla projection matrices on a polygon.
%   PI0KM1GRAD = COMPUTEPI0KM1GRAD(VEMDATA, VANDERBOUNDARY,
%   BOUNDARYWEIGHTSNORMAL, AREA, CELLDERIVATIONMATRICES,
%   CELLMATRIXH_CHOLKM1) returns a 2-by-1 cell containing the matrix
%   representation of the k-1 projection of x and y derivatives for VEM
%   whose parameters are in VEMDATA. VANDERBOUNDARY is the Vandermonde
%   matrix of monomials of order k on the boundary, BOUNDARYWEIGHTSNORMAL
%   is the 2-by-1 cell containing the vectors of boundary quadrature
%   weights scaled by the components of the outward normal, AREA is the
%   area of the polygon, CELLDERIVATIONMATRICES is the 2-by-1 cell
%   containing the derivation matrices for monomials of order k on the
%   polygon, CELLMATRIXH_CHOLKM1 is the Cholesky decomposition of the mass
%   matrix of scaled monomials of order k-1 on the polygon.
%
%   See also INITVEMDATA

Pi0km1Grad = cell(2,1);
for d = 1:2
    cellMatrixE = [ ...
        VanderBoundary(:,1:VemData.N_km1)'*diag(BoundaryWeightsNormal{d}), ...
        -area*cellDerivationMatrices{d}(1:VemData.N_km2,1:VemData.N_km1)' ];
    Pi0km1Grad{d} = cellMatrixH_Cholkm1\(cellMatrixH_Cholkm1'\cellMatrixE);
end

end

