function [PiNabla,cellMatrixG,cellMatrixB] = computePiNabla( ...
    VemData, area, laplacianMatrix, ...
    VanderInternalDerivatives, InternalWeights, ...
    VanderBoundaryDerivatives, BoundaryWeightsNormal, ...
    VanderBoundary, BoundaryWeights, VanderInternal )
%COMPUTEPINABLA Pi^nabla_k projection matrix on a polygon.
%   PINABLA = COMPUTEPINABLA(VEMDATA, AREA, LAPLACIANMATRIX,
%   VANDERINTERNALDERIVATIVES, INTERNALWEIGHTS, VANDERBOUNDARYDERIVATIVES,
%   BOUNDARYWEIGHTSNORMAL, VANDERBOUNDARY, BOUNDARYWEIGHTS, VANDERINTERNAL)
%   returns the Pi^nabla_{k} projection matrix for VEM whose parameters are
%   in VEMDATA. AREA is the area of the polygon, LAPLACIAMATRIX matrix
%   representing the laplacian operator for monomials of order k on the
%   polygon, VANDERINTERNALDERIVATIVES is the Vandermonde matrix of
%   derivatives of monomials on internal quadrature points, INTERNALWEIGHTS
%   is the vector of internal quadrature points, VANDERBOUNDARYDERIVATIVES
%   is the Vandermonde matrix of derivatives of monomials on boundary
%   quadrature points, BOUNDARYWEIGHTSNORMAL is the 2-by-1 cell containing
%   the vectors of boundary quadrature weights scaled by the components of
%   the outward normal, VANDERBOUNDARY is the Vandermonde matrix of scaled
%   monomials on the boundary, BOUNDARYWEIGHTS is the vector of boundary
%   quadrature weights, VANDERINTERNAL is the Vandermonde matrix of scaled
%   monomials on internal quadrature points.
%
%   [PINABLA,CELLMATRIXG,CELLMATRIXB] = COMPUTEPINABLA(...) also returns
%   the auxiliary matrices G and B.
%
%   See also INITVEMDATA

cellMatrixG = ...
    VanderInternalDerivatives{1}'*diag(InternalWeights)*VanderInternalDerivatives{1} + ...
    VanderInternalDerivatives{2}'*diag(InternalWeights)*VanderInternalDerivatives{2};
cellMatrixB = [ ...
    VanderBoundaryDerivatives{1}'*diag(BoundaryWeightsNormal{1}) + ...
    VanderBoundaryDerivatives{2}'*diag(BoundaryWeightsNormal{2}) , ...
    -area*laplacianMatrix(1:VemData.N_km2,:)'];

if VemData.order == 1
    cellMatrixG(1,:) = BoundaryWeights'*VanderBoundary;
    cellMatrixB(1,:) = BoundaryWeights';
else
    cellMatrixG(1,:) = InternalWeights'*VanderInternal;
    cellMatrixB(1,size(VanderBoundary,1)+1) = area;
end

PiNabla = cellMatrixG\cellMatrixB;


end
