function [cellStiffness, cellRightHandSide ] = buildLocalSystem(mesh, ...
    cellId, VemData, refTriQuad, forcingTerm, diffusionParameter)
%BUILDLOCALSYSTEM matrix and right-hand side for the VEM discretization of
%the Laplace problem with homogeneous boundary conditions
%   [CELLSTIFFNESS, CELLRIGHTHANDSIDE] = buildLocalSystem(MESH, CELLID,
%   VEMDATA, REFTRIQUAD, FORCINGTERM, DIFFUSIONPARAMETER) returns the local
%   matrix and right-hand side on the polygon whose position in MESH is
%   CELLID. VEMDATA is the structure describing the cell-independent
%   properties of the method, REFTRIQUAD is the matrix containing
%   quadrature points and weights on the reference triangle,
%   FORCINGTERM(X,Y) and DIFFUSIONPARAMETER(X,Y) must return the values of
%   the forcing term and the diffusion parameter in a given set of points
%   (X,Y)
%
%   See also INITVEMDATA, TRIGAUSSPOINTS


order = VemData.order;
alphas = VemData.monomialExponents;
numVertices = length(mesh.polygon(cellId).vertices);
nDOF = order*numVertices + VemData.N_km2;

% Compute cell dependent quantities
x_vert = vertcat(mesh.vertex(mesh.polygon(cellId).vertices).x);
y_vert = vertcat(mesh.vertex(mesh.polygon(cellId).vertices).y);
diam = sqrt(max(diff(x_vert(nchoosek(1:numVertices,2)),1,2).^2 + ...
    diff(y_vert(nchoosek(1:numVertices,2)),1,2).^2));
area = polyarea(x_vert,y_vert);
centroid = computeCellCentroid(x_vert,y_vert,area);
cellDerivationMatrices = { (1/diam)*VemData.derivationMatrices{1} ; ...
    (1/diam)*VemData.derivationMatrices{2} };
cellLaplacianMatrix = (1/(diam*diam))*VemData.laplacianMatrix;

% Compute local quadrature points and weights
N_quadnodes = size(refTriQuad,1);
local_N_quadnodes = numVertices*N_quadnodes;
InternalQuadPoints = zeros(local_N_quadnodes,2);
InternalWeights = zeros(local_N_quadnodes,1);
BoundaryQuadPoints = [x_vert,y_vert;zeros(numVertices*(order-1),2)];
GLw = VemData.gaussLobattoQuad([1 end 2:end-1],2); % permutation of reference weights to respect boundary dofs ordering
BoundaryWeightsNormal = cell(2,1);
[BoundaryWeights, BoundaryWeightsNormal{:}] = deal(zeros(numVertices*order,1));
for e = 1:numVertices
    edgeVertices = mesh.edge(mesh.polygon(cellId).edges(e)).vertices;
    A = [mesh.vertex(edgeVertices(1)).x ; mesh.vertex(edgeVertices(1)).y];
    B = [mesh.vertex(edgeVertices(2)).x ; mesh.vertex(edgeVertices(2)).y];
    % map reference Gauss Lobatto points on current edge
    v_e = B-A;
    h_e = norm(v_e);
    n_e = [v_e(2);-v_e(1)]/h_e;
    n_e = n_e*(-1)^(dot(A-centroid,n_e)<0);
    edgeInternalPointsIndices = numVertices+(e-1)*(order-1)+(1:order-1);
    edgePointsIndices = [e, mod(e,numVertices)+1, edgeInternalPointsIndices];
    mappedEdgePoints = (A+B)'/2 + VemData.gaussLobattoQuad(:,1)*(B-A)'/2;
    mappedEdgeWeights = h_e/2*GLw;
    BoundaryWeights(edgePointsIndices) = ...
        BoundaryWeights(edgePointsIndices) + mappedEdgeWeights;
    BoundaryWeightsNormal{1}(edgePointsIndices) = ...
        BoundaryWeightsNormal{1}(edgePointsIndices) + mappedEdgeWeights*n_e(1);
    BoundaryWeightsNormal{2}(edgePointsIndices) = ...
        BoundaryWeightsNormal{2}(edgePointsIndices) + mappedEdgeWeights*n_e(2);
    BoundaryQuadPoints(edgeInternalPointsIndices,:) = mappedEdgePoints(2:end-1,:);
    % map reference Gauss points on the subtriangle corresponding to
    % current edge
    subTriangle = [A, B, centroid];
    subTriPointsIndices = (e-1)*N_quadnodes+(1:N_quadnodes);
    InternalQuadPoints(subTriPointsIndices,:) = ...
        mapFromRefTri(refTriQuad(:,1:2)', subTriangle)';
    InternalWeights(subTriPointsIndices) = ...
        polyarea(subTriangle(1,:),subTriangle(2,:))*refTriQuad(:,3);
end

% Compute local Vandermonde matrices
VanderInternal = Vander(order, alphas, centroid, diam, InternalQuadPoints);
VanderBoundary = Vander(order, alphas, centroid, diam, BoundaryQuadPoints);
VanderInternalDerivatives = VanderDerivatives(VanderInternal, alphas, diam);
VanderBoundaryDerivatives = VanderDerivatives(VanderBoundary, alphas, diam);

% Compute PiNabla projector matrix
cellPiNabla = computePiNabla(VemData, area, cellLaplacianMatrix, ...
    VanderInternalDerivatives, InternalWeights, ...
    VanderBoundaryDerivatives, BoundaryWeightsNormal, ...
    VanderBoundary, BoundaryWeights, VanderInternal);

% Compute monomials mass matrix and its Cholesky factor
cellMatrixH = VanderInternal'*diag(InternalWeights)*VanderInternal;
cellMatrixH_Cholkm1 = chol(cellMatrixH(1:VemData.N_km1,1:VemData.N_km1));

% Compute L^2 projection on k-1 polynomials
cellPi0km1 = computePi0km1(VemData, area, cellMatrixH, cellPiNabla, ...
    cellMatrixH_Cholkm1);

% Compute L^2 projection of derivatives on k-1 polynomials
cellPi0km1Grad = computePi0km1Grad(VemData, VanderBoundary, ...
    BoundaryWeightsNormal, area, cellDerivationMatrices, cellMatrixH_Cholkm1);

% Compute PiNabla projector returning VEM dofs
cellPiNabla_VemDofs = computePiNablaVemDofs(VemData, VanderBoundary, area, ...
    cellMatrixH, cellPiNabla);

% Coefficient values in quadrature points
F = forcingTerm(InternalQuadPoints(:,1),InternalQuadPoints(:,2));
cellDiffusionInQuadraturePoints = diffusionParameter(InternalQuadPoints(:,1), ...
    InternalQuadPoints(:,2));

% Stiffness matrix
cellDiffusionWeights =  cellDiffusionInQuadraturePoints.*InternalWeights;
cellMatrixH_Diffusion = VanderInternal(:,1:VemData.N_km1)' * ...
    diag(cellDiffusionWeights) * VanderInternal(:,1:VemData.N_km1);
cellStiffness_proj = cellPi0km1Grad{1}'*cellMatrixH_Diffusion*cellPi0km1Grad{1} + ...
    cellPi0km1Grad{2}'*cellMatrixH_Diffusion*cellPi0km1Grad{2};
cellStiffness_stab = max(cellDiffusionInQuadraturePoints) * ...
    (eye(nDOF,nDOF)-cellPiNabla_VemDofs)'*(eye(nDOF,nDOF)-cellPiNabla_VemDofs);
cellStiffness = cellStiffness_proj + cellStiffness_stab;

% Right-hand side
cellRightHandSide = cellPi0km1'*VanderInternal(:,1:VemData.N_km1)' * ...
    diag(InternalWeights)*F;

end
