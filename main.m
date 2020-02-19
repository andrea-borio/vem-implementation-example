addpath mesh quadrature vem
%% Problem data
meshName = 'polymesher'; % see input folder
f = @(x,y) 32*(x.*(1-x)+y.*(1-y));
mu =  @(x,y) ones(size(x));
order = 5;

%% Initialize VEM
VemData = InitVemData(order);

%% Load mesh
load(fullfile('input',sprintf('%s.mat',meshName)),'mesh');
dofHandler = CreateDofHandler(mesh, VemData);

%% Get reference quadrature points
GaussTriangleQuad = TriGaussPoints(2*order+2);

%% Assemble matrix
triplets = zeros(3,dofHandler.numContributes);
rhs = zeros(dofHandler.numDofs,1);
pos = 1;
for e = 1:mesh.NP
    [cellStiffness,cellRightHandSide] = buildLocalSystem(mesh,e,VemData,...
        GaussTriangleQuad,f,mu);
    cellDofsCheck = dofHandler.cellGlobalDofIndices{e}>0;
    cellStiffnessDofsPositions = cellDofsCheck&cellDofsCheck';
    [dofRows, dofCols] = find(cellStiffnessDofsPositions);
    numCellTriplets = nnz(cellStiffnessDofsPositions);
    triplets(:,pos:pos+numCellTriplets-1) = [ ...
        dofHandler.cellGlobalDofIndices{e}(dofRows),...
        dofHandler.cellGlobalDofIndices{e}(dofCols),...
        cellStiffness(cellStiffnessDofsPositions)]';
    rhs(dofHandler.cellGlobalDofIndices{e}(cellDofsCheck)) = ...
        rhs(dofHandler.cellGlobalDofIndices{e}(cellDofsCheck)) + ...
        cellRightHandSide(cellDofsCheck);
    pos = pos+numCellTriplets;
end
A = sparse(triplets(1,:),triplets(2,:),triplets(3,:), ...
    dofHandler.numDofs, dofHandler.numDofs);
u = A\rhs;

%% Plot solution
sol = zeros(mesh.NV+mesh.NE+mesh.NP,1);
sol([dofHandler.pointDofIndices>0;dofHandler.edgeDofIndices(:)>0;...
    dofHandler.cellDofIndices(:)>0]) = u;
plotsolution(mesh,sol,VemData);
