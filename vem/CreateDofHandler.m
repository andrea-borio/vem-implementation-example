function dofHandler = CreateDofHandler(mesh, VemData)
%CREATEDOFHANDLER creates the structure to handle degrees of freedom.
%   DOFHANDLER = CREATEDOFHANDLER(MESH, VEMDATA) returns the degrees of
%   freedom handler for the mesh defined by MESH and the VEM space defined
%   by VEMDATA
%
%   See also INITVEMDATA

numCellInternalDofs = VemData.N_km2;
numEdgeInternalDofs = VemData.order-1;

dofHandler = struct('numDofs', 0, 'numDirichletDofs', 0, 'numContributes', 0, ...
    'numDirichletContributes', 0, ...
    'cellDofIndices', zeros(numCellInternalDofs, mesh.NP), ...
    'edgeDofIndices', zeros(numEdgeInternalDofs, mesh.NE), ...
    'pointDofIndices', zeros(mesh.NP,1), ...
    'cellGlobalDofIndices', {cell(mesh.NP,1)});

for i = 1:mesh.NV
    if mesh.vertex(i).marker == 0
        dofHandler.numDofs = dofHandler.numDofs+1;
        dofHandler.pointDofIndices(i) = dofHandler.numDofs;
    else
        dofHandler.numDirichletDofs = dofHandler.numDirichletDofs+1;
        dofHandler.pointDofIndices(i) = -dofHandler.numDirichletDofs;
    end
end

if VemData.order>1
    for i = 1:mesh.NE
        if mesh.edge(i).marker == 0
            for j = 1:numEdgeInternalDofs
                dofHandler.numDofs = dofHandler.numDofs+1;
                dofHandler.edgeDofIndices(j,i) = dofHandler.numDofs;
            end
        else
            for j = 1:numEdgeInternalDofs
                dofHandler.numDirichletDofs = dofHandler.numDirichletDofs+1;
                dofHandler.edgeDofIndices(j,i) = -dofHandler.numDirichletDofs;
            end
        end
    end
    for i = 1:mesh.NP
        for j = 1:numCellInternalDofs
            dofHandler.numDofs = dofHandler.numDofs+1;
            dofHandler.cellDofIndices(j,i) = dofHandler.numDofs;
        end
    end
end

for i = 1:mesh.NP
    numCellVertices = length(mesh.polygon(i).vertices);
    numCellBoundaryDofs = VemData.order*numCellVertices;
    numCellDofs = numCellBoundaryDofs+numCellInternalDofs;
    dofHandler.cellGlobalDofIndices{i} = zeros(numCellDofs,1);
    dofHandler.cellGlobalDofIndices{i}(1:numCellVertices) = ...
        dofHandler.pointDofIndices(mesh.polygon(i).vertices);
    if VemData.order > 1
        dofHandler.cellGlobalDofIndices{i}(numCellVertices+1:numCellBoundaryDofs) = ...
            dofHandler.edgeDofIndices(:,mesh.polygon(i).edges);
        dofHandler.cellGlobalDofIndices{i}(numCellBoundaryDofs+1:numCellDofs) = ...
            dofHandler.cellDofIndices(:,i);
    end
end

cellDofChecks = cell2mat(cellfun(@(x)sum([x>0,x<0]), ...
    dofHandler.cellGlobalDofIndices,'UniformOutput',false));
dofHandler.numContributes = sum(cellDofChecks(:,1).^2);
dofHandler.numDirichletContributes = sum(prod(cellDofChecks,2));



end
