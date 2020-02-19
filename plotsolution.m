function h = plotsolution(mesh,sol,VemData)
%PLOTSOLUTION plot the solution on the given mesh for specified VEM
%settings
%   h = plotsolution(mesh,sol,VemData) returns the handle to the generated
%   figure.

h = figure;
hold on

numEdgeInternalDofs = VemData.order-1;
numCellInternalDofs = VemData.N_km2;
T = zeros(sum(cellfun(@length,{mesh.polygon.vertices})),3);
X = [vertcat(mesh.vertex.x);zeros(mesh.NP,1)];
Y = [vertcat(mesh.vertex.y);zeros(mesh.NP,1)];
Z = [sol(1:mesh.NV);zeros(mesh.NP,1)];
previous_triangles = 0;
edges_to_plot = true(mesh.NE,1);
for i=1:mesh.NP
    local_nodes = [mesh.polygon(i).vertices,mesh.polygon(i).vertices(1)];
    for k = 1:length(mesh.polygon(i).edges)
        edgeIndex = mesh.polygon(i).edges(k);
        if edges_to_plot(edgeIndex)
            edges_to_plot(edgeIndex) = false;
            line_nodes = mesh.edge(edgeIndex).vertices;
            plot3(X(line_nodes),Y(line_nodes),sol(line_nodes),'k');
        end
    end
    centroid = computeCellCentroid(X(mesh.polygon(i).vertices),...
        Y(mesh.polygon(i).vertices),...
        polyarea(X(mesh.polygon(i).vertices),Y(mesh.polygon(i).vertices)));
    X(mesh.NV+i) = centroid(1);
    Y(mesh.NV+i) = centroid(2);
    T(previous_triangles+(1:length(local_nodes)),:) = ...
        [local_nodes',local_nodes([2:end,1])' , ...
        (mesh.NV+i)*ones(length(local_nodes),1)];
    if VemData.order == 1
        Z(mesh.NV+i) = mean(sol(local_nodes));
    else
        Z(mesh.NV+i) = sol(mesh.NV+mesh.NE*numEdgeInternalDofs+(i-1)*numCellInternalDofs+1);
    end
    previous_triangles = previous_triangles + length(local_nodes);
end
trisurf(T,X,Y,Z,'linestyle','none');
view(25,45)
axis tight
