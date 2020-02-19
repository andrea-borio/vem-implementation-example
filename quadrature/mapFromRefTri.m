function [coords,jacobian_matrix] = mapFromRefTri (hat_coords, vertexes)
%mapFromRefTri map from reference triangle
%   COORDS = mapFromRefTri (HAT_COORDS, VERTEXES) maps the given vector of
%   2D column vectors (i.e. a 2-by-N matrix) HAT_COORDS, defined on the
%   reference triangle, on the triangle defined by the 2-by-3 matrix
%   VERTEXES. 
%
%   [COORDS,JACOBIAN_MATRIX] = mapFromRefTri (HAT_COORDS, VERTEXES) also
%   returns the jacobian matrix of the transformation.

jacobian_matrix = [vertexes(1,mod(1,3)+1)-vertexes(1,1) ...
    vertexes(1,mod(2,3)+1)-vertexes(1,1); ...
    vertexes(2,mod(1,3)+1)-vertexes(2,1) ...
    vertexes(2,mod(2,3)+1)-vertexes(2,1) ];

coords = zeros(size(hat_coords));
for i = 1: size(hat_coords,2)
    coords(:,i) = vertexes(:,1) + jacobian_matrix * hat_coords(:,i);
end

end
