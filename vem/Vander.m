function VandermondeMatrix = Vander(order, exponents, centroid, diameter, points)
%VANDER Vandermonde matrix of scaled monomials
%   VANDERMONDEMATRIX = VANDER(ORDER, EXPONENTS, CENTROID, DIAMETER,
%   POINTS) computes the Vandermonde matrix of the monomials of order ORDER
%   whose exponents are in EXPONENTS, evaluated at POINTS. CENTROID is the
%   centroid of the polygon and DIAMETER its diameter.
%
%   See also MONOMIALEXPONENTS

numMonomials = size(exponents,1);
numPoints = size(points,1);

VPartial_X = ones(numPoints,order+1);
VPartial_X(:,2) = (points(:,1)-centroid(1))./diameter;
VPartial_Y = ones(numPoints,order+1);
VPartial_Y(:,2) = (points(:,2)-centroid(2))./diameter;
for j = 3:order+1
    VPartial_X(:,j) = VPartial_X(:,j-1).*VPartial_X(:,2);
    VPartial_Y(:,j) = VPartial_Y(:,j-1).*VPartial_Y(:,2);
end

VandermondeMatrix = zeros(numPoints, numMonomials);
for j = 1:numMonomials
    VandermondeMatrix(:,j) = VPartial_X(:,1+exponents(j,1)) .* ...
        VPartial_Y(:,1+exponents(j,2));
end

