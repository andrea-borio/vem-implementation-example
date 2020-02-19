function centroid = computeCellCentroid(x,y,area)
%COMPUTECELLCENTROID computes the centroid of a polygon given its vertices
%and area
%   computeCellCentroid(x,y,area)

x_1 = x([2:end 1]);
y_1 = y([2:end 1]);
centroid = [sum((x+x_1).*(x.*y_1-x_1.*y))/(6*area) ; 
    sum((y+y_1).*(x.*y_1-x_1.*y))/(6*area)];

end
