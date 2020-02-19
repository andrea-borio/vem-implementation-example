function VemData = InitVemData(order)
%INITVEMDATA global parameters for the VEM discretization of the Laplace
%problem.
%   VEMDATA = INITVEMDATA(ORDER) returns a structure containing the
%   properties required by the VEM discretization of given ORDER that are
%   independent of the polygons.

VemData = struct('order', order, ...
    'monomialExponents', monomialExponents(order), ...
    'derivationMatrices', { cell(2,1) }, ...
    'laplacianMatrix', [], ...
    'N_k', 0, ...
    'N_km1', 0, ...
    'N_km2', 0, ...
    'gaussLobattoQuad', lglnodes(2*order-1));

VemData.N_k = size(VemData.monomialExponents,1);
VemData.N_km1 = VemData.N_k - order - 1;
VemData.N_km2 = VemData.N_km1 - order;

% reference derivation matrices (to be multiplied by 1/h_E, for all E).
VemData.derivationMatrices = { zeros(VemData.N_k) , zeros(VemData.N_k) };
for m=2:VemData.N_k
    if VemData.monomialExponents(m,1)>0
        VemData.derivationMatrices{1}(m - VemData.monomialExponents(m,1) - ...
            VemData.monomialExponents(m,2),m) = VemData.monomialExponents(m,1);
    end
    if VemData.monomialExponents(m,2)>0
        VemData.derivationMatrices{2}(m - VemData.monomialExponents(m,1) - ...
            VemData.monomialExponents(m,2) - 1,m) = VemData.monomialExponents(m,2);
    end
end

% reference laplacian matrix (to be multiplied by 1/h_E^2, for all E).
VemData.laplacianMatrix = VemData.derivationMatrices{1}^2 + ...
    VemData.derivationMatrices{2}^2;

end
