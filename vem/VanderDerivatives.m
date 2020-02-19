function VandermondeMatrices = VanderDerivatives(Vander, exponents, diameter)
%VANDERDERIVATIVES Vandermonde matrices of derivatives of scaled monomials
%   VANDERMONDEMATRICES = VANDERDERIVATIVES(VANDER, EXPONENTS, DIAMETER)
%   returns a 2-by-1 cell array containing the Vandermonde matrices of
%   derivatives of scaled monomials whose Vandermonde matrix is VANDER and
%   whose exponents are in EXPONENTS. DIAMETER is the diameter of the
%   polygon
%
%   See also MONOMIALEXPONENTS, VANDER

[numQuadraturePoints, numMonomials] = size(Vander);
VandermondeMatrices = cell(2,1);
[VandermondeMatrices{:}] = deal(zeros(numQuadraturePoints,numMonomials));
jDerivatives = [0;0];
for j = 2:numMonomials
    jDerivatives(1) = j - sum(exponents(j,:));
    jDerivatives(2) = jDerivatives(1) - 1;
    for d = 1:2
        if exponents(j,d)>0
            VandermondeMatrices{d}(:,j) = (1/diameter)*exponents(j,d)*...
                Vander(:,jDerivatives(d));
        end
    end
end



end
