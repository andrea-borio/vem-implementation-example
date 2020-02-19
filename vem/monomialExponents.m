function expo = monomialExponents(order)
%MONOMIALEXPONENTS
%   EXPO = monomialExponents(ORDER) generates the table of exponents of the
%   2D monomials of given ORDER

numMonomials = (order+1)*(order+2)/2;
expo=zeros(numMonomials,2);

for i=2:numMonomials
    if expo(i-1,1)==0
        expo(i,1)=expo(i-1,2)+1;
    else
        expo(i,:)=[expo(i-1,1)-1 expo(i-1,2)+1];
    end
end

end
