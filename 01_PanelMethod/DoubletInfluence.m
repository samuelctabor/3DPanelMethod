function [C,Phi] = DoubletInfluenceNew( r_S, n,A )
%DoubletInfluence Function for computing effect of 3D point doublet.
%   r_S is the vector from doublet to point.
%   nD is the doublet heading vector (only x is supported).
    if nargin<3
        A=1;
    end

    P = r_S;
    PN = dot(P,n,2);
    
    C = A*(3*repmat(PN,1,3).*P - rownorm(P).^2.*n)./(rownorm(P).^5);
    
    Phi = A*PN./rownormV(P).^3;
    
    
    % 
    %Phi = -A*
end

