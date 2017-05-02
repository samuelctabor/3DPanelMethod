function [ C,Phi ] = SourceInfluence( r_S,A )
%SourceInfluence Effect of source singularity
%   The singularity type is a 3D point source.
    if nargin<2
        A=1;
    end
    C = A*r_S./(rownorm(r_S).^3);
    Phi = A*1./rownormV(r_S);
end

