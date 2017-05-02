function Cp=CpFromSurfaceGradient(Elements,Nodes,Gradient,Ufs)
% CpFromSurfaceGradient Calculates the surface Cp from the surface
% circulation gradient and freestream velocity.
%   Ul = Ufs.l + 4*pi*dmu/dl
%   Um = Ufs.m + 4*pi*dmu/dm
%   Ut = sqrt(Ul^2+Um^2)
%   Cp = 1 - Ut^2/Ufs^2

    [Vn,~,Pc] = ComputeVn(Nodes,Elements);
    Vl = (Nodes(Elements.P1,:) + Nodes(Elements.P2,:))/2 - Pc;
    Vl = Vl./rownorm(Vl);
    Vm = cross(Vn, Vl);
    
    Ufs_vec = repmat(Ufs',[length(Elements.P1),1]);
    
    Ufs_l = dot(Vl,Ufs_vec,2);
    Ufs_m = dot(Vm,Ufs_vec,2);
    
    U_l = Ufs_l - 4*pi*Gradient(:,1);
    U_m = Ufs_m - 4*pi*Gradient(:,2);
    
    Cp = 1 - (U_l.^2 + U_m.^2)/(norm(Ufs)^2);
end