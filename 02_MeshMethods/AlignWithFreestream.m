function GeometryOut = AlignWithFreestream(Geometry,Nodes,Vfs)
    %AlignWithFreestream Re-order nodes so that the l vector is aligned
    %with freestream.
    %   Detailed explanation goes here
    
    Vfs = repmat(Vfs',size(Geometry.P1,1),1);
    
    [Vn,~,Pc] = ComputeVn(Nodes,Geometry);
    Vl = (Nodes(Geometry.P1,:) + Nodes(Geometry.P2,:))/2 - Pc;
    Vl = Vl./rownorm(Vl);
    Vm = cross(Vn, Vl);
    
    IsOK = abs(dot(Vl,Vfs,2)) > abs(dot(Vm,Vfs,2));
    
    GeometryOut = Geometry;
    GeometryOut.P1(~IsOK) = Geometry.P2(~IsOK);
    GeometryOut.P2(~IsOK) = Geometry.P3(~IsOK);
    GeometryOut.P3(~IsOK) = Geometry.P4(~IsOK);
    GeometryOut.P4(~IsOK) = Geometry.P1(~IsOK);
    
    GeometryOut.Patch = [Geometry.P1(:),Geometry.P2(:),Geometry.P3(:),Geometry.P4(:)];
    
end

