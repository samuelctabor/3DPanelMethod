function GeometryOut = ReverseNormals(Geometry)
    %ReverseNormals Reverse the normals of all elements
    
    GeometryOut.P1 = Geometry.P2;
    GeometryOut.P2 = Geometry.P1;
    GeometryOut.P3 = Geometry.P4;
    GeometryOut.P4 = Geometry.P3;
    
    GeometryOut.Patch = [GeometryOut.P1(:),GeometryOut.P2(:),GeometryOut.P3(:),GeometryOut.P4(:)];
end

