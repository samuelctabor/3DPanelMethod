function [Geometry, nxp, nyp] = PointsToGeometry(Nodes, n1, n2, ClosedDim)
    % Convert an array of nodes to an appropriate surface definition.
    % Nodes: Nx3 array of nodal xyz coordinates.
    % n1: Number of rows.
    % n2: Number of columns. n1xn2=N.
    % ClosedDim: Indicator of dimension to close the surface in. 
    %           0: Don't close. All elements have unique nodes.
    %           1: Close in row direction. A row of elements will link
    %           the first and last rows of nodes.
    %           2: Close in column direction. A column of elements will
    %           link the first and last columns of nodes.
    
    NodeIdx = reshape(1:length(Nodes),[n1,n2]);
    
    Geometry.P3 = NodeIdx(1:n1-1,1:n2-1);
    Geometry.P4 = NodeIdx(1:n1-1,2:n2  );
    Geometry.P1 = NodeIdx(2:n1  ,2:n2  );
    Geometry.P2 = NodeIdx(2:n1  ,1:n2-1);
    
    if ClosedDim==1
        Geometry.P3 = [Geometry.P3; NodeIdx(n1,1:n2-1)];
        Geometry.P4 = [Geometry.P4; NodeIdx(n1,2:n2  )];
        Geometry.P1 = [Geometry.P1; NodeIdx(1 ,2:n2  )];
        Geometry.P2 = [Geometry.P2; NodeIdx(1 ,1:n2-1)];
        nxp = n1;
        nyp = n2-1;
    elseif ClosedDim==2
        Geometry.P3 = [Geometry.P3, NodeIdx(1:n1-1,n2)];
        Geometry.P4 = [Geometry.P4, NodeIdx(1:n1-1,1 )];
        Geometry.P1 = [Geometry.P1, NodeIdx(2:n1  ,1 )];
        Geometry.P2 = [Geometry.P2, NodeIdx(2:n1  ,n2)];
        nxp = n1-1;
        nyp = n2;
    else
        nxp = n1-1;
        nyp = n2-1;
    end
    
    Geometry.Patch = [Geometry.P1(:), Geometry.P2(:), Geometry.P3(:), Geometry.P4(:)];
    
    Geometry.P1 = Geometry.P1(:);
    Geometry.P2 = Geometry.P2(:);
    Geometry.P3 = Geometry.P3(:);
    Geometry.P4 = Geometry.P4(:);
    
    Geometry.n1 = nxp;
    Geometry.n2 = nyp;

%     vn = cross(Geometry.P3-Geometry.P1,Geometry.P4-Geometry.P2);
%     Geometry.Area = rownormV(vn)/2;
%     Geometry.vn = vn./rownorm(vn);
% 
%     Geometry.PC = (Geometry.P1+Geometry.P2+Geometry.P3+Geometry.P4)/4;
%     
%     Geometry.n1 = n1;
%     Geometry.n2 = n2;
end

