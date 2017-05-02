function [Gradient] = SurfaceGradientNodes( Geometry, Nodes, Scalar, Wake )
    %SurfaceGradient Compute surface gradient by Green-Gauss theorem
    %   Volume integral of scalar gradient is equal to the surface integral
    %   of the scalar
    % Note that this method can run into problems in the following
    % circumstances:
    % - Edges of connected regions. 
    % - Areas of high aspect ratio cells e.g. leading edge. Loss of
    % precision due to difference in weights being much larger in spanwise
    % diction than foil ordinate direction, while the opposite is true
    % for the scalar.
    
    nElements = size(Geometry.P1,1);
    nNodes = size(Nodes,1);
    
    [Vn,Area,Pc] = ComputeVn(Nodes,Geometry);
    Vl = (Nodes(Geometry.P1,:) + Nodes(Geometry.P2,:))/2 - Pc;
    Vl = Vl./rownorm(Vl);
    Vm = cross(Vn, Vl);
    
    % First find the values on the nodes.
    % Inverse-distance weighted average of neighbouring elements.    

    nMaxElems = 10;
    ConnectedElems = zeros(nNodes,nMaxElems);
    nConnectedElems = zeros(nNodes,1);
    for iNode = 1:nNodes
        idx = find(Geometry.Patch==iNode);
        [elem,~] = ind2sub(size(Geometry.Patch),idx);
        elem = unique(elem);
        ConnectedElems(iNode,1:length(elem)) = elem;
        nConnectedElems(iNode) = length(elem);
    end
    
    %
    % Find the distances.
    %
    % Form arrays of node and element positions and take distance.
    PElem = Pc(ConnectedElems(ConnectedElems~=0),:);
    AllNodes = repmat((1:nNodes)',1,nMaxElems);
    PNode = Nodes(AllNodes(ConnectedElems~=0),:);
    Dist = zeros(size(ConnectedElems));
    Dist(ConnectedElems~=0) = rownormV(PElem - PNode);
    
    %
    % Find the nodal values of the scalar.
    %
    % Compute weights based on inverse distance averaging.
    Weights = zeros(size(Dist));
    Weights(ConnectedElems~=0) = (1./Dist(ConnectedElems~=0));
    Weights = Weights./repmat(sum(Weights,2),1,nMaxElems);
    % Create array of element values by indexing into supplied vector.
    TempElemScalar = zeros(size(ConnectedElems));
    TempElemScalar(ConnectedElems~=0) = Scalar(ConnectedElems(ConnectedElems~=0));
    % Multiply by weights and add up to get nodal values.
    NodeValues = sum(Weights .* TempElemScalar,2);

    
    Gradient = zeros(nElements,2);
    % Surface integrate around each element. Contributions from edges on
    % trailing edge are ignored.
    IsTE = ismember(1:nNodes,Geometry.Patch) & ismember(1:nNodes,Wake.Patch);
    for iElem = 1:nElements
        % Form DCM
        DCM = [Vl(iElem,:);Vm(iElem,:);Vn(iElem,:)];
        % Edge normal and face value
        Integral = [0;0;0];
        nodeValues = NodeValues(Geometry.Patch(iElem,:));
        
        % Check for TE
        for iEdge = 1:4
            i1 = iEdge;
            i2 = mod(iEdge,4)+1;
            N1 = Geometry.Patch(iElem,i1);
            N2 = Geometry.Patch(iElem,i2);
            if N1~=N2
                if (IsTE(N1) && IsTE(N2))
                    nodeValues([i1,i2]) = Scalar(iElem);
                end
            end
        end
        
        % Do calc
        for iEdge = 1:4
            i1 = iEdge;
            i2 = mod(iEdge,4)+1;
            N1 = Geometry.Patch(iElem,i1);
            N2 = Geometry.Patch(iElem,i2);
            if N1~=N2
                EdgeValue = (nodeValues(i1) + nodeValues(i2))/2;
                
                % Resolve into element coordinate system
                Delta = Nodes(N2,:) - Nodes(N1,:);
                
                LocalDelta = DCM*Delta';
                EdgeNormal = cross(LocalDelta, [0;0;1]);
                
                Integral = Integral + EdgeNormal*EdgeValue;
            end
        end

        Gradient(iElem,:) = Integral(1:2) / Area(iElem);
    end
end

