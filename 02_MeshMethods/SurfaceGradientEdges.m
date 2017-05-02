function [Gradient] = SurfaceGradientEdges( Geometry, Nodes, Scalar, Wake)
    %SurfaceGradient Compute surface gradient by Green-Gauss theorem
    %   Volume integral of scalar gradient is equal to the surface integral
    %   of the scalar
    
    % Find the values on the edges.
    % First we need to construct the edges.
    
    nElements = size(Geometry.P1,1);
    nNodes = size(Nodes,1);
    
    EdgeNodes = zeros(1000,2);
    ElemEdges = zeros(size(Geometry.P1));
    Start = 1;
    IsTri = arrayfun(@(i) length(unique(Geometry.Patch(i,:))), 1:nElements) == 3; 
    
    [Vn,Area,Pc] = ComputeVn(Nodes,Geometry);
    Vl = (Nodes(Geometry.P1,:) + Nodes(Geometry.P2,:))/2 - Pc;
    Vl = Vl./rownorm(Vl);
    Vm = cross(Vn, Vl);
    
    for iEl = 1:nElements
        
        NBNodeSet3 = [1,2; 2,3; 3,1];
        NBNodeSet4 = [1,2; 2,3; 3,4; 4,1];
        
        if IsTri(iEl)
            % 3 nodes
            NBNodeSet = NBNodeSet3;
        else
            NBNodeSet = NBNodeSet4;
        end
        
        NewEdges = reshape(Geometry.Patch(iEl,NBNodeSet')', size(NBNodeSet'))';
        ElemEdges(iEl,1:length(NewEdges)) = [Start:Start+length(NewEdges)-1];
        EdgeNodes(Start:Start+length(NewEdges)-1,:) = NewEdges;
        Start = Start+length(NewEdges);
    end
    EdgeNodes(Start:end,:)=[];
    
    %Remove duplicates
    [EdgeNodes, IdxSort] = sort(EdgeNodes,2);
    [UniqueEdgeNodes,~,IC] = unique(EdgeNodes,'rows');
    
    ElemUniqueEdges = zeros(size(ElemEdges));
    ElemUniqueEdges(ElemEdges~=0) = IC(ElemEdges(ElemEdges~=0));

    NodeID = 1:nNodes;
    IsNodeTE = ismember(NodeID,Geometry.Patch) & ismember(NodeID,Wake.Patch);
    IsEdgeTE = ismember(UniqueEdgeNodes(:,1),NodeID(IsNodeTE)) & ismember(UniqueEdgeNodes(:,2),NodeID(IsNodeTE));

    % Now we know how everything connects.
    
    % For now assume the edge value is just the average of the face values
    EdgeValues         = zeros(size(UniqueEdgeNodes,1),1);
    for iEdge = 1:length(UniqueEdgeNodes)
        % Find the adjoining elements.
        Elems = any(ElemUniqueEdges==iEdge,2);

        % Find the edge value.
        if sum(Elems)==1
            EdgeValues(iEdge) = Scalar(Elems);
        else
            % This is the distance-weighted average, where the distace is that
            % from the element centres to a certain point on the edge.
            % This point is the intersection of the edge and a plane defined
            % by the element centres and the average of the two element normals.
            meanNormal = mean(Vn(Elems,:));
            planeNormal = cross(meanNormal,diff(Pc(Elems,:)));
            pc = Pc(Elems,:);
            pIntersection = IntersectPlaneLine(planeNormal,pc(1,:),diff(Nodes(UniqueEdgeNodes(iEdge,:),:)),Nodes(UniqueEdgeNodes(iEdge,1),:));

            % Find distances of each centre to point on edge.
            d = rownormV(Pc(Elems,:) - [pIntersection;pIntersection]);

            EdgeValues(iEdge)  = sum(flip(d/sum(d)).*Scalar(Elems));
        end
    end
    
    Delta = Nodes(UniqueEdgeNodes(:,2),:) - Nodes(UniqueEdgeNodes(:,1),:);

    Gradient = zeros(nElements,2);
    Cont = zeros(nElements,4);
    % Surface integrate around each element.
    for iElem = 1:nElements
        % Form DCM
        DCM = [Vl(iElem,:);Vm(iElem,:);Vn(iElem,:)];
        % Edge normal and face value
        Integral = [0;0;0];
        for iEdge = 1:size(ElemUniqueEdges,2)
            EdgeId = ElemUniqueEdges(iElem,iEdge);
            if EdgeId~=0
                % Resolve into element coordinate system
                LocalDelta = DCM*Delta(EdgeId,:)';
                EdgeNormal = cross(LocalDelta, [0;0;1]);
                
                % Check direction
                Reversed = UniqueEdgeNodes(EdgeId,1) ~= Geometry.Patch(iElem,iEdge);
                Sign = 1 - 2*Reversed;
                
                if IsEdgeTE(EdgeId)
                    edgeValue = Scalar(iElem);
                else
                    edgeValue = EdgeValues(EdgeId);
                end
                
                Integral = Integral + Sign*EdgeNormal*edgeValue;
                Cont(iElem,iEdge) = Sign*EdgeNormal(1)*EdgeValues(EdgeId);
            end
        end

        Gradient(iElem,:) = Integral(1:2) / Area(iElem);
    end
end

function pIntersect = IntersectPlaneLine(vPlane, pPlane, vLine, pLine)
    % See https://en.wikipedia.org/wiki/Line%E2%80%93plane_intersection
    d = dot(pPlane - pLine, vPlane) ./ dot(vLine, vPlane);
    pIntersect = pLine + d*vLine;
end
