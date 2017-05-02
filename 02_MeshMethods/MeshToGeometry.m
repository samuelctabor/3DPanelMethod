function [Geometry, Wake, NodesOut] = MeshToGeometry( Nodes, Elements )
    %MeshToGeometry Convert an imported mesh to a represtation suitable for
    %solver.
    %   Detailed explanation goes here
    
    if (0)
        % Only required if the geometry is not properly connected.
        [Nodes, Elements] = MergeNodes(Nodes,Elements,'wake','upper');
        [Nodes, Elements] = MergeNodes(Nodes,Elements,'wake','lower');
    end
    
    Surfaces = {Elements(~strcmp({Elements.Name},'wake')).Name};
    Geometry = SurfacesToGeometry(Nodes,Elements,Surfaces);
    
    Wake = SurfacesToGeometry(Nodes,Elements,{'wake'});
    
    NodesOut = [Nodes.X, Nodes.Y, Nodes.Z];
    
    % Identify the geometry elements that form the upper TE.
    [Wake.UpperTE, Wake.LowerTE] = FindTE(Geometry,Wake,NodesOut);
end

function Geometry=SurfacesToGeometry(Nodes,Elements,Surfaces)
    
    AllNodes = zeros(10000,4);
    Start = 1;
    
    for i=1:length(Surfaces)
        Idx = strcmp({Elements.Name},Surfaces{i});
        TheseNodes = Elements(Idx).Nodes;
        if size(TheseNodes,2)==3
            TheseNodes(:,4) = TheseNodes(:,3);
        end
        AllNodes(Start:Start+size(TheseNodes,1)-1,:) = TheseNodes;
        Start = Start + size(TheseNodes,1); 
    end
    
    AllNodes(Start:end,:) = [];
    
    X1 = Nodes.X(AllNodes(:,1));
    Y1 = Nodes.Y(AllNodes(:,1));
    Z1 = Nodes.Z(AllNodes(:,1));
    
    X2 = Nodes.X(AllNodes(:,2));
    Y2 = Nodes.Y(AllNodes(:,2));
    Z2 = Nodes.Z(AllNodes(:,2));
    
    X3 = Nodes.X(AllNodes(:,3));
    Y3 = Nodes.Y(AllNodes(:,3));
    Z3 = Nodes.Z(AllNodes(:,3));
    
    X4 = Nodes.X(AllNodes(:,4));
    Y4 = Nodes.Y(AllNodes(:,4));
    Z4 = Nodes.Z(AllNodes(:,4));
    
    Geometry.P1 = AllNodes(:,1);
    Geometry.P2 = AllNodes(:,2);
    Geometry.P3 = AllNodes(:,3);
    Geometry.P4 = AllNodes(:,4);
    
    Geometry.Patch = AllNodes;
end
