function [ NodesOut, ElementsOut ] = MergeNodes(Nodes, Elements, Surf1, Surf2, tol)
    %MergeNodes Combine nodes that are coincident
    %   All nodes in Surf1 will be merged with those in Surf2 that are
    %   coindicent.
    NodeCoord = [Nodes.X,Nodes.Y,Nodes.Z];
    if nargin<5
        tol = 1e-10;
    end

    if ischar(Surf1)
        Idx1 = find(strcmp(Surf1,{Elements.Name}));
        Idx2 = find(strcmp(Surf2,{Elements.Name}));
    else
        Idx1 = Surf1;
        Idx2 = Surf2;
    end
    NodeIdx1 = unique(Elements(Idx1).Nodes);
    NodeIdx2 = unique(Elements(Idx2).Nodes);
    Nodes1 = NodeCoord(NodeIdx1(:),:);
    Nodes2 = NodeCoord(NodeIdx2(:),:);

    Dist = squareform(pdist([Nodes1;Nodes2]));

    % We only want the distance from each node in Surf1 to those in Surf2,
    % i.e. not the distances between nodes in Surf1 to those in Surf1 or
    % those in Surf2 to those in Surf2.
    SubDist = Dist(1:length(Nodes1),(length(Nodes1)+1):end);

    % Find where they are close.
    IdxClose = find(SubDist<tol);

    % Corresponding node numbers.
    [i,j] = ind2sub(size(SubDist),IdxClose);
    MergeNode1 = NodeIdx1(i);
    MergeNode2 = NodeIdx2(j);
    
    % Make sure we don't merge nodes that are already the same!
    IsSame = MergeNode1 == MergeNode2;
    MergeNode1 = MergeNode1(~IsSame);
    MergeNode2 = MergeNode2(~IsSame);

    % Remove the duplicate nodes.
    NodesOut = Nodes;

    NodesOut.Id(MergeNode1) = [];
    NodesOut.X(MergeNode1) = [];
    NodesOut.Y(MergeNode1) = [];
    NodesOut.Z(MergeNode1) = [];

    % Update all node references to reflect the change in indexing.
    % - Change references of duplicate nodes to point at partners
    % - Change all node references to reflect removal of duplicates.
    Map = 1:size(Nodes.X,1);
    Map(MergeNode1) = MergeNode2;
    [~,~,Map] = unique(Map);

    ElementsOut = Elements;
    for iElem=1:length(Elements)
        ElementsOut(iElem).Nodes = Map(ElementsOut(iElem).Nodes);

        if any(any(ElementsOut(iElem).Nodes==0))
            error('Merging error');
        end
    end
end

