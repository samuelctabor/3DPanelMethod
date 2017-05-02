function [ h ] = PlotMesh( Nodes,Elements, Scalar)
    %PlotMesh Plot a mesh consisting of nodes and elements structs.
    %   Detailed explanation goes here
    
    h = figure;
    hold on;
    
    colors = {'b','g','r','m'};
    Legend = {};
    for iEl = 1:length(Elements)
        PatchX = zeros(size(Elements(iEl).Nodes));
        PatchY = zeros(size(Elements(iEl).Nodes));
        PatchZ = zeros(size(Elements(iEl).Nodes));
        for iN = 1:size(Elements(iEl).Nodes,2)
            PatchX(:,iN) =  Nodes(1).X(Elements(iEl).Nodes(:,iN));
            PatchY(:,iN) = -Nodes(1).Z(Elements(iEl).Nodes(:,iN));
            PatchZ(:,iN) =  Nodes(1).Y(Elements(iEl).Nodes(:,iN));
        end
        
        if nargin<3
            patch(PatchX',PatchY',PatchZ',colors{iEl},'EdgeAlpha',0.5);
            Legend{iEl} = Elements(iEl).Name;
        else
            patch(PatchX',PatchY',PatchZ',Scalar);
            Legend = 'mu';
        end
        
    end
    legend(Legend);
    axis('equal');
end

