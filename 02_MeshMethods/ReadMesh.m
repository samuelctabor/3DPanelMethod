function [MeshNodes,MeshElements] = ReadMesh( FileName )
    %ReadMesh Read a .inp mesh file
    %   Detailed explanation goes here
    fid = fopen(FileName);
    Lines = textscan(fid,'%s','Delimiter','\n');
    Lines = Lines{1};
    
    % Remove comment lines (**)
    Idx = cellfun(@(s) length(s)<2 || ~strcmp(s(1:2),'**'), Lines,'UniformOutput',true);
    Lines = Lines(Idx);
    
    % Find useful lines.
    IsHeading = cellfun(@(s) ~isempty(strfind(s,'*NODE')) | ~isempty(strfind(s,'*ELEMENT')) | ~isempty(strfind(s,'*SHELL')), Lines);
    HeadingLines = find(IsHeading);
    HeadingLines(end+1) = length(Lines);
    
    iNode = 1;
    iElement = 1;
    for i=1:length(HeadingLines)-1
        if strcmp(Lines{HeadingLines(i)}, '*NODE')
            MeshNodes(iNode) = ParseNodes(Lines(HeadingLines(i)+1: HeadingLines(i+1) - 1));
            iNode = iNode+1;
        elseif ~isempty(strfind(Lines{HeadingLines(i)}, '*ELEMENT, TYPE=S3'))
            MeshElements(iElement) = ParseElements(Lines(HeadingLines(i): HeadingLines(i+1) - 1),3);
            iElement = iElement+1;
        elseif ~isempty(strfind(Lines{HeadingLines(i)}, '*ELEMENT, TYPE=S4'))
            MeshElements(iElement) = ParseElements(Lines(HeadingLines(i): HeadingLines(i+1) - 1),4);
            iElement = iElement+1;
        end
    end
    fclose(fid);
   
    % Renumber nodes
    if max(MeshNodes.Id)>length(MeshNodes.X)
        Map(MeshNodes.Id) = 1:length(MeshNodes.X);
        
        for iElem=1:length(MeshElements)
            MeshElements(iElem).Nodes = Map(MeshElements(iElem).Nodes);
        end
    end
end

function NodesStruct = ParseNodes(Lines)
    Nodes = cellfun(@(s)sscanf(s,'%i, %f, %f, %f'), Lines, 'UniformOutput',false);
    Nodes = cell2mat(cellfun(@(a) a',Nodes, 'UniformOutput',false));
    NodesStruct.Id = Nodes(:,1);
    NodesStruct.X = Nodes(:,2);
    NodesStruct.Y = Nodes(:,3);
    NodesStruct.Z = Nodes(:,4);
end

function ElementsStruct = ParseElements(Lines,nNodes)
    % Find the surface name.
    Header = Lines{1};
    substr = strsplit(Header,';');
    ElementsStruct.Name = substr{end};
    Lines = Lines(2:end);
    
    if nNodes==3
        Elements = cellfun(@(s)sscanf(s,'%i, %i, %i, %i'), Lines, 'UniformOutput',false);
    elseif nNodes==4
        Elements = cellfun(@(s)sscanf(s,'%i, %i, %i, %i, %i'), Lines, 'UniformOutput',false);
    end
    Elements = cell2mat(cellfun(@(a) a',Elements, 'UniformOutput',false));
    
    ElementsStruct.Id = Elements(:,1);
    ElementsStruct.Nodes = Elements(:,2:1+nNodes);
end
