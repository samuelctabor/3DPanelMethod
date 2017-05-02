function TestMesh( )
    %TestMesh Test function for reading external meshes and running
    %solution.
    %   Detailed explanation goes here
    
    close all;

    % Read mesh.
    [Nodes,Elements] = ReadMesh('G:/ST2_files/Info/2410_6.inp');
    Elements = Elements(cellfun(@(s)~strcmp(s,'tip'),{Elements.Name}));
    PlotMesh(Nodes,Elements);

    % Create geometry and wake.
    [Geometry,Wake, Nodes] = MeshToGeometry(Nodes,Elements);
    
    % Change axes.
    Nodes = ([1, 0, 0;...
              0, 0,-1;...
              0, 1, 0] * Nodes')';
         
    % Scale
    Nodes = Nodes*1/(max(Nodes(Geometry.Patch(:),1))-min(Nodes(Geometry.Patch(:),1)));
    Nodes(:,2) = -Nodes(:,2)*20/(max(Nodes(Geometry.Patch(:),2))-min(Nodes(Geometry.Patch(:),2)));
         
    % Move wake nodes far downstream.
    MoveNode = ismember(1:size(Nodes,1),Wake.Patch) & ~ismember(1:size(Nodes,1),Geometry.Patch);
    Nodes(MoveNode,1) = Nodes(MoveNode,1) + 1000;
    
    %
    % Set up BCs and run solution.
    %
    alpha = deg2rad(5);
    
    Vfs = [cos(alpha), 0, sin(alpha);...
           0,          1, 0;...
           sin(alpha), 0, cos(alpha)] * [1;0;0];
       
    % Align elements with freestream
    Geometry = AlignWithFreestream(Geometry,Nodes,Vfs);
    Geometry = ReverseNormals(Geometry);
    
    Solution = RunSolution(Geometry, Vfs, Nodes, Wake);
   
    Gradient = SurfaceGradientEdges(Geometry,Nodes,Solution.mu,Wake);
    %Gradient = SurfaceGradientNodes(Geometry,Nodes,Solution.mu,Wake);
    Cp       = CpFromSurfaceGradient(Geometry,Nodes,Gradient,Vfs);
   
    
        
    [Vn,~,Pc] = ComputeVn(Nodes,Geometry);
    Vl = (Nodes(Geometry.P1,:) + Nodes(Geometry.P2,:))/2 - Pc;
    Vl = Vl./rownorm(Vl);
    Vm = cross(Vn, Vl);
    ElemIdx = Pc(:,2)>9.9 & Pc(:,2)<10.1;
    %ElemIdx = Pc(:,2)>10 & Pc(:,2)<10.4;
    
    figure,
    subplot(3,1,1);
    plot(Pc(ElemIdx,1),Pc(ElemIdx,3),'bo'); title('Co-ords');
    %text(Pc(ElemIdx,1),Pc(ElemIdx,3),arrayfun(@(s) sprintf(' %i',s),Elems(ElemIdx),'UniformOutput',false));
    
    subplot(3,1,2);
    plot(Pc(ElemIdx,1),Solution.mu(ElemIdx),'bo'); title('\mu');
    %text(Pc(ElemIdx,1),Solution.mu(ElemIdx),arrayfun(@(s) sprintf(' %i',s),Elems(ElemIdx),'UniformOutput',false));
    
    subplot(3,1,3);
    plot(Pc(ElemIdx,1),Cp(ElemIdx),'bo'); title('Cp');
    %text(Pc(ElemIdx,1),Solution.mu(ElemIdx),arrayfun(@(s) sprintf(' %i',s),Elems(ElemIdx),'UniformOutput',false));
    

    PlotPanelCoordinateSystems(Geometry,Nodes,ElemIdx);
    
    ConfigurationSurfacePlot(Geometry,Nodes,Wake,Solution.mu);
    ConfigurationSurfacePlot(Geometry,Nodes,Wake,Gradient);
    ConfigurationSurfacePlot(Geometry,Nodes,Wake,Cp);
end

