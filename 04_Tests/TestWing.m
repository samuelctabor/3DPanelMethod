function TestWing
    %Test_Wing Summary of this function goes here
    %   Detailed explanation goes here
    addpath('G:/ST2_files/MATLAB/FoamCutter');
    close all;
    
    %
    % Create the wing geometry.
    %
    
    L = 20;
    ny = 20;
    Downstream = 500;    
    alpha = degtorad(5);
    Chord = 1;
    
    Foil.designation = '2410';
    Foil.n           = 100;
    Foil.HalfCosineSpacing = 1;
    Foil.wantFile = 0;
    Foil.is_finiteTE = 0;
    
    NacaFoil = naca4gen(Foil);
    
    NacaFoil.x = Chord*[NacaFoil.xU(1:(end-1));NacaFoil.xL];
    NacaFoil.z = Chord*[NacaFoil.zU(1:(end-1));NacaFoil.zL];
    
    nx = length(NacaFoil.x)-1;
    Foil.Name = Foil.designation;
    Foil.Points = [NacaFoil.x(1:nx),NacaFoil.z(1:nx)];
    
    
    y = fliplr(L * cos(pi * (0:(ny-1)) / (ny-1)));
    
    X = repmat(Foil.Points(:,1), 1, ny);
    Y = repmat(               y, nx,1);
    Z = repmat(Foil.Points(:,2), 1, ny);
    
    % Pinch the tips closed.
%     Z(:,1) = 0;
%     Z(:,end) = 0;
    
    Nodes = [X(:),Y(:),Z(:)];
    IsClosedDim = 1;
    [Geometry, nxp, nyp] = PointsToGeometry(Nodes,nx,ny,IsClosedDim);

    % Randomise node order.
    rng(1234);
    if (0)
        for i=1:size(Geometry.Patch,1)
            idx= mod((1:4) - 1 + 1,4)+1;
            %idx= mod((1:4) + round(4*rand), 4) + 1;
            Geometry.Patch(i,:) = Geometry.Patch(i,idx);
        end
        Geometry.P1 = Geometry.Patch(:,1);
        Geometry.P2 = Geometry.Patch(:,2);
        Geometry.P3 = Geometry.Patch(:,3);
        Geometry.P4 = Geometry.Patch(:,4);
    end
    
    %
    % Create the wake.
    %
    % Set up the nodes, we'll sort out the connected elements later.
    TENodesIdx = reshape(1:length(Nodes),[nx,ny]);
    
    NewNodes = Nodes(TENodesIdx(1,:),:) + Downstream * repmat([1,0,0],ny, 1);
    NewNodesIdx = (size(Nodes,1)+1):(size(Nodes,1)+size(NewNodes,1));
    Nodes = [Nodes;NewNodes];
    
    Wake.P2 = TENodesIdx(1,1:ny-1)';
    Wake.P1 = TENodesIdx(1,2:ny  )';
    Wake.P3 = NewNodesIdx(1:end-1)';
    Wake.P4 = NewNodesIdx(2:end)';
    Wake.Patch=[Wake.P1,Wake.P2,Wake.P3,Wake.P4];

    if (1)
        % Connect wake to surface elements by searching. This keeps it general
        % enough for use with other meshes.
        [Wake.UpperTE, Wake.LowerTE] = FindTE(Geometry, Wake, Nodes);
    else
        % Connect based on geometry structure
        Wake.UpperTE = sub2ind([nxp,nyp],   1*ones(1,nyp), 1:nyp)';
        Wake.LowerTE = sub2ind([nxp,nyp], nxp*ones(1,nyp), 1:nyp)';
    end
        
    %
    % Set up BCs and run solution.
    %
    Vfs = [cos(alpha), 0, sin(alpha);...
           0,          1, 0;...
           sin(alpha), 0, cos(alpha)] * [1;0;0];
       
    if (1)
        Solution = RunSolution(Geometry, Vfs, Nodes, Wake);
        save('WingSolution','Solution');
    else
        load('WingSolution');
    end
    
    Mu = reshape(Solution.mu,Geometry.n1,Geometry.n2);
    
    [Vn, Area] = ComputeVn(Nodes,Geometry);
    [Vl,Vm,~] = SurfaceSpeedsFromDoublets(Mu,Nodes,Geometry,Vfs);
    
    %
    % Calculate forces and Cp.
    %
    
    Vt = sqrt(Vl.^2+Vm.^2);
    q = 0.5*1.225*norm(Vfs)^2;
    %Vt=Vl;
    Cp = (1-Vt(:).^2);
    
    Forces = q * repmat(Area,1,3) .* Vn .* repmat(-Cp,1,3);
    Force = sum(Forces,1);
   
    % TEMP
    %Gradient = SurfaceGradientNodes(Geometry,Nodes,Solution.mu,Wake);
    Gradient = SurfaceGradientEdges( Geometry,Nodes,Solution.mu, Wake);
    Cp       = CpFromSurfaceGradient(Geometry,Nodes,Gradient,Vfs);
    %

    %
    % Show patch plot of circulation.
    %
    ConfigurationSurfacePlot(Geometry,Nodes,Wake,Solution.mu);
    
    %
    % Show patch plot of Cp.
    %
    ConfigurationSurfacePlot(Geometry,Nodes,Wake,-Cp);

    
    %
    % Create comparisons with Xfoil and Lifting-Line Theory.
    %
    CpArr=reshape(Cp,nxp,nyp);
    PC = (Nodes(Geometry.P1,:)+Nodes(Geometry.P2,:)+Nodes(Geometry.P3,:)+Nodes(Geometry.P4,:))/4;
    pcx = reshape(PC(:,1),[nxp,nyp]);
    
    XfoilDat = importdata('P:\ProjA\48731\48731-001\TECHNICAL\Unversioned code\2D_Aerofoil_Analysis\XFoil data\NACA2410.txt');
    
    % Integrate to find lift.
    % Xfoil - GUI reports 0.8443
    Xfoil.x  = XfoilDat.data(:,1);
    Xfoil.Cp = XfoilDat.data(:,2);
    Xfoil.Cl = trapz(Xfoil.x,Xfoil.Cp);

    PanelMethod.x  = pcx(:,10);
    PanelMethod.Cp = CpArr(:,10);
    PanelMethod.Cl = trapz(PanelMethod.x,PanelMethod.Cp);
    
    figure;
    
    axes('Position',[0.1,0.6,0.5,0.35]);
    hold on;
    plot(PanelMethod.x,-PanelMethod.Cp,'b-o');
    
    plot(Xfoil.x,-Xfoil.Cp,'r');
    legend('3D panel','Xfoil');
    grid on;
    xlabel('x');
    ylabel('-Cp');
    
    axes('Position',[0.1,0.13,0.5,0.35]);
    plot(NacaFoil.x,NacaFoil.z,'b-');
    grid on;
    xlabel('x');
    ylabel('y');
    
    axes('Position',[0.7,0.2,0.25,0.75]);
    bar([0.8443,Xfoil.Cl,PanelMethod.Cl]);
    set(gca,'XTickLabel',{'Xfoil GUI', 'Xfoil \int_{s}Cp', 'Panel Method'},...
            'XTickLabelRotation',45);
    ylabel('Cl');
    
    % Note that panel method will be slightly lower due to 3D effects. This 
    % difference should asymptote to zero as aspect ratio is increased.

    % LLT

    n=10;
    theta = (1:2*n)' * pi/(2*n+1);
    c_llt = ones(1,20);
    % Zero-lift and lift curve slope from Xfoil.
    alpha_eff = alpha + deg2rad(2.112);
    clalpha = 6.769;
    
    [Cl3D_Cl2D]= CalculateCoefficients(theta,c_llt,L,2*pi);
    cl = Cl3D_Cl2D.*clalpha*(repmat(alpha_eff,size(theta)));
    LLT.Position = L*cos(theta);
    LLT.Lift     = cl;
    
    Cl = trapz(PanelMethod.x,CpArr);
    Y = PC(sub2ind([nxp,nyp],ones(1,nyp),1:nyp),2);
    
    figure;hold on;grid on;
    plot(Y,Cl,'b-o')
    plot(LLT.Position,LLT.Lift,'r-x')
    set(gca,'YLim',[0,1]);
    
    legend('3D Panel','Lifting-Line');
    xlabel('Spanwise position');
    ylabel('Local Cl');
    

    [~,~,Pc] = ComputeVn(Nodes,Geometry);
    ElemIdx = Pc(:,2)>-0.1 & Pc(:,2)<0.1;

    PlotPanelCoordinateSystems(Geometry,Nodes,ElemIdx);
end

