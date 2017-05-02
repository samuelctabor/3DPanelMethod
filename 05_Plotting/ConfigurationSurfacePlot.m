function ConfigurationSurfacePlot(Geometry,Nodes,Wake,SurfaceField,ElemIdx,PlotWake)
    %ConfigurationSurfacePlot Patch plot of a quantity on the configuration
    %surface.
    
    if nargin<6 || isempty(PlotWake)
       PlotWake = false; 
    end
    if nargin<5 || isempty(ElemIdx)
        ElemIdx = true(size(Geometry.P1));
    end
    
    figure,hold on;
    x = Nodes(:,1); y=Nodes(:,2); z=Nodes(:,3);
    
    patch(x(Geometry.Patch(ElemIdx,:))', y(Geometry.Patch(ElemIdx,:))', z(Geometry.Patch(ElemIdx,:))', SurfaceField(ElemIdx));
    
    if PlotWake
        patch([Nodes(Wake.P1,1),Nodes(Wake.P2,1),Nodes(Wake.P3,1),Nodes(Wake.P4,1)]',...
              [Nodes(Wake.P1,2),Nodes(Wake.P2,2),Nodes(Wake.P3,2),Nodes(Wake.P4,2)]',...
              [Nodes(Wake.P1,3),Nodes(Wake.P2,3),Nodes(Wake.P3,3),Nodes(Wake.P4,3)]',...
              'r');
    end
    
    title('\mu');
    colorbar;
    colormap jet;
    
    %axis equal;
    set(gca,'XLim',[min(x(Geometry.P1))-1,max(x(Geometry.P1))+1]);
    grid on;
end

