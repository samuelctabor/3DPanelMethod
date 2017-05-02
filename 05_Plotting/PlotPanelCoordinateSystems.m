function [ output_args ] = PlotPanelCoordinateSystems(Geometry, Nodes, ElemIdx, scale, showNumbers)
    %UNTITLED2 Summary of this function goes here
    %   Detailed explanation goes here
        % Plot surface normals
    
    if nargin<3 || isempty(ElemIdx)
        ElemIdx = true(Geometry.P1);
    end
    
    if nargin<4 || isempty(scale)
        scale = 0.01;
    end
    
    if nargin<5 || isempty(showNumbers)
        showNumbers = false;
    end
    
    x = Nodes(:,1); y=Nodes(:,2); z=Nodes(:,3);
    
    [Vn,~,Pc] = ComputeVn(Nodes,Geometry);
    Vl = (Nodes(Geometry.P1,:) + Nodes(Geometry.P2,:))/2 - Pc;
    Vl = Vl./rownorm(Vl);
    Vm = cross(Vn, Vl);
    
    
    figure,hold on;
    patch(x(Geometry.Patch(ElemIdx,:))', y(Geometry.Patch(ElemIdx,:))', z(Geometry.Patch(ElemIdx,:))',[0.5,0.5,0.5],'FaceAlpha',0.5);
    
    if showNumbers
        labels = arrayfun(@(i) sprintf(' %i',i), 1:length(Geometry.P1), 'UniformOutput',false);
        text(mean(x(Geometry.Patch),2),mean(y(Geometry.Patch),2),mean(z(Geometry.Patch),2),labels);
    end
    
    plot3([Pc(ElemIdx,1),Pc(ElemIdx,1) + Vn(ElemIdx,1)*scale]',...
          [Pc(ElemIdx,2),Pc(ElemIdx,2) + Vn(ElemIdx,2)*scale]',...
          [Pc(ElemIdx,3),Pc(ElemIdx,3) + Vn(ElemIdx,3)*scale]','r-o');
      
    plot3([Pc(ElemIdx,1),Pc(ElemIdx,1) + Vl(ElemIdx,1)*scale]',...
          [Pc(ElemIdx,2),Pc(ElemIdx,2) + Vl(ElemIdx,2)*scale]',...
          [Pc(ElemIdx,3),Pc(ElemIdx,3) + Vl(ElemIdx,3)*scale]','b-o');
      
    plot3([Pc(ElemIdx,1),Pc(ElemIdx,1) + Vm(ElemIdx,1)*scale]',...
          [Pc(ElemIdx,2),Pc(ElemIdx,2) + Vm(ElemIdx,2)*scale]',...
          [Pc(ElemIdx,3),Pc(ElemIdx,3) + Vm(ElemIdx,3)*scale]','g-o');
    axis equal;
    
end

