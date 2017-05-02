function [ h ] = TestSphere()
    %Test_Sphere Test function comparing 3D panel prediction to exact
    %solution for potential flow around a sphere.
    %   Detailed explanation goes here
    doplot = 1;
    
    nx = 20;
    nth= 20;

    rx=1;
    ry=1;

    bodyfunc2 = @(x) sqrt(1-x.^2);
    bodyfunc = @(x) bodyfunc2(x./rx);

    th = linspace(-pi/2,pi/2,nx)';

    xControl = rx*sin(th);

    R = ry*bodyfunc(xControl);

    Px = repmat(xControl,1,nth);

    th = linspace(0,2*pi,nth+1)-pi/nth;
    th = th(1:nth);
    Py = R*cos(th);
    Pz = R*sin(th);

    Nodes = [Px(:),Py(:),Pz(:)];
    Geometry = PointsToGeometry(Nodes,nx,nth,2);
    
    [Vn,~,PC] = ComputeVn(Nodes,Geometry);

    if doplot
        figure;
        x = Nodes(:,1); y = Nodes(:,2); z = Nodes(:,3);
        h=patch(x(Geometry.Patch)',...
                y(Geometry.Patch)',...
                z(Geometry.Patch)',...
                zeros(length(Geometry.P1),1)');
        set(h,'FaceAlpha',0.2);
        hold on;
        poutside = PC + 0.1*Vn;
        plot3(poutside(:,1),poutside(:,2),poutside(:,3),'ro');
        scale = 0.1;
        plot3([PC(:,1),PC(:,1)+scale*Vn(:,1)]',...
              [PC(:,2),PC(:,2)+scale*Vn(:,2)]',...
              [PC(:,3),PC(:,3)+scale*Vn(:,3)]','r-');
        xlabel('x');
        ylabel('y');
        zlabel('z');

        for iPanel = 1:length(PC)
            text(PC(iPanel,1),...
                 PC(iPanel,2),...
                 PC(iPanel,3),...
                 num2str(iPanel));
        end
    end

    
    Vfs = [1;0;0];
    Solution = RunSolution(Geometry, Vfs, Nodes);

    
    nPanel = size(PC,1);


    Vfs_temp = repmat(Vfs',nPanel,1);

    h=patch(x(Geometry.Patch)',...
            y(Geometry.Patch)',...
            z(Geometry.Patch)',...
            Solution.mu);
    set(gcf,'Renderer','zbuffer');
    title('mu');

    Mu = reshape(Solution.mu,Geometry.n1,Geometry.n2);
    [Vl,Vm,l] = SurfaceSpeedsFromDoublets(Mu,Nodes,Geometry,Vfs);
    m = cross(l,Vn,2);
    Vt = sqrt(Vl.^2+Vm.^2);
    %Vt=Vl;

    % Compare to the single doublet solution.
    r_eff = mean(rownormV(PC));
    mu_single = 2*pi*norm(Vfs)*r_eff^3 / (4*pi);
    C_theory = DoubletInfluence(PC,repmat(-Vfs',nPanel,1));


    % Plot 4*pi*grad(mu)
    px = dot(Vfs_temp,PC,2);

    th = acosd(px);

    % Plot |Vinf|
    VdotL.FarField         = abs(dot(Vfs_temp,l,2) + dot(Vfs_temp,m,2));
    VdotL.PointDoublet     = rownormV(Vfs_temp+C_theory*mu_single);
    VdotL.Theory           = 3/2 * norm(Vfs) * sind(th);
    VdotL.FFPlusGradMu     = Vt(:);
    
    figure;
    hold on;
    plot(th,VdotL.FarField,'m--',...
         th,VdotL.FFPlusGradMu,'g-o',...
         th,VdotL.Theory,'r',...
         th,VdotL.PointDoublet,'r^');
    
    title('Vt');
    legend('|v_{\infty}.l|','|v_{\infty}.l| + 4\pi\nabla\mu','Theory','Point doublet');
end

