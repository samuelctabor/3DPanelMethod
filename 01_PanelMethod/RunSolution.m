function Solution = RunSolution(Geometry, Vfs, Nodes, Wake)
    % Source/doublet panel method.
    nPanel = size(Geometry.P1,1);

    Pdoublet = zeros(nPanel,nPanel);
    Psource  = zeros(nPanel,nPanel);

    
    fprintf('Computing AICs . . .\n');
    fprintf('%i panels in model\n', nPanel);
    
    tic;
    
    vn  = ComputeVn(Nodes,Geometry);

    for iPanel = 1:nPanel
        if mod(iPanel/round(nPanel/10),1)==0
            fprintf('%i%%\n', round(100*iPanel/nPanel));
        end
            
        % For every panel, compute the influence of all panels on its control
        % point.
        PC = mean(Nodes([Geometry.P1(iPanel),Geometry.P2(iPanel),Geometry.P3(iPanel),Geometry.P4(iPanel)],:));
        
        
        Pinside = PC - 5*eps*vn(iPanel,:);
        
        rinside = repmat(Pinside,nPanel,1);

        [~,~,P_si,P_di] = ...
            SourceDoubletPanelInfluence4023(rinside, ...
                                            Nodes(Geometry.P1,:),...
                                            Nodes(Geometry.P2,:),...
                                            Nodes(Geometry.P3,:),...
                                            Nodes(Geometry.P4,:),...
                                            vn);

                                        
        if any(~isfinite(P_si)) || any(~isfinite(P_di))
            error('Stop');
        end

                                        
        Psource(iPanel,:)  = P_si;
        Pdoublet(iPanel,:) = P_di;
        
        % Add the influence of the wake
        if nargin>3
            wakeVn  = ComputeVn(Nodes,Wake);
            
            rinside = repmat(Pinside,size(Wake.P1,1),1);

            [~,~,~,P_di_w] = ...
            SourceDoubletPanelInfluence4023(rinside, ...
                                            Nodes(Wake.P1,:),...
                                            Nodes(Wake.P2,:),...
                                            Nodes(Wake.P3,:),...
                                            Nodes(Wake.P4,:),...
                                            wakeVn);
                                        
            Pdoublet(iPanel,Wake.UpperTE) = Pdoublet(iPanel,Wake.UpperTE) + ...
                                            P_di_w';
            Pdoublet(iPanel,Wake.LowerTE) = Pdoublet(iPanel,Wake.LowerTE) - ...
                                            P_di_w';

            Pwake(iPanel,Wake.UpperTE) = P_di_w';
            Pwake(iPanel,Wake.LowerTE) = P_di_w';

        end
    end
    toc;

    % We need zero peturbation potential on the inside of the surface.
    % So Psource*sigma + Pdoublet*mu = 0
    % Pdoublet*mu = - Psource*sigma
    % Pdoublet*mu = - Psource*-Vinf.n
    %          mu = Pdoublet\Psource*Vinf.n

    % The source panel provides the jump in velocity.
    % The doublet panel provides the jump in potential.

    Vbc = -dot(repmat(Vfs',nPanel,1),vn,2);
    
    sigma_dirichlet = Vbc/(4*pi);

    fprintf('Solving system . . .\n');

    mu_dirichlet = Pdoublet\(-Psource*sigma_dirichlet);

    Solution.sigma = sigma_dirichlet;
    Solution.mu = mu_dirichlet;
    Solution.Vbc = Vbc;
end