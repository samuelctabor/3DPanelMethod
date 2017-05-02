function [ vl,vm, l, m ] = SurfaceSpeedsFromDoublets( mu,Nodes,Elements,Vfs)
%SurfaceSpeedsFromDoublets 
%   The surface speed is d mu / d s.
    p1 = Nodes(Elements.P1,:);
    p2 = Nodes(Elements.P2,:);
    p3 = Nodes(Elements.P3,:);
    p4 = Nodes(Elements.P4,:);

    % Find the l and m distances between panel centres.
    %nl = size(mu,1);
    %nm = size(mu,2);
    nPanel = numel(mu);
    
    n = cross(p3-p1,...
              p4-p2);
    n = n./rownorm(n);
    
    % Should really move this stuff elsewhere!
        rc = (p1+p2+p3+p4)/4;

        r1 = p1-rc;
        r2 = p2-rc;
        
        midp = (r1+r2)/2;
        % Move onto plane through rc normal to n 
        %(remove component parallel to normal)
        midp = midp - n.*repmat(dot(midp,n,2),1,3);
        l = midp./rownorm(midp);

        m = cross(n,l,2);
        
    %
    
    PC = rc;


    nl = size(mu,1);
    nm = size(mu,2);
    
    delP = zeros(nl,nm);
    delQ = zeros(nl,nm);
    
    
    XC = reshape(PC(:,1),nl,nm);
    YC = reshape(PC(:,2),nl,nm);
    ZC = reshape(PC(:,3),nl,nm);
    
    % l direction.
    dXC = XC(1:nl-1,:)-XC(2:nl,:);
    dYC = YC(1:nl-1,:)-YC(2:nl,:);
    dZC = ZC(1:nl-1,:)-ZC(2:nl,:);
    
    h = sqrt(dXC.^2+dYC.^2+dZC.^2);
    
    h1= h(1:nl-2,:);
    h2= h(2:nl-1,  :);
   
    delP(2:nl-1,:) = - h2./(h1.*(h1+h2)) .* mu(1:nl-2,:) - ...
                     (h1-h2)./(h1.*h2) .* mu(2:nl-1,:) + ...
                      h1./(h2.*(h1+h2)).* mu(3:nl,:);
    
    h1= h(1,:);
    h2= h(2,:);
    
    delP(1,:) = -  (2*h1+h2)./(h1.*(h1+h2)).* mu(1,:) + ...
                 (h1+h2)./(h1.*h2)       .* mu(2,:) - ...
                 h1./((h1+h2).*h2)       .* mu(3,:);
                 
    h1= h(nl-1,:);
    h2= h(nl-2,:);
    
    delP(nl,:) =  + (2*h1+h2)./(h1.*(h1+h2)).* mu(nl,:) - ...
                 (h1+h2)./(h1.*h2)        .* mu(nl-1,:) + ...
                 h1./((h1+h2).*h2)        .* mu(nl-2,:);
    
    
    % m direction.
    dXC = XC(:,1:nm-1)-XC(:,2:nm);
    dYC = YC(:,1:nm-1)-YC(:,2:nm);
    dZC = ZC(:,1:nm-1)-ZC(:,2:nm);
    
    h = sqrt(dXC.^2+dYC.^2+dZC.^2);
    
    h1= h(:,1:nm-2);
    h2= h(:,2:nm-1);
   
    delQ(:,2:nm-1) = - h2./(h1.*(h1+h2)) .* mu(:,1:nm-2) - ...
                     (h1-h2)./(h1.*h2) .* mu(:,2:nm-1) + ...
                      h1./(h2.*(h1+h2)).* mu(:,3:nm);
    
    h1= h(:,1);
    h2= h(:,2);
    
    delQ(:,1) = -  (2*h1+h2)./(h1.*(h1+h2)).* mu(:,1) + ...
                 (h1+h2)./(h1.*h2)       .* mu(:,2) - ...
                 h1./((h1+h2).*h2)       .* mu(:,3);
                 
    h1= h(:,nm-1);
    h2= h(:,nm-2);
    
    delQ(:,nm) =   (2*h1+h2)./(h1.*(h1+h2)).* mu(:,nm) - ...
                 (h1+h2)./(h1.*h2)        .* mu(:,nm-1) + ...
                 h1./((h1+h2).*h2)        .* mu(:,nm-2);
    
    % velocities
    % Side 2 
    ps2 = (p3+p4)./2 - rc;
    
    %vl = -4*pi*(rownormV(ps2).*delP(:) + dot(ps2,m,2).*delQ(:))./dot(ps2,l,2);
    vl = -4*pi*delP;
    vl = reshape(vl,nl,nm);
    vm = -4*pi*delQ;
    
    Vfs_temp = repmat(Vfs',nPanel,1);
    vl = vl + reshape(dot(Vfs_temp,l,2),nl,nm);
    vm = vm + reshape(dot(Vfs_temp,m,2),nl,nm);
end

