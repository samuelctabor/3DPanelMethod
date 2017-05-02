function [ C_s, C_d, Phi_s, Phi_d ] = SourceDoubletPanelInfluence4023( Pp, p1, p2, p3, p4, n)

%SourceInfluence Effect of source singularity
    
    % rp is vector to point of iterest (global frame).
    % nv is normal vector (global frame).
    % A is area.
    % r1-4 are corner positions (local frame).
    
    % Panel needs to be flat and corner points are x1-4 and y1-4.
    
    % p vectors are panel corners relative to point of interest.
    
 
    nP = size(Pp,1);
    
    % Find l. 
    % Centre position.
    pc = (p1+p2+p3+p4)/4;

    rmid = (p1 + p2)/2 - pc;
    % Move onto plane through rc normal to n 
    %(remove component parallel to normal)
    rmid = rmid - n.*repmat(dot(rmid,n,2),1,3);
    l = rmid./rownorm(rmid);

    m = cross(n,l,2);

    PC = cat(3,p1,p2,p3,p4,p1);
    
    PJK = Pp-pc;
    pjk = rownormV(PJK);
    
    C_s = zeros(nP,3);
    Phi_s = zeros(nP,1);
    
    C_d = zeros(nP,3);
    Phi_d = zeros(nP,1);
    
    Area = rownormV(cross(p3-p1,p4-p2,2))/2;
    %IsClose = pjk < 5*sqrt(Area);
    IsClose = true(size(pjk));
    
    [C_s(IsClose,:), C_d(IsClose,:), Phi_s(IsClose,:), Phi_d(IsClose,:)] = ...
        NearField(PC(IsClose,:,:), Pp(IsClose,:), pc(IsClose,:),l(IsClose,:), m(IsClose,:), n(IsClose,:));
    
%     
%     [C_s(~IsClose,:), C_d(~IsClose,:), Phi_s(~IsClose,:), Phi_d(~IsClose,:)] = ...
%         FarField(Pp(~IsClose,:), pc(~IsClose,:), n(~IsClose,:),Area(~IsClose));
    
end

function [ C_s, C_d, Phi_s, Phi_d ] = NearField(PC, Pp, pc, l, m, n)
   
    PJK = Pp-pc;
    PN = dot(PJK,n,2);
    
    nPanel = size(PC,1);
    
    C_s   = zeros(nPanel,3);
    C_d   = zeros(nPanel,3);
    phi_s = zeros(nPanel,1);
    phi_d = zeros(nPanel,1);
    
    for i=1:4
        a = Pp-PC(:,:,i);
        b = Pp-PC(:,:,i+1);
        s = a - b;  % See fig 14

        A = rownormV(a);
        B = rownormV(b);
        S = rownormV(s);                
        SM = dot(s,m,2);
        SL = dot(s,l,2);
        AM = dot(a,m,2);
        AL = dot(a,l,2);
        Al = AM.*SL - AL.*SM;
        PA = PN.^2 .* SL + Al.*AM;
        PB = PA - Al .* SM;
        
        GL = (1./S) .* log(abs((A+B+S)./(A+B-S)));
        
        
        RNUM = SM.*PN .* (B.*PA-A.*PB);
        DNOM = PA.*PB + PN.^2.*A.*B.*SM.^2;
        Cjk = atan2(RNUM,DNOM);
        
        isClose = abs(PN)<10*eps;
        
        % If the point really is on the panel plane, let it be slightly above.
        panelSide = sign(PN);
        panelSide(panelSide==0) = 1;
        
        h = cross(a,s);
        % By my calculation this should be negative so that stripside is
        % positive when Pp lies on the right of s.
        stripSide = -sign(dot(h,n,2));
        stripSide(stripSide==0) = 1;
        
        direction = panelSide.*stripSide;

        
        % Within the strip
        Cjk(isClose & DNOM<-eps )     = pi   * direction(isClose & DNOM<-eps);
        % On the edge of the strip
        Cjk(isClose & abs(DNOM)<=eps) = pi/2 * direction(isClose & abs(DNOM)<=eps);
        % Outside the strip
        Cjk(isClose & DNOM>eps )      = 0    * direction(isClose & DNOM>eps);
        
        % Short side correction.
        idx = find(abs(S)<eps);
        
        % Velocity vector contribution
        %Ci_d = cross(a,b,2).*repmat((A+B)./(A.*B.*(A.*B + dot(a,b,2))),1,3);
        
        %onLine = abs(A.*B.*(A.*B + dot(a,b,2)))<10*eps;
        
        %Ci_d(onLine,:)=0;
        %Ci_d(idx,:) = zeros(numel(idx),3);

        %C_d = C_d + Ci_d;

        phi_d(:,i) = Cjk;
        phi_d(idx,i) = zeros(numel(idx),1);
        
        
        %Ci_s = (repmat(SM.*GL,1,3).*l - repmat(SL.*GL,1,3).*m) + repmat(Cjk,1,3).*n;
        %Ci_s(idx,:) = zeros(numel(idx),3);
        
        %C_s = C_s + Ci_s;
        
        phi_s(:,i) = (Al.*GL - PN.*Cjk);
        phi_s(idx,i) = zeros(numel(idx),1);
    end

    Phi_s = sum(phi_s,2);
    Phi_d = sum(phi_d,2);
end

function [ C_s, C_d, Phi_s, Phi_d ] = FarField(Pp, pc, n, Area)
    
    PJK = Pp-pc;
    PN = dot(PJK,n,2);
    
    C_d = repmat(Area,1,3).*(3*repmat(PN,1,3).*PJK - rownorm(PJK).^2.*n)./(rownorm(PJK).^5);
    Phi_d = Area.*PN./rownormV(PJK).^3;
    
    C_s = repmat(Area,1,3).*PJK./(rownorm(PJK).^3);
    Phi_s = Area./rownormV(PJK);
end
    
