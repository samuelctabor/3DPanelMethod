% Test SourcePanel.
close all;


% P1 = [ 1,-1,0]/2;
% P2 = [ 1, 1,0]/2;
% P3 = [-1, 1,0]/2;
% P4 = [-1,-1,0]/2;
P1 = [-1,  1,0]/2;
P2 = [-1, -1,0]/2;
P3 = [ 1, -1,0]/2;
P4 = [ 1,  1,0]/2;
   

vn = cross(P3-P1,P4-P2);
Area = rownormV(vn)/2;
vn = vn./rownorm(vn);

PC = (P1+P2+P3+P4)/4;

Scan = [-10:0.01:10];

%Scan = [0:0.01:10];

% Z transect.
Origin = [0,0,0];
ScanV = [0,0,1];
VelV = [0,0,1];

% X transect
%Origin = [0,0,1];
%ScanV = [1,0,0];

% XY transect
%Orgin = [0,0,1];
%ScanV = [0,1/sqrt(2),1/sqrt(2)];


X = Origin(1) + Scan*ScanV(1);
Y = Origin(2) + Scan*ScanV(2);
Z = Origin(3) + Scan*ScanV(3);

% Temp: single panel and point.
Type = 0;
switch Type
    case 0
        X = 0;
        Y = 0;
        Z = 1;
    case 1
        l = -[cosd(-20),0,sind(-20)];
        P = [0.5,0,0] - l;
        X = P(1);
        Y = P(2);
        Z = P(3);
    case 2
        X = 0;
        Y = 0;
        Z = -1e-20;
    case 3
        X = 0;
        Y = 0;
        Z = 1e-20;
    case 4 
        X = 0;
        Y = 0;
        Z = -1;
end
% 
R = Scan;
n = numel(X);
 
[Cs,Cd,Ps,Pd] = SourceDoubletPanelInfluence4023([X(:),Y(:),Z(:)],repmat(P1,n,1),repmat(P2,n,1),repmat(P3,n,1),repmat(P4,n,1),repmat(vn,n,1));
% dot(Cs,l,2)
% dot(Cd,l,2)q
[Cd_s,Pd_s] = DoubletInfluenceNew([X(:),Y(:),Z(:)],repmat(vn,n,1),Area);
[Cs_s,Ps_s] = SourceInfluence([X(:),Y(:),Z(:)],Area);

X = [PC(1)-5*eps*vn(1); PC(1)+5*eps*vn(1)];
Y = [PC(2)-5*eps*vn(2); PC(2)+5*eps*vn(2)];
Z = [PC(3)-5*eps*vn(3); PC(3)+5*eps*vn(3)];
n = numel(X);

Rc = Z(:);
[Cs_c,Cd_c,Ps_c,Pd_c] = SourceDoubletPanelInfluence4023([X(:),Y(:),Z(:)],repmat(P1,n,1),repmat(P2,n,1),repmat(P3,n,1),repmat(P4,n,1),repmat(vn,n,1));

VelVec = repmat(VelV,numel(R),1);

figure,hold on;
plot(R,Ps,'b',R,Ps_s,'r');
plot(Rc,Ps_c,'b:');
plot([-10,0 ],[Ps_c,Ps_c],'b:');
plot([ 0 ,10],[Ps_c,Ps_c],'b:');
title('Source Phi');xlabel('z');ylabel('\Phi');
set(gca,'YLim',[-10,0]);

figure,hold on;
title('Source Vel');xlabel('R');ylabel('V');set(gca,'YLim',[-10,10]);
plot(R,dot(Cs,VelVec,2),'b',R,dot(Cs_s,VelVec,2),'r');
plot([-10,0 ],[Cs_c(1,3),Cs_c(1,3)],'b:');
plot([ 0 ,10],[Cs_c(2,3),Cs_c(2,3)],'b:');

figure,hold on;
title('Doublet Phi');xlabel('R');ylabel('\Phi');set(gca,'YLim',[-10,10]);
plot(R,Pd,'b',R,Pd_s,'r');
plot([-10,0 ],[Pd_c(1),Pd_c(1)],'b:');
plot([ 0 ,10],[Pd_c(2),Pd_c(2)],'b:');

figure,hold on;
title('Doublet Vel');xlabel('R');ylabel('V');set(gca,'YLim',[-10,10]);
plot(R,dot(Cd,VelVec,2),'b',R,dot(Cd_s,VelVec,2),'r');
plot([-10,0 ],[Cd_c(1,3),Cd_c(1,3)],'b:');
plot([ 0 ,10],[Cd_c(2,3),Cd_c(2,3)],'b:');

