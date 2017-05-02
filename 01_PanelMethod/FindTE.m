function [ upperTE, lowerTE ] = FindTE(Geometry, Wake, Nodes)
    %FindTE Discover the elements that a wake panel connects to.
    %   Search for elements that share two nodes with the wake panel.
    %   distinguish between 'upper TE' and 'lower TE' based on element
    %   normals.
    VnGeo  = ComputeVn(Nodes,Geometry);
    VnWake = ComputeVn(Nodes,Wake);
    
    upperTE = zeros(size(Wake.P1));
    lowerTE = zeros(size(Wake.P1));
    
    for iW = 1:length(Wake.P1)
        IsConnected = sum(ismember([Geometry.P1,Geometry.P2,Geometry.P3,Geometry.P4],...
                 [Wake.P1(iW), Wake.P2(iW), Wake.P3(iW), Wake.P4(iW)]),2)==2;
             
        if sum(IsConnected)~=2
            error('There should be two geometry panels connected to a wake panel');
        end
        
        IsConnectedIndex = find(IsConnected);
        
        % It's the upper TE if the element normals point in roughly the
        % same direction.
        VnWakeTmp = repmat(VnWake(iW,:),[length(IsConnectedIndex),1]);
        IsUpper = dot(VnWakeTmp,VnGeo(IsConnected,:),2)>0;
        
        upperTE(iW) = IsConnectedIndex(IsUpper);
        lowerTE(iW) = IsConnectedIndex(~IsUpper);
    end
end

