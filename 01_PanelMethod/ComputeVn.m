function [vn, area, Pc] = ComputeVn(Nodes,Elems)
    p1 = Nodes(Elems.P1,:);
    p2 = Nodes(Elems.P2,:);
    p3 = Nodes(Elems.P3,:);
    p4 = Nodes(Elems.P4,:);
    
    vn  = cross(p3-p1,...
                p4-p2,2);
    area = rownormV(vn)/2;
    vn = vn./rownorm(vn);
    
    Pc = (p1+p2+p3+p4)/4;
end