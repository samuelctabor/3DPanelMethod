function [ output_args ] = CompareH1Hstar()
    %UNTITLED3 Summary of this function goes here
    %   Detailed explanation goes here
    
    H= (1:10)';
    
    figure,hold on;
    plot(H,CalculateHstar(H),'b');
    plot(H,CalculateH1(H),'r');
    
    legend('H* (Drela)','H1 (Veldman)');
end

% From Drela.
function Hstar = CalculateHstar(H)
    Hstar = zeros(size(H));
    Hk = CalculateHk(H);
    Hstar(Hk<=4) = 1.515 + 0.076*( 4-Hk(Hk<=4)).^2./Hk(Hk<=4);
    Hstar(Hk>4)  = 1.515 + 0.040*(-4+Hk(Hk >4)).^2./Hk(Hk>4);
end

function [Hk] =CalculateHk(H)
    Hk = H; % Low M
end

% From Veldman.
function [H1, dH1dH] =CalculateH1(H)
    % Head's shape factor, H1 = f(H)
    dhtdH = 0*H;
    
    htCand = [H, 0.5*(H-2.732) + 2.732];
    [ht,Idx] = min(htCand,[],2);

    dhtdH(Idx==1) = 1;
    dhtdH(Idx==2) = 0.5;
    
    H1 = 0*H;
    dH1dH = 0*H;
    
    H1(H<=4) = 0.5*ht(H<=4).*(ht(H<=4)+2)./(ht(H<=4)-1);
    H1(H>4)  = 1.75 + 5.52273*ht(H>4) ./(ht(H>4)+5.818181);
    
    dH1dH(H<=4) = 0.5 * (ht(H<=4).^2 - 2*ht(H<=4) - 2)./(ht(H<=4)-1).^2   .* dhtdH(H<=4);
    dH1dH(H>4)  = 5.52273 * 5.81818 ./ (ht(H>4)+5.81818).^2  .* dhtdH(H>4);
end