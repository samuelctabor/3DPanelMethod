function [ output_args ] = TestBL_VeldmanTurbulent()
    %TestBL Testa= an integral boundary layer formulation.
    %   This is based on "A simple interactio law for viscous-inviscid
    %   interaction", Arthur E. P. Veldman
    % See also "An interactive boundary layer modelling methodology for
    % aerodynamic flows", L. Smith
    
    
    % 
    %close all;
    
    Fluid.rho = 1.225;
    Fluid.mu = 1.8e-5;
    
    n=1000;
    x = (1:n)'/(n);
    
    if (1)
        % Leading-edge clustering.
        x = 1-cos(x*pi/2);
    end
    
    dx = mean(diff(x))*ones(size(x));
    
    Ue = 1*ones(size(x));
    
    dUedx = zeros(size(x)) - 0.0;
    
    if (0)
    % We could initiate the solution using dstar0 = 1.7208 sqrt(nu*x0/Ue),
        % This doesn't work, probably because it gives a shape factor of
        % 2.59 (laminar)
        delta1 = 1.7208*sqrt(Fluid.mu*x(1)/(Fluid.rho*Ue(1)));
        % theta0 = 0.664*sqrt(nu*x0/Ue)
        theta(1) = 0.664*sqrt(Fluid.mu*x(1)/Fluid.rho*Ue(1));
        H(1) = delta1/theta(1);
    else
        %H = 2.59*ones(size(x));
        H = 1.4*ones(size(x));
        %theta = 1e-3*ones(size(x));
        theta = 1e-3*ones(size(x));
    end
    Cf = zeros(size(x));
    
    figure('Name','Veldman');

    
    h1 = subplot(4,1,1);
    hold on; grid on;
    title(h1,'H');
    
    h2 = subplot(4,1,2);
    hold on; grid on;
    title(h2,'\theta');
    
    h3 = subplot(4,1,3);
    hold on; grid on;
    title(h3,'C_f');
    
    h4 = subplot(4,1,4);
    hold on; grid on;
    title(h4,'\delta^*');
    
    xSep = NaN;
    
    for i=2:length(x)

        % What station to use for B and R terms.
        % i-1 for the upstream station.
        % Could eventually use the same station with iteration.
        Idx = i-1; 
        B = CalculateB(H(Idx), theta(Idx));
        R = CalculateR(H(Idx), theta(Idx), dUedx(i), Ue(i), Fluid);

        % Phi [theta; H]
        dPhidx = B\R;

        theta(i) = theta(i-1) + dPhidx(1)*dx(i);  
        H(i)     = H(i-1)     + dPhidx(2)*dx(i);
        Cf(i) = CalculateCf(theta(i), H(i), Ue(i), Fluid);
        
        if Cf(i)<=5e-4 && isnan(xSep)
            fprintf('Separation at %f!\n', x(i));
            xSep = x(i);
            %break;
        end
    end
    
    deltaStar = theta.*H;
    delta = 8*deltaStar; % Assuming a 1/7 profile
    
    plot(h1,x,H,'b.-');
    plot(h1,[xSep,xSep],get(h1,'YLim'),'r:');
    
    plot(h2,x,theta,'r');
    plot(h2,[xSep,xSep],get(h2,'YLim'),'r:');
    
    plot(h3,x,Cf,'g');
    plot(h3,[xSep,xSep],get(h3,'YLim'),'r:');
    
    plot(h4,x,deltaStar,'g');
    plot(h4,[xSep,xSep],get(h4,'YLim'),'r:');
    
    
    % Cf approximate formulae
    % Incropera 7.35
    ReX = Fluid.rho * x .* Ue / Fluid.mu;
    plot(h3, x, 0.0592*ReX.^(-1/5),'k:'); 
    
    % delta approximate forula
    % Incropera 7.36
    delta = 0.37*x.*ReX.^(-1/5);
    %Assuming a 1/7 profile http://personalpages.manchester.ac.uk/staff/david.d.apsley/lectures/turbbl/integral.pdf
    deltaStar2 = delta/8;  
    plot(h4, x, deltaStar2,'k:');        
end

function B=CalculateB(H,theta)
    % d(States)/dx.
    [H1, dH1dH] = CalculateH1(H);
    B = [1, 0; H1, theta*dH1dH];
end

function R = CalculateR(H, theta, dUedx, Ue, Fluid)
   % RHS vector.
   
   H1 = CalculateH1(H);
   
   % Entrainment coefficient, Ce = f(H1)
   Ce = 0.0306*(H1 - 3.0).^(-0.6169);
   
   cf = CalculateCf(theta, H, Ue, Fluid);
   
   R(1,1) =  0.5*cf - theta./Ue * (2+H) .* dUedx;
   R(2,1) = Ce - H1.*theta.*dUedx./Ue;
end

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

function cf = CalculateCf(theta, H, Ue, Fluid)
    % Skin friction coef Cf = f(theta,H,Ue)
    ReTheta = max(Fluid.rho * Ue .* theta / Fluid.mu, 12);
    
    cf0 = 0.01013./(log10(ReTheta) - 1.02) - 0.00075;
    h0 = 1.0 - 6.55*sqrt(0.5*cf0);

    cf = cf0 .* (0.9./(H.*h0 - 0.4) - 0.5);
end