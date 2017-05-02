function TestBL_Drela()
    %TestBL Testa= an integral boundary layer formulation.
    %  Drela 
    
    % [1] G:\ST2_files\MATLAB\PotentialFlow\BoundaryLayers\Drela.pdf
    % [2] G:\ST2_files\MATLAB\PotentialFlow\BoundaryLayers\Veldman_interactionLaw.pdf
    
    close all;
    
    Fluid.rho = 1.225;
    Fluid.mu = 1.8e-5;
    
    x = linspace(0,10,10000)';
    dx = mean(diff(x))*ones(size(x));
    
    Ue = 1*ones(size(x));
    
    dUedx = zeros(size(x));
    
    H = 1.4*ones(size(x));
    theta = 1e-3*ones(size(x));
    
    figure('Name','Drela');
        
    h1 = subplot(3,1,1);
    hold on; grid on;
    title(h1,'H');
    
    h2 = subplot(3,1,2);
    hold on; grid on;
    title(h2,'\theta');
    
    h3 = subplot(3,1,3);
    hold on; grid on;
    title(h3,'Cf');
    

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
    end

    
    plot(h1,x,H,'b');
    plot(h2,x,theta,'r');
    plot(h3,x,Cf,'g');
    
    % Cf approximate formulae
    % http://www.cfd-online.com/Wiki/Skin_friction_coefficient
    ReX = Fluid.rho * x .* Ue / Fluid.mu;
    %plot(h3, x, 0.0576*ReX.^(-1/5),'k:');
    plot(h3, x, 0.0592*ReX.^(-1/5),'k:');
    
end

function B = CalculateB(H,theta)
    % d(States)/dx.
    % Drela's states are 
    % theta, b.l. momentum thickness
    % m, mass defect = ue * delta*
    % n or Ctau, laminar disturbance wave amplitude or turbulent shear stess coefficient
    
    % But I'm going to use:
    % theta
    % H
    % n or Ctau.
    % As it is more consistent with Veldman and with the formulation in [1]
    
    [~, dHstardH] = CalculateHstar(H);
    
    B = [1, 0; 0, theta*dHstardH];
    
end

function R = CalculateR(H, theta, dUedx, Ue, Fluid)
   % RHS vector.
    
   Cf = CalculateCf(theta,H,Ue,Fluid);
   
   Hk = CalculateHk(H);
   Hstar = CalculateHstar(Hk);
   Ctau = CalculateCtauEq(Hstar, Hk, H);
   CD = CalculateCD(Cf,H,Ctau);

   
   R = [Cf/2 - (2+H)*(theta/Ue)*dUedx;...
       2*CD - Hstar*Cf/2 - Hstar*(1-H)*theta/Ue * dUedx];
       
end

function [Hk] = CalculateHk(H)
    % Equation 9.
    Hk = H; % Low M
end

function [Hstar, dHstardH] = CalculateHstar(H)
    % Equation 10.
    Hstar = zeros(size(H));
    Hk = CalculateHk(H);
    Hstar(Hk<=4) = 1.515 + 0.076*( 4-Hk(Hk<=4)).^2./Hk(Hk<=4);
    Hstar(Hk>4)  = 1.515 + 0.040*(-4+Hk(Hk >4)).^2./Hk(Hk>4);
    
    dHstardH = (2/Hk)*(Hk-4) - (1/Hk^2)*(Hk-4)^2;
    dHstardH(Hk<=4) = 0.076*dHstardH(Hk<=4);
    dHstardH(Hk>4)  = 0.040*dHstardH(Hk>4);
end

function Hstarstar = CalculateHstarstar(H)
    %Equation 13 under low Me.
    Hstarstar = 0;
end

function cf = CalculateCf(theta, H, Ue, Fluid)
    % Skin friction coef Cf = f(theta,H,Ue)
    ReTheta = max(Fluid.rho * Ue .* theta / Fluid.mu, 12);
    Hk = CalculateHk(H);
    
    Fc = 1; % Low M
    cf = 0.3 * exp(-1.33*Hk) * (log10(ReTheta/Fc))^(-1.74-0.31*Hk) + ...
        0.00011*(tanh(4-Hk/0.0875) - 1);
   
    cf = cf/Fc;
end

function CtauEq = CalculateCtauEq(Hstar, Hk, H)
    
    Us = CalculateUs(Hstar,Hk,H);
    CtauEq = Hstar * (0.015/(1-Us)) * (Hk-1)^3/(Hk^2*H);
end


function Us = CalculateUs(Hstar,Hk,H)
    Us = (Hstar/2)*(1-(4/3)*(Hk-1)/H);
end

function CD = CalculateCD(Cf, H, Ctau)
    Hk = CalculateHk(H);
    Hstar = CalculateHstar(Hk);
    Us = CalculateUs(Hstar,Hk,H);
    CD = (Cf/2) * Us + Ctau * (1-Us);
end
