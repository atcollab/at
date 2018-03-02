function varargout = hassinskyfit(sigma0, q, R, L, frf, h, V, phis, varargin)

% [z lambda sigma mu] = HASSINSKYFIT(SIGMA0, R, L)
% This function will generate a charge distribution for a given SIGMA0 
% (zero current bunch length), q (bunch total charge in Coloumbs and can be
% a vector), using a generalised impedance model with R (OHms) and L
% (Henry) according to the Haissinski equation.
%
% SIGMA0 (meters, scalar)  : zero current bunch length
% q      (coloumbs, vector): total charge in bunch
% R      (Ohms, scalar)    : ring real effective resistance
% L      (Henry, scalar)   : ring real effecitve impedance
% frf    (Hz, scalar)      : RF frequency
% h      (integer, scalar) : harmonic number of the ring
% V      (Volts, scalar)   : cavity voltage
% phis   (radians, scalar) : synchronous phase
%
% Output:
%
% z      (meters, vector)     : longitudinal position
% lambda (normalised, vector) : charge distribution along z.
% sigma  (meters, vector)     : bunch length for each q using fitgaussian.
% mu     (meters, vector)     : mu of gaussian fit from fitgaussian.
%
% WARNING: The gaussian fit uses the following command:
%                  fitgaussian(I2,'plot',0,'bg_grad',0,'DC',0);
%          and will thus be a sigma in "points/pixels", no background
%          gradient and zero DC. Will fit the area, sigma, assymetry factor
%          and the centre of gravity for the gaussian.
%
% See also FITGAUSSIAN
%
% 02/03/2010 Eugene

DEBUG = 0;
TOL = 1e-10;
MAXIT = 500;

if DEBUG
%     sigma0 = 20.3e-12*3e8; %natural bunch length in m (145.5/5.6)*e-12*3e8
%     R=1600; %resistance (in Ohm)
%     L=25e-9; %Inductance (in nanoHenry)
%     frf = 499.672e6; %RF frequency
%     V=3000000; %RF voltage
%     q = ([1:10]/1000)/1.388e6;
%     h = 360;
%     phis = 2.83964;
end

c = PhysConstant.speed_of_light_in_vacuum.value;
dz = sigma0/50;
z = [-sigma0*6:dz:sigma0*6];

% initialise
sigmas = zeros(size(q));
mus = zeros(size(q));
lambda = zeros(length(z),length(q));

for i=1:length(q)
    mu = 0;
    K = -c^2*q(i)/(2*pi*frf*V*cos(phis)*(sigma0^2));
     
    if 0
        % test solution to Hassinski equation based on Oide and Yokoya.
        % The method has a different method of finding the solution to the
        % Haissinski equation where you integrate from "infinity" at the
        % head the bunch towards the tail. This avoids the a fixed point
        % type algorithm to find the equilibrium solution. This bit has yet
        % to be tested and debugged fully
        % Eugene 2/3/2010
        warning('This algorithm not tested. Use with caution.');
        I2 = zeros(size(z));
        A = 1/(sigma0*sqrt(2*pi)); % Initial normalisation factor
        sumint = 0;
        while abs(sumint-1) > 0.0001
            for j=(length(z)-1):-1:1
                % Assume that the solution for I2 is zero for suffidciently
                % large z.
                intR = R*sum(I2(j:end))*dz;
                intL = c*L*sum(gradient(I2(j:end),dz))*dz;
                I2(j) = A*(exp(-(z(j)-mu).^2./(2*sigma0^2))).*(exp((-K.*(intR + intL))));
            end
            sumint = (sum(I2)*dz);
            A = A/sumint;
        end
        converge = 1;
    else
        % Fixted point algorithm based on Rohan's old code that has been
        % put into functinal form. 
        if i==1
            % initial distribution, subsequent starting distributions set at
            % the bottom of the loop.
            %         kl = R/2/sqrt(pi)/sigma0;
            mu = 0; %c*kl*q(i)/(2*pi*frf*V*cos(phis));
            I = (1/sqrt(2*pi)/sigma0)*exp(-(z-mu).^2/(2*sigma0^2));
            I = I./(sum(I)*dz); % Normalise
            Ip = -((z-mu)/sigma0^2).*I;
        else
            I = I_stable;
            Ip = Ip_stable;
        end
        
        converge = 0;
        itno = 1;
        while ~converge
            intR = R*cumsum(I)*dz;
            intL = c*L*cumsum(Ip)*dz;
            I2 = (1/(sigma0*sqrt(2*pi)))*(exp(-(z-mu).^2./(2*sigma0^2))).*(exp((-K.*(intR + intL))));
            
            %Renormalise the function
            I2 = I2./(sum(I2)*dz);
            if DEBUG
                fprintf('Iteration: %03d, STD: %e\n',itno,std(I2(:)-I(:)));
            end
            stderror = std(I2(:)-I(:));
            if stderror < TOL; %convegence condition
                converge=1;
                % plot(I2,'r'); %plot final distribution in red
            else
                I = I2;
                Ip = gradient(I2,dz);
            end
            itno = itno + 1;
            
            %         if DEBUG && rem(itno,10) == 0
            %             figure(98);
            %             plot(I2);
            %         end
            
            if itno > MAXIT
                fprintf('WARNING: Did not converge during fit of Hassinsky''s equation!!\n');
                fprintf('q: %g, std error: %g\n',q(i),stderror);
                figure(99);
                plot(z,I);
                break;
            end
        end
    end
    
    lambda(:,i) = I2;
    
    if nargout > 2
        if converge
            % units of meters
            temp = fitgaussian(I2,'plot',1,'bg_grad',0,'DC',0,'scale',dz);
            sigmas(i) = temp.sigma;
            mus(i) = temp.mu - z(end);

            % Set I and Ip for the next sequence and only keep "stable"
            % solutions that converge. Otherwise small instabilities will
            % build up too quickly for subsequent solutions with higher
            % currents, q.
            %         kl = R/2/sqrt(pi)/sigmas(i-1);
            mu = 0; %c*kl*q(i)/(2*pi*frf*V*cos(phis));
            I_stable = I2;
            Ip_stable = gradient(I2,dz);
        else
            % Will try subset (plot envelope). There are "solutions" where
            % a sinusoidal pattern is observed at the head of the bunch.
            % This part of the code attempts to plot an envelope gaussian
            % function by looking for the peaks of the oscilating function
            % and only fitting to those points.
            ii = find(gradient((lambda(:,i))) > 0); 
            ii = ii(abs(diff(ii)) > 1);
            if ~isempty(ii) && length(ii) > 1
                 % it usually can't find where the gaussian tail of the
                 % bunch starts. So assume the distance from the last peak
                 % to the start of the tail of the gaussian bunch is the
                 % same as the distance between the head of the gaussian
                 % part of the bunch and the first peak.
                 ii(end+1) = ii(end) + (ii(2)-ii(1))*2; % Generate last point as it can't seem to find the last turn around point.
                 %  datai = [1:ii(1), ii(2:end-1)', ii(end):length(z)];
                 datai = [1:ii(1), ii(end):length(z)];
            else
                % sometimes convertance isn't met however the nongaussian
                % element is not strong so a standard gaussian fit is
                % sufficient.
                datai = [1:length(z)]';
            end
            
            z_subset = z(datai);
            I2_subset = lambda(datai,i);
            
            temp = fitgaussian(I2_subset,'plot',0,'bg_grad',0,'DC',0,'scale',1,'x',z_subset);
            % pause;
            % must have reached iteration limit without convergence
            % therefore cannot trust fitdata
            sigmas(i) = temp.sigma;
            mus(i) = temp.mu - z(end);
        end
    else
        I_stable = I2;
        Ip_stable = gradient(I2,dz);
    end
end

varargout{1} = z;
varargout{2} = lambda;
varargout{3} = sigmas;
varargout{4} = mus;


return



% Loop over currents
current = [];
sigmas = [];
I2 = [];
initstep = 0.5;
itno = 0;
for loop = 1:endloop
if itno >= 600
    initstep = initstep/2;
end
% in mA
if loop == 1
    current(loop) = 0.001;
else
    current(loop) = current(loop-1) + initstep;
end
%current(loop) = initstep*(loop-1) +0.001; 
%I = zeros(1001,1); %initial charge distribution
%I2 =zeros(1001,1); %final charge distribution 
%Ip = zeros(1001,1); %deriviative of charge distribution
c=3e8;
Frf = 499.672e6; %RF frequency
V=3000000; %RF voltage
Q=(current(end)/1000)/1.3880e6; %convert current to charge 

%I2p/I2 = -z/sigma^2 + (c^2*Q*(R*I + c*L*Ip))/(2*pi*Frf*V*cosd(158)*sigma^2)

stop = 0;
itno=0;
%iterate!
SIGMAx = sigma;
K = -c^2*Q/(2*pi*Frf*V*cosd(158)*(SIGMAx^2));
%set the initial distribution as a gaussian
if loop==1
    Ax = 1;
    MUx = 0; %centre it
    dz = SIGMAx/50;
    z = [-SIGMAx*6:dz:SIGMAx*6];
    I = 1 * (1/(SIGMAx*sqrt(2*pi)))*exp(-(((z-MUx).^2)/(2*SIGMAx^2)));
    % Normalise the initial gaussian
    I = I./(sum(I)*dz);
    % figure(100); hold on; plot(I);  %plot initial distribution
    % make derivative disribution 'Ip'
    Ip = -((z-MUx)/SIGMAx^2).*I;
else
    I = gaussianY;
    Ip = gradient(I,dz);
end

while stop==0;
    itno=itno+1;
    
    for a=1:length(z) %do the I(z'-z)dz itergrations first
        intR = sum(R*I(1:a))*dz;
        intL = sum(c*L*Ip(1:a))*dz;
        %the Hassinski equation:
        I2(a) = (1/(SIGMAx*sqrt(2*pi)))*(exp(-(z(a)-MUx)^2/(2*SIGMAx^2)))*(exp((-K*(intR + intL))));
    end
    %Renormalise the function 
    I2 = I2./(sum(I2)*dz);
 
    %figure(300); 
    %plot(I2,'g'); %plot iterations in green
    %hold on;
    
    fprintf('Iteration: %03d, STD: %e\n',itno,std(I2(:)-I(:)));
    if std(I2(:)-I(:)) < 1e-13; %convegence condition  
        stop=1;
       % plot(I2,'r'); %plot final distribution in red
        itno; %print out the iteration number
    else
       I = I2;
       Ip = gradient(I2,dz);
%         Ip(1)=0;
%        for a=1:1001
%           I(a) = I2(a);
%           if a>1&&a<1001; %take trapezoidal derivative
%               Ip(a) = 0.5*((I2(a)-I2(a-1)) + (I2(a+1)-I2(a)) );
%           else if a==1 %cant do this on the edges
%               Ip(a) = I2(a+1) - I2(a); 
%               else 
%                   Ip(a) = I2(a) - I2(a-1);
%               end
%           end
%        end
    end
    if itno==500; %stop after 100 iterations.
        stop =1; 
    end
end
% Fit a gaussian
temp = fitgaussian(I2,'plot',0,'bg_grad',0,'DC',0);
    %gaussianY = ones(sizey,1);
    xx = [1:length(I2)];
    gaussianY = temp.area * exp(-0.5*((xx-temp.mu)./((1+sign(xx-temp.mu)*temp.Assym_factor)*temp.sigma)).^2) / sqrt(2*pi*temp.sigma^2) + temp.bg_gradient*xx + temp.DC;
    sigmas(loop) = temp.sigma/c*1e12; % convert into ps
    %means(i) = Estimates(2);
    
    figure(201);
    clf;
    plot(z,I2,'b-'); hold on;
    plot(z,gaussianY,'-r');

end % for each current
save([pref '_hass'],'fit*','cal*','ave*','measured*','current','sigmas');

measuredsigma = avesigmas;
vectormeas = [avecurrents; measuredsigma];
vectorsim = [current; sigmas];
vectorerror = [avecurrents; measurederrorsigma];

figure(950); xlabel('Current (mA)'); ylabel('Bunch Sigma (arb. units)');hold on; 
plot(vectorsim(1,:),vectorsim(2,:),color);
%errorbar(vectormeas(1,17:19),vectormeas(2,17:19),vectorerror(2,17:19),'.'); 
errorbar(vectormeas(1,:),vectormeas(2,:),vectorerror(2,:),['.',color]);
%errorbar(vectormeas(1,17:19),vectormeas(2,17:19),vectorerror(2,17:19),'.'); 
 