% Function to calculate average energy and magnetization for the 1D Ising model
% T = Temperature, N = linear lattice size, J = Ising coupling.
function [E,M, Elist, Mlist, M_thermalized, E_thermalized, firstThermalValueIndex, Chi_m C_v]=ising2D(T,N,J,plot_flag)
%%  Initial configuration
grid = sign(.5-rand(N,N)); % Random initial configuration
%%  Initiation
t =5e5*N*N; % Number of steps
Elist=zeros(t,1);
Mlist=zeros(t,1);
%periodic boundary conditions
x = circshift(grid,[0 1])+circshift(grid,[0 -1])+ circshift(grid,[1 0])+ circshift(grid,[-1 0]);
Energy = -J*sum(sum(grid.*x)); %initial Energy
Magnet=sum(sum(grid)); %initial magnetization

%2D
trials_x = randi(N,t,1); %cheaper to generate all at once
trials_y = randi(N,t,1); %cheaper to generate all at once

%%  Metropolis algorithm
for i=1:t
   
    %2D
    s_x= trials_x(i);
    s_y= trials_y(i);
    % account for four directions
    % if on left edge
    if s_x ~= 1 
        west=grid(s_x-1, s_y);
    else 
        west=grid(N,s_y);
    end 
     % if on right edge
    if s_x ~= N 
        east=grid(s_x+1, s_y);
    else 
        east=grid(1,s_y);
    end 
    % if on top edge
    if s_y ~= 1 
        north=grid(s_x, s_y-1);
    else 
        north=grid(s_x,N);
    end 
    % if on bottom edge
    if s_y ~= N 
        south=grid(s_x, s_y+1);
    else 
        south=grid(s_x,1);
    end 
    dE=2*J*grid(s_x,s_y)*(west+north+east+south);  % change in energy
    p = exp(-dE/T);
    % Acceptance test (including the case dE<0).
    if rand <= p
        grid(s_x,s_y) = -1*grid(s_x,s_y); 
        Energy=Energy+dE;
        Magnet=Magnet+2*grid(s_x,s_y);
    end
    % Update energy and magnetization.
    Mlist(i) =Magnet;
    Elist(i) =Energy;
end
%% Display time series of energy and magnetization
Elist(Elist==0)=[];
Mlist(Mlist==0)=[];
Mlist=abs(Mlist);
Mlist=Mlist/(N*N);
Elist=Elist/(N*N);    %normalize.

%%  Magnetization and energy density
% Obtain list of all configurations before thermalization.
Elist(Elist==0)=[];Mlist(Mlist==0)=[];
M_thermalized = Mlist(t/2:end); % t/2 is well-after thermalization
E_thermalized = Elist(t/2:end);

% Average over post thermalization configurations.
M=sum(Mlist)/length(M_thermalized);
E=sum(Elist)/length(E_thermalized);

%% Find thermalization time 

% compute the mode of the smoothed E_thermalized values
thermalizedValue = mode(E_thermalized);

% find energy difference between computed energies and thermalized value 
energyDifference = abs(Elist - thermalizedValue);

% determine where this difference is less than a specified tolerance
thermalValues = energyDifference < 1e-12;

% the first occurance of these values will be the minimum themalization time
firstThermalValueIndex = find(thermalValues,1);

%% Physical quantities of interest

% take our energy and magnetization values after thermalization point 
Elist2 = Elist(firstThermalValueIndex:end);
Mlist2 = Mlist(firstThermalValueIndex:end);

% Magnetic susceptibility
Chi_m = (sum(Mlist2.^2)/t - sum(Mlist2)^2/(t^2))/T;

% Heat capacity
C_v = (sum(Elist2.^2)/t - sum(Elist2)^2/(t^2))/(T^2);

% Magnetization
M=sum(Mlist2)/(length(Mlist2)*N^2); 

% Energy 
E=sum(Elist2)/(length(Elist2)*N^2);
end
