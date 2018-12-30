ising2D(2.15, 50, 1, 1)

% Function to calculate average energy and magnetization for the 1D Ising model
% T = Temperature, N = linear lattice size, J = Ising coupling.
function [E,M]=ising2D(T,N,J,plot_flag)
%%  Initial configuration
grid = sign(.5-rand(N,N)); % Random initial configuration
%%  Initiation
t =1e4*N; % Number of steps
Elist=zeros(t,1);
Mlist=zeros(t,1);
x = circshift(grid,[0 1])+circshift(grid,[0 -1])+ circshift(grid,[1 0])+ circshift(grid,[-1 0]);
Energy = -J*sum(sum(grid.*x)) %initial Energy
Magnet=sum(sum(grid)); %initial magnetization
trials_x = randi(N,t,1); %cheaper to generate all at once
trials_y = randi(N,t,1); %cheaper to generate all at once

%%  Metropolis algorithm
for i=1:t
    s_x= trials_x(i);
    s_y= trials_y(i);
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
    % Refresh display of spin configuration every N trials.
%             if mod(i,N)==0 && plot_flag==1
%                 bar(grid); drawnow;
%             end
%     
end
%% Display time series of energy and magnetization
Elist(Elist==0)=[];Mlist(Mlist==0)=[];
Mlist=abs(Mlist);
Mlist=Mlist/(N*N);
Elist=Elist/(N*N);    %normalize.
if plot_flag==1
    figure; 
    subplot(2,1,1)
    plot(Elist)
    title('Elist')
    subplot(2,1,2)
    plot(Mlist)
    title('Mlist')
end
figure;
image(grid,'CDataMapping','scaled')

% TODO: after performing this procedure, Energyf should be the same as the last
% element in the Elist... this is not the case. Find out why.
xf = circshift(grid,[0 1])+circshift(grid,[0 -1])+ circshift(grid,[1 0])+ circshift(grid,[-1 0]);
Energyf = -J*sum(sum(grid.*xf))
Elistend = Elist(end)
%%  Magnetization and energy density
% Eliminate all configurations before thermalization.
Mlist(1:50*N)=[];
Elist(1:50*N)=[];
% Average over post thermalization configurations.
M=sum(Mlist)/numel(Mlist);
E=sum(Elist)/numel(Elist);
end
