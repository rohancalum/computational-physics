% Function to calculate average energy and magnetization for the 1D Ising model
% T = Temperature, N = linear lattice size, J = Ising coupling.
function [E,M,chi,Cv, Elist, Mlist]=ising1D(T,N,J,nn,a,plot_flag)
%%  Initial configuration
grid = sign(.5-rand(N,1)); % Random initial configuration
%%  Initiation
iter = 2e4;
t =iter*N;

Elist=zeros(t,1);
Mlist=zeros(t,1);
Cv=zeros(t,1);
chi=zeros(t,1);
pos = 1:N;

%% Initial Energy
Energy = 0;
if (nn == 1)
    Energy = -J*sum(grid.*circshift(grid,1));
else
    for i = 1:N
        %fprintf('\ni: %d\n', i);
        for j = 1:i-1
            d = j - i;
            if (d < 0)
                d = min(abs(d),abs(d+N));
            end
            if (d > 0)
                d = min(abs(d),abs(d-N));
            end
            if (j ~= i)
                Energy = Energy  - J*(d)^(-a)*grid(i)*grid(j);
            end
        end
    end
end

Magnet=sum(grid); %initial magnetization
trials = randi(N,t,1); %cheaper to generate all at once

%%  Metropolis algorithm
for i=1:t,
    s=trials(i);
    if s~=1; left=grid(s-1);else left=grid(N);end
    if s~=N; right=grid(s+1);else right=grid(1);end
    
    dE=0;  % change in energy
    old_Energy = Energy;
    if (nn == 1)
        dE=2*J* grid(s)*(left+right);  % change in energy
    else
        for j = 1:N
            d = j - s;
            if (d < 0)
                d = min(abs(d),abs(d+N));
            end
            if (d > 0)
                d = min(abs(d),abs(d-N));
            end
            if (j ~= s)
                dE = dE + 2*J*(d)^(-a)*grid(j)*grid(s);
            end
        end
        
        %% Not any faster in practice
        % d_array = pos - s;
        % d_array(d_array < 0) = min(abs(d_array(d_array < 0)), abs(d_array(d_array < 0) + N));
        % d_array(d_array > 0) = min(abs(d_array(d_array > 0)), abs(d_array(d_array > 0) - N));
        % d_array(s) = Inf;
        
        %dE = sum(2*J*((d.^(-a)).*grid)*grid(s));
    end
    
    p = exp(-dE/T);
    % Acceptance test (including the case dE<0).
    if rand <= p,
        grid(s) = -1*grid(s);
        Energy=Energy+dE;
        Magnet=Magnet+2*grid(s);
    end
    % Update energy and magnetization.
    Mlist(i) =Magnet;
    Elist(i) =Energy;
    %Refresh display of spin configuration every N trials.
    %         if mod(i,N)==0 && plot_flag==1;
    %             bar(grid); drawnow;
    %         end
end

%% Double check that the Energy is what it should be
%% If your energy graph has a sudden dip at the end, somthing went wrong
Energy = 0;
if (nn == 1)
    Energy = -J*sum(grid.*circshift(grid,1));
else
    for i = 1:N
        %fprintf('\ni: %d\n', i);
        for j = 1:i-1
            d = j - i;
            if (d < 0)
                d = min(abs(d),abs(d+N));
            end
            if (d > 0)
                d = min(abs(d),abs(d-N));
            end
            if (j ~= i)
                Energy = Energy  - J*(d)^(-a)*grid(i)*grid(j);
            end
        end
    end
end
Elist(end) = Energy;

%% Display time series of energy and magnetization
Elist(Elist==0)=[];
Mlist(Mlist==0)=[];
% Eliminate all configurations before thermalization.
Elist_trunc = Elist(t/2:end);
Mlist_trunc = Mlist(t/2:end);
% Magnetic susceptibility
chi=(sum(Mlist_trunc.^2)/numel(Mlist_trunc)-sum(Mlist_trunc)^2/numel(Mlist_trunc)^2)/T/N;
% Heat capacity
Cv=(sum(Elist_trunc.^2)/numel(Elist_trunc)-sum(Elist_trunc)^2/numel(Elist_trunc)^2)/(T^2)/N;

%normalize
Mlist=Mlist/N;
Elist=Elist/N;

Mlist_trunc = Mlist_trunc/N;
Elist_trunc = Elist_trunc/N;

if plot_flag==1
    figure;
    subplot(2,1,1)
    plot(Elist)
    subplot(2,1,2)
    plot(Mlist)
end
%%  Magnetization and energy density
% Average over post thermalization configurations.
M=sum(Mlist_trunc)/numel(Mlist_trunc);
M=sum(abs(Mlist_trunc))/numel(abs(Mlist_trunc));
E=sum(Elist_trunc)/numel(Elist_trunc);
end