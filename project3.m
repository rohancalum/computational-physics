% Variable mappings: 

% theta     ->      y1
% v         ->      y2*l
% alpha     ->      y2_prime

%% QUESTION 1: Unforced

% defining constants
A = 0;
m = 1;
l = 1;
g = 1;
w = 2/3;

% initial conditions 
%    theta, v
y0 = [0.2,  0];

% time interval 
% period 
T = 2*pi*sqrt(l/g);
tspan = [0 10*T];

nu_vals = [1 5 10];

for i = 1:length(nu_vals)
    % defining pair of 1st order equations
    fun = @(t,y) [y(2), (A/m)*sin(w*t)-(g/m)*sin(y(1)) - (nu_vals(i)/m)*y(2)];
    [y,t] = RK4(fun,tspan,y0);
    % plots
    figure;
    plot(t,y(:,1), '-r', 'Linewidth', 1.5);
    xlabel('time (seconds)', 'Interpreter', 'Latex');
    xlim([0 10*T]);
    ylabel('$\theta$ (rads)', 'Interpreter', 'Latex');
    title(['$\nu$ = ', num2str(nu_vals(i))], 'Interpreter', 'Latex');
    figure;
    plot(l*y(:,2),y(:,1),'-r', 'Linewidth', 1.5);
    xlabel('$v$ (rads/s)', 'Interpreter', 'Latex');
    ylabel('$\theta$ (rads)', 'Interpreter', 'Latex');
    title(['$\nu$ = ', num2str(nu_vals(i))], 'Interpreter', 'Latex');
end 

%% QUESTION 2: Driving Force
close all
% defining constants
m = 1;
l = 1;
nu = 1/2; 
g = 1;
w = 2/3;

% initial conditions 
%    theta, v
y0 = [0.2,  0];

% period 
T = 2*pi/w;
% time interval 
tspan = [0 300*T];
n=12;
A_vals = [0.5,1.2];

for i = 1:length(A_vals)
    
    % defining pair of 1st order equations
    fun = @(t,y) [y(2); (A_vals(i)/m)*sin(w*t)-(g/m)*sin(y(1)) - (nu/m)*y(2)];
    %[y,t] = RK4(fun,tspan,y0);
    options = odeset('RelTol',10^-n,'AbsTol',10^-n);
    [t,y] = ode45(fun,tspan,y0, options);
    
    
    ywrapped = y - 2*pi*floor( (y+pi)/(2*pi) ); 
      % plots
    figure;
    plot(t,ywrapped(:,1), '-r');
    xlim([0 10*T]);
    xlabel('time (seconds)', 'Interpreter', 'Latex');
    ylabel('$\theta$ (rads)', 'Interpreter', 'Latex');
    title(['$A$ = ', num2str(A_vals(i))], 'Interpreter', 'Latex');
    figure;
    plot(ywrapped(:,1), l*y(:,2), '-r');
    xlim([-pi pi]);
    ylabel('$v$ (m/s)', 'Interpreter', 'Latex');
    xlabel('$\theta$ (rads)', 'Interpreter', 'Latex');
    title(['$A$ = ', num2str(A_vals(i))], 'Interpreter', 'Latex');
end 

%% QUESTION 3:

% defining constants
m = 1;
l = 1;
nu = 1/2; 
g = 1;
w = 2/3;

% initial conditions 
%    theta, v
y0 = [0.2,  0];

% period 
T = 2*pi/w;
% time interval 
tspan = [0 300*T];

A_vals = [1.2];

for i = 1:length(A_vals)
    % defining pair of 1st order equations
    fun = @(t,y) [y(2), (A_vals(i)/m)*sin(w*t)-(g/m)*sin(y(1)) - (nu/m)*y(2)];
    [y,t] = RK4(fun,tspan,y0);
    % plots
    figure;
    plot(t,y(:,1));
    xlabel('time (seconds)', 'Interpreter', 'Latex');
    ylabel('$\theta$ (rads)', 'Interpreter', 'Latex');
    title(['A = ', num2str(A_vals(i))], 'Interpreter', 'Latex');
    figure;
    plot(y(:,1),l*y(:,2));
    xlim([-pi, pi])
    ylabel('$v$ (m/s)', 'Interpreter', 'Latex');
    xlabel('$\theta$ (rads)', 'Interpreter', 'Latex');
    title(['A = ', num2str(A_vals(i))], 'Interpreter', 'Latex');
end 

%% Question 4: A = 1.35, 1.44, 1.465

clear all
% defining constants
m = 1;
l = 1;
nu = 1/2; 
g = 1;
w = 2/3;

% initial conditions 
%    theta, v
y0 = [0.2,  0];

% period 
T = 2*pi/w;
% time interval 
tspan = [0 300*T];
n=11;
A_vals = [1.35,1.44,1.465];

for i = 1:length(A_vals)
    fun = @(t,y) [y(2); (A_vals(i)/m)*sin(w*t)-(g/m)*sin(y(1)) - (nu/m)*y(2)];
    options = odeset('RelTol',10^-n,'AbsTol',10^-n);
    [t,y] = ode45(fun,tspan,y0, options);
    
    ywrapped = y - 2*pi*floor( (y+pi)/(2*pi) ); 
    figure;
    plot(t,ywrapped(:,1), '-r');
    xlim([0 10*T]);
    xlabel('time (seconds)', 'Interpreter', 'Latex');
    ylabel('$\theta$ (rads)', 'Interpreter', 'Latex');
    title(['$A$ = ', num2str(A_vals(i))], 'Interpreter', 'Latex');
    figure;
    plot(ywrapped(:,1), l*y(:,2), '-r');
    xlim([-pi pi]);
    ylabel('$v$ (m/s)', 'Interpreter', 'Latex');
    xlabel('$\theta$ (rads)', 'Interpreter', 'Latex');
    title(['$A$ = ', num2str(A_vals(i))], 'Interpreter', 'Latex');
end 

%% Question 5: Poincare plots
clear all
% defining constants
m = 1;
l = 1;
nu = 1/2; 
g = 1;
w = 2/3;

% initial conditions 
%    theta, v
y0 = [0.2,  0];

% period 
T = 2*pi/w;
% time interval 
tspan = [0 300*T];
A_vals = [1.465]%0.5, 1.2,1.35,1.44,1.465]; %q2

NVALS = [10, 100, 1000, 2000];
for i = 1:length(A_vals)
    for j = 1:length(NVALS)
        t_vals =[0];
        for k = 1:NVALS(j)
            t = 2*pi*k/w;
            t_vals = [t_vals;t];
        end 
        fun = @(t,y) [y(2); (A_vals(i)/m)*sin(w*t)-(g/m)*sin(y(1)) - (nu/m)*y(2)];
        options = odeset('RelTol',10^-11,'AbsTol',10^-11);
        [t,y] = ode45(fun,t_vals,y0, options);
        ywrapped = y - 2*pi*floor((y+pi)/(2*pi));
    
    hold on
    subplot(2,2,j)
    plot(ywrapped(:,1), l*y(:,2), 'o');
    xlim([-pi pi]);
    ylabel('$v$ (m/s)', 'Interpreter', 'Latex');
    xlabel('$\theta$ (rads)', 'Interpreter', 'Latex');
    title(['n = ', num2str(NVALS(j)), ' and A = ', num2str(A_vals(i))], 'Interpreter', 'Latex');  
    end 
        figure;
        plot(0:1:2000,ywrapped(:,1), 'o');
        xlabel('n', 'Interpreter', 'Latex');
        ylabel('$\theta(2 \pi n / \omega)$ (rads)', 'Interpreter', 'Latex');
        title(['A = ', num2str(A_vals(i))], 'Interpreter', 'Latex'); 
    
end 

    
%% QUESTION 6: Large n behaviour 

% Poincare plots for 0.5 ? A ? 1.2 and 1.35 ? A ? 1.5

largeN = 1100;
t_vals =[0];
for k = 1000:1:largeN
    t = 2*pi*k/w;
    t_vals = [t_vals;t];
end

Avals1 = 0.5:0.01:1.2;
Avals2 = 1.35:0.01:1.5;

thetaMatrix = zeros(length(t_vals),length(Avals2));

for i = 1:length(Avals2)
    fun = @(t,y) [y(2); (Avals2(i)/m)*sin(w*t)-(g/m)*sin(y(1)) - (nu/m)*y(2)];
    options = odeset('RelTol',10^-11,'AbsTol',10^-11);
    [t,y] = ode45(fun,t_vals,y0, options);
    ywrapped = y - 2*pi*floor((y+pi)/(2*pi));
    yAtN = ywrapped(:,1);
    thetaMatrix(:,i) = yAtN
end 


for column = 1:length(Avals2)
    hold on 
    plot(Avals2(column)*ones(length(t_vals)-1,1), thetaMatrix(2:end,column),'.');
end 
xlabel('A', 'Interpreter', 'Latex');
ylabel('$\theta (2 \pi n / \omega)$', 'Interpreter', 'Latex');

    