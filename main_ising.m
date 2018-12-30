%% Question 2 
N = 10;
J = 1;
plot_flag = 1;   

% averaging
a = [];
for i = 1:100 
   [E, M, Elist, Mlist, M_thermalized, E_thermalized,firstThermalValueIndex] = ising2D(T,N,J,plot_flag);
   a = [a;firstThermalValueIndex];
end 
b = mean(a)

%% Question 3
t_vals = (2:0.02:2.35);
C_v_vals = [];
Chi_m_vals = [];
Evals = [];
Mvals = [];

% iterate over time values 
for i = 1:length(t_vals)
    [E, M, Elist, Mlist, M_thermalized, E_thermalized,firstThermalValueIndex, C_v, Chi_m] = ising2D(t_vals(i),N,J,plot_flag);
    C_v_vals = [C_v_vals;C_v];
    Chi_m_vals = [Chi_m_vals;Chi_m];
    Evals = [Evals;E];
    Mvals = [Mvals;M];
    fprintf(1,'Iteration #: %d \n',i)
end 
    
subplot(2,2,1)
plot(t_vals, Chi_m_vals, '^')
title('$\Chi_m$(T)', 'Interpreter','latex')
xlabel('T')

subplot(2,2,2)
plot(t_vals, C_v_vals,'^')
title('C_v(T)')
xlabel('T')

subplot(2,2,3)
plot(t_vals, Evals,'^')
title('E(T)')
xlabel('T')

subplot(2,2,4)
plot(t_vals, Mvals,'^')
title('M(T)')
xlabel('T')

%% Question 4
plot(t_vals, Chi_m_vals, '^')
title('$\Chi_m$(T)', 'Interpreter','latex')
xlabel('T')