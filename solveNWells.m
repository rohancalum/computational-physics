function [funcNWells, funcValsNWells] = solveNWells(psi_init, psiPrime_init, N)
% define constants
wellWidth = 0.6;
spacing = 0.2;
m = 1; 
V = 10;
hbar2 = 0.076199682;
alpha = @(E) sqrt(2*m*E/hbar2);
beta = @(E) sqrt(2*m*(V - E)/hbar2);
evals = linspace(0,V,10000);

 for i = 1:N  
    % initial 
    psi = @(E) psi_init(E);
    psiPrime = @(E) psiPrime_init(E);
    
    % compute psi(x) and psi'(x) at end of well
    psiWell_i = @(E) (cos(alpha(E).*(wellWidth)).*psi(E) + (1./alpha(E)).*sin(alpha(E).*(wellWidth)).*psiPrime(E));
    psiPrimeWell_i = @(E) (-alpha(E).*sin(alpha(E).*(wellWidth)).*psi(E) + cos(alpha(E).*(wellWidth)).*psiPrime(E));
    
    % propagate psi(x) and psi'(x) through forbidden region of well
    psiForbidden = @(E) cosh(beta(E).*(spacing)).*psiWell_i(E) + (1./beta(E)).*sinh(beta(E).*(spacing)).*psiPrimeWell_i(E);
    psiPrimeForbidden = @(E) beta(E).*sinh(beta(E).*(spacing)).*psiWell_i(E) + cosh(beta(E).*(spacing)).*psiPrimeWell_i(E);  
    
    % set psi(x) and psi'(x) to initial psi and psi', repeat N times
    psi_init = @(E) psiForbidden(E);
    psiPrime_init = @(E) psiPrimeForbidden(E);  
 end
 
 % extract value of psi(x) and psi'(x) BEFORE forbidden region of well N
 funcNWells = @(E) psiPrimeWell_i(E) +  beta(E).*psiWell_i(E);
 funcValsNWells = funcNWells(evals);
end 

