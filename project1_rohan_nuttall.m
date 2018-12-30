clear all
% defining constants
S = 0.2;
a = 0.3; 
m = 1; 
V = 10; 
hbar2 = 0.076199682;
alpha = @(E) sqrt(2*m*E/hbar2);
beta = @(E) sqrt(2*m*(V - E)/hbar2);
ENERGY = linspace(0,V,1000);

% analytical solution to find roots of (even and odd)
f_even = @(E) beta(E).*cos(alpha(E)*a) - alpha(E).*sin(alpha(E)*a);
f_odd = @(E) alpha(E).*cos(alpha(E)*a) + beta(E).*sin(alpha(E)*a); 
f_odd_Prime = @(E) beta(E).^2.*cos(alpha(E)*a) - alpha(E).*sin(alpha(E)*a);
f_even_Prime = @(E) -beta(E).*alpha(E).*sin(alpha(E)*a) - alpha(E).^2.*cos(alpha(E)*a);

% plot even function 
figure(1)
subplot(2,1,1)
plot(ENERGY,f_even(ENERGY), ENERGY, zeros(length(ENERGY),1), 'r')
ylabel('f(E)_{even}')   
xlabel('E [arb. units]')
ylim([-15, 15])

% plot odd function
subplot(2,1,2)
plot(ENERGY,f_odd(ENERGY), ENERGY, zeros(length(ENERGY),1), 'r')
ylabel('f(E)_{odd}')
xlabel('E [arb. units]')
ylim([-15, 15])

% find roots (using Newton's method) of analytical solution based on guesses from plots
guess1_even = 1;
guess2_even = 7;
guess_odd = 4;
root1_even = newton(f_even, f_even_Prime, guess1_even)
root2_even = newton(f_even, f_even_Prime, guess2_even)
root2_odd  = newton(f_odd, f_odd_Prime, guess_odd)

% QUESTION 1

% define constants
x = 2*a;
b = 0;
wellWidth = x-b;

% set old psi vector as a function of energy 
psi_b = @(E,dS) 1;
psiPrime_b = @(E,dS) beta(E);

% compute new psi vector as a function of energy
psiWell_1 = @(E,dS) (cos(alpha(E).*(wellWidth)).*psi_b(E) + (1./alpha(E)).*sin(alpha(E).*(wellWidth)).*psiPrime_b(E));
psiPrimeWell_1 = @(E,dS) (-alpha(E).*sin(alpha(E).*(wellWidth)).*psi_b(E) + cos(alpha(E).*(wellWidth)).*psiPrime_b(E));

% define f(E, dS)
f = @(E, dS) psiPrimeWell_1(E) + beta(E).*psiWell_1(E);

% plot numerical solution of f(E) for a single well 
figure(2)
plot(ENERGY,f(ENERGY), ENERGY, zeros(length(ENERGY),1), 'r') 
ylabel('f(E) = \psi''(c) + \beta*\psi(c)')
xlabel('E [arb. units]')

% compute zeros of f(E) numerically for single well (with guesses)
root1_WELL1 = fzero(f,1)
root2_WELL1 = fzero(f,3)
root3_WELL1 = fzero(f,7)

% compute numerical accuracy9
accuracy1 = (root1_WELL1/root1_even)*100
accuracy2 = (root2_WELL1/root2_odd)*100
accuracy3 = (root3_WELL1/root2_even)*100


% QUESTION 2 

% set initial well spacing
s = 0.2;

% compute psi and psi' in forbidden region 
psiForbidden = @(E,dS) cosh(beta(E).*(dS)).*psiWell_1(E,dS) + (1./beta(E)).*sinh(beta(E).*(dS)).*psiPrimeWell_1(E,dS);
psiPrimeForbidden = @(E,dS) beta(E).*sinh(beta(E).*(dS)).*psiWell_1(E,dS) + cosh(beta(E).*(dS)).*psiPrimeWell_1(E,dS);

% compute psi and psi' in allowed region of well 2
psiWell_2 = @(E,dS) (cos(alpha(E).*(wellWidth)).*psiForbidden(E,dS) + (1./alpha(E)).*sin(alpha(E).*(wellWidth)).*psiPrimeForbidden(E,dS));
psiPrimeWell_2 = @(E,dS) (-alpha(E).*sin(alpha(E).*(wellWidth)).*psiForbidden(E,dS) + cos(alpha(E).*(wellWidth)).*psiPrimeForbidden(E,dS));

% define f(E) for double well as a function of energy and well spacing 
fDoubleWell = @(E, dS) psiPrimeWell_2(E,dS) + beta(E).*psiWell_2(E,dS);

% plot f(E) of the double well with a well spacing of s
plot(ENERGY,fDoubleWell(ENERGY,s), ENERGY, zeros(length(ENERGY),1), 'r')
xlabel('Energy [arb. units]')
ylabel('f(E)')
ylim([-50,500])

% redefine f(E) in terms of one variable
fDoubleWell_s = @(E) fDoubleWell(E,s);

% numerical roots of f(E) for double well 
root1_WELL2 = fzero(fDoubleWell_s,0.7)
root2_WELL2 = fzero(fDoubleWell_s,0.8)
root3_WELL2 = fzero(fDoubleWell_s,2.5)
root4_WELL2 = fzero(fDoubleWell_s,3.0)
root5_WELL2 = fzero(fDoubleWell_s,6.0)
root6_WELL2 = fzero(fDoubleWell_s,6.5)
root7_WELL2 = fzero(fDoubleWell_s,9.7)

% compute lowest two energies as a function of well spacing
slist = linspace(0, 0.2, 100);
lowestTwoEnergies_candidates = zeros(2,length(slist));
lowestTwoEnergies_roots = zeros(2,length(slist));
for i = 1:length(slist)
    func = @(E) fDoubleWell(E, slist(i));
    funcToEval = fDoubleWell(ENERGY, slist(i));
    possibleRoots = [];
    % bisection root finding method
    for x = 1:length(ENERGY)-1
        if funcToEval(x)*funcToEval(x+1) < 0
            possibleRoots = [possibleRoots;ENERGY(x)];
        end
    end
    % load lowest two possible roots into candidates matrix
    lowestTwoEnergies_candidates(:,i) = possibleRoots(1:2);
    
    % compute computed energy values for each well spacing
    lowestTwoEnergies_root1 = fzero(func, lowestTwoEnergies_candidates(1,i));
    lowestTwoEnergies_root2 = fzero(func, lowestTwoEnergies_candidates(2,i));
    
    % load computed energy values into roots matrix
    lowestTwoEnergies_roots(1,i) = lowestTwoEnergies_root1;
    lowestTwoEnergies_roots(2,i) = lowestTwoEnergies_root2;

end 

lowestTwoEnergies_candidates;

%investigate why the roots of the 2nd energy level jump down from 0.73 to
%0.70 in one step around a spacing of ~0.18nm. Must be something off in the
%root finding method...

lowestTwoEnergies_roots;

% plot 1st lowest two energy
figure()
plot(slist, lowestTwoEnergies_roots(1,:))
hold on 
% plot 2nd lowest two energy
plot(slist, lowestTwoEnergies_roots(2,:))
xlabel('Spacing (nm)')
ylabel('Energy (eV)')

% plot dE = E_2 - E 1
figure()
plot(slist, lowestTwoEnergies_roots(2,:) - lowestTwoEnergies_roots(1,:))
xlabel('Spacing (nm)')
ylabel('Energy Difference E_2 - E_1 (eV)')

% QUESTION 3 

% defining constants
N=4;
psi_b = @(E) 1;
psiPrime_b = @(E) beta(E);
evals = linspace(0,V,10000);

% call solveNWells function
[funcNWells, funcValsNWells] = solveNWells(psi_b, psiPrime_b, N);

% plot f(E) for N-wells
figure()
plot(evals, funcValsNWells, evals, zeros(length(evals),1), 'r')
title(sprintf('F(E) of the N=%d well',N))
ylim([-50,500])

% search for roots
possibleRoots = [];
for x = 1:length(evals)-1
    if funcValsNWells(x)*funcValsNWells(x+1) < 0
        possibleRoots = [possibleRoots;evals(x)];
    end
end

% compute roots
actualRoots = [];
for i = 1:length(possibleRoots)
    actualRoots = [actualRoots; fzero(funcNWells,possibleRoots(i))];
end 
actualRoots