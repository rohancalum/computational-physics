format long

%% Question 1

X = 0:1/10000:1;

% set optimal h
h = (10^-16)^(1/7)

% function to take derivative
FofX = @(x) sin(x.^2);
FofXPrime = @(x) 2*x.*cos(x.^2);
Fprime = FofXPrime(X);

% compute P_6(x)
p_x = [];
for x = 0:1/10000:1
    
    % 7-points to take derivative at
    xi = [x-3*h, x-2*h, x-h, x+h, x+2*h, x+3*h];
    
    % Lagrange interpolant coefficients from derivation
    L1 = (-1/60);
    L2 = (3/20);
    L3 = (-3/4);
    L4 = (3/4);
    L5 = (-3/20);
    L6 = (1/60);
    p = (1/h)*(L1*FofX(xi(1)) + L2*FofX(xi(2)) + L3*FofX(xi(3)) + L4*FofX(xi(4)) + L5*FofX(xi(5)) + L6*FofX(xi(6))); 
    p_x = [p_x,p];
end

figure(1)
subplot(3,1,1)
plot(X, Fprime)
title('Analytical Derivative: f`(x)= 2xcos(x)')

subplot(3,1,2)
plot(X,p_x)
title('Numerical Derivative: P_6(x)')

subplot(3,1,3)
error = abs(Fprime - p_x);
plot(X, error)
title('Error: f`(x) - P_6(x)')


%% Question 2

% set bounds
a = 0 ;
b = 1 ;

% function to integrate
f = @(x) sin(x)

% romberg approximation elements
R11 = ((b-a)/2)*(f(a)+f(b))
R21 = ((b-a)/4)*(f(a)+f(b)) + ((b-a)/2)*(f((a+b)/2))
R31 = ((b-a)/8)*(f(a)+f(b)) + ((b-a)/4)*(f((3*a+b)/4) + f((a+b)/2) + f((a+3*b)/4))
R22 = R21 + (R21-R11)/3;
R32 = R31 + (R31-R21)/3;

% 6th order accurate approximation
R33 = R32 + (R32-R22)/15

% value of analytical function
I = -cos(1) + cos(0)

% error
error =  abs(I-R33)













