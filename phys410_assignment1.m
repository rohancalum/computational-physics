
% Question 1

set step size to search for vicinity of roots. 
format long 
step = 1/50;
guesses = [];
for i = -1:step:1-step
    if F(i)*F(i+step)>0
    else
        guesses = [guesses;i,i+step]; 
    end 
end 
guesses

%
newtonInput = [];
for i = 1:length(guesses)
    newtonInput = bisection(guesses(i,1), guesses(i,2));
    roots = [roots; newton(F, Fprime, newtonInput)]
end 
newtonInput


function root=newton(F, Fprime,guess)
iter=1;
maxiter=1000;
x=guess;
accuracy=1;
tolerance=1e-12;
while accuracy>tolerance && iter<maxiter;
    x=x-feval(F,x)/feval(Fprime,x);
    accuracy=abs(feval(F,x));
end
root=x;
end

Simple bisection for finding roots of functions of a single variable.
Use as starting point, add comments and elaborations of basic code as needed.
FofX the function whose root we are finding (example below), to be included in the Matlab path/
function root=bisection(xL,xR)
input checking
if F(xL)*F(xR)>0; disp('Initial interval contains even number of roots');
    return;
end
tolerance = 1e-12;
accuracy=100*tolerance;
while accuracy > tolerance
    xmiddle=(xL+xR)/2;
    if F(xL)*F(xmiddle)<0
        xR=xmiddle;
    else
        xL=xmiddle;
    end
    accuracy=abs(F(xmiddle));
end
root=xmiddle;
end


function out = F(x)
out = (6425*x^8 - 12012*x^6 + 6930*x^4 - 1260*x^2 + 35)/128;
end

function out = Fprime(x)
out = (6425*x^6 - 9009*x^4 + 3465*x^2 - 315)/16;
end 

    


