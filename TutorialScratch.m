k = -3;
f = @(s) exp(1i*k*s);
s = -1000:1:1000;
Q = integral(f, 25, 25)

%This can be shown analytically by directly solving for integer bounds
%that are equal, but how do we bridge the gap to infinite bounds? 

%%
%Exercise 5
p = 16;
f = @(x) (x.^(-0.5));
n = 20;

%acquire nodes x and weights w for Gauss quadrature of p points
I = 0;
for i=0:1:n
    left_endpt = (0.5/2^i);
    right_endpt = (1.0/(2^i));
    [x,w] = gauss(p);
    intvlChngFacOne = (right_endpt - left_endpt)/2; 
    intvlChngFacTwo = (right_endpt + left_endpt)/2;
    I = I + intvlChngFacOne*sum(w.*f(intvlChngFacOne*x+intvlChngFacTwo));
end

fprintf("%f",I) 

%Compute the actual integrand and the answer is 2.
%%
%HELP

%Do not know how to find the alpha for semiaxes sinh alpha and cosh alpha
alpha = 1 ;
%Check the convergence rate with N = np
exp(-2*alpha*p*n);
%what does it mean that convergence is algebraic of the high order 2p if f
% is in C^2p ? 
%%
%Exercise 6


k = 2 ;
f = @(s) exp(1i*k*s);
n = 20;
p = 16;
%%
%composite Gauss-legendre rules
I = 0;
i = 0:1:20-1 ;
left_endpt = (0.5./(2.^i));
right_endpt = (1.0./(2.^i));
[x,w] = gauss(p);

intvlChngFacOne = ((right_endpt - left_endpt)/2)'; 
intvlChngFacTwo = ((right_endpt + left_endpt)/2)';
y= intvlChngFacTwo*zeros(1,p);
z = f(intvlChngFacOne*(x')+y) ;
w = zeros(n,1)*w';
sumOverP = sum(w.*z,2); 
IGauss = sum(intvlChngFacOne.*sumOverP);
fprintf("%f",IGauss) 

%%
%PTR rules 
N = 16*20;
w = 2*pi/N;
j = 1:N;
sj = 2*pi*j/N ; 
I_PTR = sum(w*f(sj));

%%
G = smoothstar(0.3,5)
G = curvquad(G,'ptr',100,16)


%%
sigma = @(s) 1;
G = smoothstar(0.4,5);
G = curvquad(G,'panel', 70, 10);
targets = (-2:0.5:2);

u = evalDLP_mine(sigma,G,targets)






