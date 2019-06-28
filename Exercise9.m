%Exercise 9 Done in 2D (no conversion to complex)

%This does not work. 

%choosing interval [-1,1] 

G = curveCircle(1,400); %function on a circle of radius 1
x=-1.5:0.05:1.5;
y=-1.5:0.05:1.5; %sampling points inside and outside of omega 
sigma = @(x) 1; %constant density should yield -1 in omega, 0 outside 
result = evaluateDLP2d(sigma,G, [x;y]);


%%
