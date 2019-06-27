%% 
%Exercise 6
k = 1 ;
f = @(s) exp(1i*k*s);
p = 16; %points
n = 50; %panels 


fprintf("For N = %d\n", p*n);
%%

%Interval [0,1000]
%composite Gauss-legendre rules
i = 0:1:n-1 ;
left_endpt = (20*i);
right_endpt = (20*(i+1));
[x,w] = gauss(p);

intvlChngFacOne = ((right_endpt - left_endpt)/2)'; 
intvlChngFacTwo = ((right_endpt + left_endpt)/2)';

sumOverPanels = 0;

for j=1:1:n %looping over panels
    f_eval = f(intvlChngFacOne(j)*(x)+intvlChngFacTwo(j)) ;    
    f_eval = w.*f_eval;
    sumOverPoints = intvlChngFacOne(j)*sum(f_eval); 
    sumOverPanels = sumOverPanels+sumOverPoints;
end   

IGauss = sumOverPanels;
fprintf("Gauss-Legendre Estimate: %f\n",IGauss) 



%%
I_actual = integral(f,0,1000);
fprintf("Integral Answer: %f\n",I_actual) 

%%
%%
%PTR rules 
k=1;
epsilon = 1e-2;
while (1)
    N = p*n*k;
    w = 1000/N;
    j = 1:N;
    sj = 1000*j/N ; 
    I_PTR = sum(w*f(sj));
    if (abs(IGauss-I_PTR) < epsilon)
        break;
    end
    k = k + 1;
end

fprintf("%d N required before PTR converges to N point Gauss Legendre Estimation with epsilon %f \n",k, epsilon); 
fprintf("PTR Rules Estimate: %f\n",I_PTR)


