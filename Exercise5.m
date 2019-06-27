
%%
%Exercise 5
p = 16;
f = @(x) (x.^(-0.5));
err = 0*8;
Ns = 20:10:100;
alpha_pred = 0*length(Ns); 
alpha_obs = 0*(length(Ns));
j=1;
for n= Ns
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

    fprintf("Estimation: %f",I) 

    %Actual integrand analytically evaluates to 2
    actual = 2;
        
    err(j) = abs((I - actual));
    alpha_obs(j) = log(abs(I-actual))/(-2*n*p);
    alpha_pred(j) = (log(2) + log(sqrt(2^(n/2))))/(2*n*p);
    fprintf("predicted_alpha: %f\n",alpha_pred(j));
    fprintf("observed_alpha: %f\n", log(abs(I-actual))/(-2*n*p)); 
    j = j + 1;

end

figure;
semilogy(Ns, alpha_obs, 'b');
hold on; plot(Ns, err, 'r');
hold on;  plot(Ns, alpha_pred, 'y');
xlabel('Ns'); title('Convergence Rate'); ylabel('Quadrature Error');




