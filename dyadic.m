% Solution to dyadic refinement exercise. Barnett 6/5/14
clear;
f = @(s) s.^(-1/2);  % singular function
I = 2;   % its exact integral
p=16;
[x w] = gauss(p); x = (x+3)/2; w = w/2; % shift to [1,2]
ns = 5:5:100; Is = 0*ns;
for i=1:numel(ns), n=ns(i); N = p*n;   % convergence loop
  for j=n:-1:1, Is(i) = Is(i) + sum(w.*f(x/2^j))/2^j; end % loop over panels
  fprintf('N=%d\t I_N = %.16g\n',N,Is(i))
end
figure; semilogy(ns,abs(Is-I),'+-'); xlabel('N'); ylabel('quadrature error');
% Note we take advantage of high relative accuracy in nodes close to zero

