% convergence of PTR. Also serves as Matlab warm up.
% Solutions. Barnett 6/5/14

a = 0.2;
f = @(s) a./(a^2+sin(s/2).^2);  % note elementwise ops
s=0:1e-3:2*pi; figure; plot(s,f(s),'-'); xlabel('s'); title('integrand');


Ns = 10:10:100; Is = 0*Ns;   % convergence study; Is will be data array
for i=1:numel(Ns), N = Ns(i);
  w = 2*pi/N*ones(1,N); s = 2*pi/N*(1:N);
  Is(i) = sum(w.*f(s));                     % do quadrature
  fprintf('N=%d\t I_N = %.16g\n',N,Is(i))
end

al = 2*imag(asin(1i*a)); % imaginary strip width of analyticity of f
figure; semilogy(Ns,abs(Is-Is(end)),'+-'); hold on; plot(Ns,exp(-al*Ns),'r-');

xlabel('N'); ylabel('quadrature error'); title('convergence of PTR');
v = axis; v(3:4)=[1e-16 1]; axis(v);
