% test the DLP Laplace 2D evaluator evalDLP, convergence at a point. 
% Barnett 6/21/14

clear
G = smoothstar(.3,5);                 % curve shape
p = 3;          % # nodes per panel 
t = 0.2+-.3i;    % target
Ns = 20:20:500; e = 0*Ns;
for i=1:numel(Ns), N = Ns(i);
  G = curvquad(G, 'panel', N,p);   % curve quadrature rule
  sigma = ones(size(G.s));         % density; note N may not be what requested
  u = evalDLP(t,G,sigma);               % use the evaluator you will write
  e(i) = u+1;  % error
end
figure; loglog(Ns,abs(e),'+-'); hold on; plot(Ns,(Ns/25).^(-2*p),'r-');
plot(Ns,(Ns/25).^-p,'m-');
axis tight; v = axis; v(3) = 1e-16; axis(v);
title('interior Laplace DLP ptwise error convergence');
% not clearly order 2p (predicted) unless p around 6. Large p it's closer to p.
% small p (<5) it becomes exponential to do with PTR.
