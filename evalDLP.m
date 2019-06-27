function u = evalDLP(t,G,sigma)
% EVALDLP - evaluate 2D Laplace double-layer potential on target points
%
% u = evalDLP(t,G,sigma) where t is a list of target points (points in the
%  complex plane), G is a curve struct (see curvquad.m) with N quadrature
%  nodes, and sigma is a N-element column vector, returns the potential
%  at the target points evaluated using the quadrature rule in G.
%  A direct sum is used, vectorized over target points only.
%
% Barnett 6/6/14

N = numel(G.x);
u = 0*t;        % initialize output
for j=1:N
  d = t-G.x(j);  % displacement of targets from jth src pt
  u = u + (G.w(j) * G.sp(j) * sigma(j)) * real(conj(G.nx(j)).*d)./abs(d).^2;
end
u = u/(2*pi);
