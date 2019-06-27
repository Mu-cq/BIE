function u = evalDLPhelm(t,G,sigma,k,eta)
% EVALDLPHELM - evaluate 2D Helmholtz CFIE potential on target points
%
% u = evalDLPhelm(t,G,sigma,k,eta) where t is a list of target points
%  (points in the complex plane), G is a curve struct (see curvquad.m) with
%  N quadrature nodes, sigma is a N-element column vector, k is the wavenumber,
%  and eta (~k) controls the amount of SLP in the CFIE, returns the potential
%  at the target points evaluated using the quadrature rule in G.
%  A direct sum is used, vectorized over target points only.
%
% Barnett 6/8/14

N = numel(G.x);
u = 0*t;        % initialize output
for j=1:N
  d = t-G.x(j);  % displacement of targets from jth src pt
  kr = k*
  u = u + k*ct*besselh(1,kr) + 1i*eta*besselh(0,kr);
  
  (G.w(j) * G.sp(j) * sigma(j)) * real(conj(G.nx(j)).*d)./abs(d).^2;
end
u = (1i/4)*u;
