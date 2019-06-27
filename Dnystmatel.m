function a = Dnystmatel(G,i,j)
% DNYSTMATEL - return one element of 2D Laplace double-layer Nystrom matrix D
%
% Aij = Dnystmatel(G,i,j) returns A_{ij} given i,j, and a curve struct G
%  containing a boundary quadrature rule. No jump relation included for i=j.
% Barnett 6/7/14. Note this is tutorial, rather than efficient, code.

if i~=j
  d = G.x(i)-G.x(j); a = (1/(2*pi)) * real(conj(G.nx(j))*d)/abs(d)^2;
else
  a = -G.cur(j)/(4*pi);   % diagonal limit formula
end
%if strcmp(type,'exterior')
%    a = a + 1;
%end
a = a * G.sp(j)*G.w(j);   % speed weights
