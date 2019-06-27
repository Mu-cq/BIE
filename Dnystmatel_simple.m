function a = Dnystmatel(G,i,j)
% DNYSTMATEL - return one element of 2D Laplace double-layer Nystrom matrix D
if i~=j
  d = G.x(i)-G.x(j); a = (1/(2*pi)) * real(conj(G.nx(j))*d)/abs(d)^2;
else
  a = -G.cur(j)/(4*pi);   % diagonal limit formula
end
a = a * G.sp(j)*G.w(j);   % speed weights
