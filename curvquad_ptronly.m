function G = curvquad(G,rule,N)
% CURVQUAD - set up underlying quadrature for a closed curve struct
if strcmp(rule,'ptr')
  G.s = 2*pi*(1:N)'/N; G.w = 2*pi/N*ones(N,1);
end
s = G.s; G.x = G.Z(s); G.sp = abs(G.Zp(s)); G.nx = -1i*G.Zp(s)./G.sp;
G.cur = -real(conj(G.Zpp(s)).*G.nx)./G.sp.^2;  % curvature
