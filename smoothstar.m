function G = smoothstar(a,w)
% SMOOTHSTAR - set up [0,2pi) parametrization of smooth star closed curve
R = @(t) 1 + a*cos(w*t); Rp = @(t) -w*a*sin(w*t); Rpp = @(t) -w*w*a*cos(w*t);
G.Z = @(t) R(t).*exp(1i*t); G.Zp = @(t) (Rp(t) + 1i*R(t)).*exp(1i*t);
G.Zpp = @(t) (Rpp(t) + 2i*Rp(t) - R(t)).*exp(1i*t);
