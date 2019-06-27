function [u ux uy do] = lapDevalclose(x, s, tau, side) 
% LAPDEVALCLOSE - evaluate potential and field for DLP on global quadr curve
%
% u = lapDevalclose(x,s,tau,side) returns potentials at targets x due to DLP
%  with density tau living on curve s with global periodic trapezoid rule.
%  side controls whether targets are inside (default) or outside.
%  A scheme based on Helsing/Ioakimidis globally compensated quadr is used.
%
% Our definition of the DLP on curve \Gamma is, associating R^2 with the complex
% plane, with x = (x_1,x_2) a point in the complex plane,
%
%    u(x) = Re (-1/2\pi i) \int_\Gamma \tau(y) / (y-x) dy
%
% [u ux uy] = lapDevalclose(x,s,tau,side) also returns field (the two first
%  partial derivatives of potential u)
%
% inputs:
% x = M-by-1 list of targets in complex plane
% s = curve object containing s.x source pts, s.nx outward unit normals, s.cur 
%     curvatures, s.w weights (2pi/N), s.sp speeds.
%     s.a interior point far from bdry (only used for exterior case).
% tau = double-layer density values at nodes
% side = 'i','e' to indicate targets are interior or exterior to the curve.
%
% outputs:
% u     = potential values at targets x (M-by-1)
% ux,uy = (optional) first partials of potential at targets
% do    = (optional) diagnostic outputs:
%   do.vb : complex values on bdry (ie, v^+ for side='e' or v^- for side='i')
%
% complexity O(N^2) for evaluation of v^+ or v^-, plus O(NM) for globally
% compensated quadrature to targets
%
% Barnett 10/9/13. Switched to cauchycompeval 10/23/13.
% Edit CBMS (make s.w not speedweights, just weights) 6/7/14

if nargin<4, side='i'; end
wantder = nargout>1;

% Helsing step 1: eval bdry limits at nodes of v = complex DLP(tau)...
vb = 0*s.x;                % will become v^+ or v^- bdry data of holom func
taup = perispecdiff(tau);  % parameter-deriv of tau
cw = 1i*s.nx.*s.w.*s.sp;   % complex speed weights
N = numel(s.x);
for i=1:N, j = [1:i-1, i+1:N];   % skip pt i
  vb(i) = sum((tau(j)-tau(i))./(s.x(j)-s.x(i)).*cw(j)) + taup(i)*s.w(i);
end
vb = vb*(1/(-2i*pi));   % prefactor
if side=='i', vb = vb - tau; end % JR's add for v^-, cancel for v^+
do.vb = vb;

% Helsing step 2: compensated close-evaluation of u = Re(v) & its deriv...
if wantder, [v vp] = cauchycompeval(x,s,vb,side);
  ux = real(vp); uy = -imag(vp);
else v = cauchycompeval(x,s,vb,side); end
u = real(v);

function fp = perispecdiff(f) % spectral differentiation sampled periodic func
N = numel(f);          % can be row or col vec
fh = fft(f);
if mod(N,2)==0, fp = ifft(1i*[0:N/2-1, 0, -N/2+1:-1]' .* fh(:));  % even case
else, fp = ifft(1i*[0:(N-1)/2, -(N-1/2):-1]' .* fh(:)); end       % odd
fp = reshape(fp, size(f));

function testperispecdiff
N = 50; tj = 2*pi/N*(1:N)';
f = sin(3*tj); fp = 3*cos(3*tj);   % trial periodic function & its deriv
norm(fp-perispecdiff(f))
