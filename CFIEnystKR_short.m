function a = CFIEnystKR(G,i,j,k,eta)
% CFIENYSTKR - ij element of 2D Helmholtz CFIE Nystrom matrix, Kapur-Rokhlin corr
g6 = [4.967362978287758 -16.20501504859126 25.85153761832639 ...
      -22.22599466791883 9.930104998037539 -1.817995878141594]; % 6th order
sw = G.sp(j)*G.w(j);                          % speed weight
N = numel(G.x); l = mod(i-j,N); if l>N/2, l=N-l; end   % index distance i to j
if l>0 & l<=6, sw = sw * (1 + g6(l)); end     % apply correction to sw
if i==j, a = 0; return; end                   % kill diagonal
d = G.x(i)-G.x(j); kr = k*abs(d);             % last 3 lines do CFIE kernel:
costheta = real(conj(G.nx(j)).*d)./abs(d);    % theta angle between x-y & ny
a = (1i/4) * (k*costheta*besselh(1,kr) - 1i*eta*besselh(0,kr)) * sw;
