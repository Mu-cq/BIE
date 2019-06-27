% make panels.eps and curvenodes.eps figures. Barnett 6/5/14
p=16;
n=5;
L = 2*pi;
l=0.1;
[x w] = gauss(p); x1 = (x+1)/2; w1 = w/2; % put on [0,1]
se = (0:n)/n*L; % panel ends
w = repmat(w1, [1 n]);  % stack copies of weights together
s = nan(n*p,1); for i=1:n, s((i-1)*p + (1:p)) = se(i)+(se(i+1)-se(i))*x1; end
% now s and w are set up.

figure; plot([0 L],[0 0],'-'); hold on; plot(s,0*s,'k.','markersize',10);
for i=1:n+1, plot([se(i) se(i)], [-l l],'k-'); end
axis equal off;
print -depsc2 ../panels.eps

G = smoothstar(.3,5); % curve
G = curvquad(G,'ptr',100);
figure; plot([G.x;G.x(1)], '-'); hold on; plot(G.x, 'k.','markersize',10);
axis equal off; print -depsc2 ../curvenodes.eps
