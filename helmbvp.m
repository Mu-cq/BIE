% Demo code for BIE CFIE solution of 2D Helmholtz BVP. Barnett 6/8/14

G = smoothstar(.3,5); N = 500; G = curvquad(G,'ptr',N);      % set up curve
k = 10; eta = k;                           % wavenumber, SLP mixing amount
f = @(z) besselh(0,k*abs(z-(0.2+0.3i)));   % known soln: interior source
rhs = 2*f(G.x);                            % bdry data at nodes
if 0
A = nan(N,N); for i=1:N, for j=1:N, % bare nystrom for testing
    if i~=j, d = G.x(i)-G.x(j); kr = k*abs(d);           % CFIE kernel...
      costheta = real(conj(G.nx(j)).*d)./abs(d);  % theta angle between x-y & ny
      A(i,j) = (1i/4) * (k*costheta*besselh(1,kr) - 1i*eta*besselh(0,kr));
    else, A(i,j) = -G.cur(j)/(4*pi); end
    A(i,j) = 2* A(i,j) * G.sp(j) * G.w(j); % speed wei and make 2(D-i.eta.S)
  end, end
else
A = nan(N,N); for i=1:N, for j=1:N, A(i,j) = 2*CFIEnystKR(G,i,j,k,eta);end,end
end
sigma = (eye(N) + A) \ rhs;        % dense solve
g = -1.6:0.02:1.6; [xx yy] = meshgrid(g); t = xx + 1i*yy;    % targets for plot
'eval', u = evalCFIEhelm(t,G,sigma,k,eta);            % evaluate soln

figure; imagesc(g,g,log10(abs(u-f(t)))); caxis([-16 0]); colorbar;
hold on; plot([G.x;G.x(1)], '-'); plot(G.x, 'k.','markersize',5); % show geom
axis xy equal tight; title('ext Dirichlet 2D Helm BVP log_{10} error')
%set(gcf,'paperposition',[0 0 4 3]); print -depsc2 ../helmbvperr.eps

figure; imagesc(g,g,real(u)); caxis(.4*[-1 1]); colorbar;
hold on; plot([G.x;G.x(1)], '-'); plot(G.x, 'k.','markersize',5); % show geom
axis xy equal tight; title('ext Dirichlet 2D Helm BVP Re u')
%set(gcf,'paperposition',[0 0 4 3]); print -depsc2 ../helmbvp.eps
