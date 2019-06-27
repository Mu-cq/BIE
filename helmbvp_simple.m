G = smoothstar(.3,5); N = 500; G = curvquad(G,'ptr',N);      % set up curve
k = 10; eta = k;                           % wavenumber, SLP mixing amount
f = @(z) besselh(0,k*abs(z-(0.2+0.3i)));   % known soln: interior source
rhs = 2*f(G.x);                            % bdry data at nodes
A = nan(N,N); for i=1:N, for j=1:N, A(i,j) = 2*CFIEnystKR(G,i,j,k,eta); end, end
sigma = (eye(N) + A) \ rhs;                % dense solve
g = -1.6:0.02:1.6; [xx yy] = meshgrid(g); t = xx + 1i*yy;    % targets for plot
u = evalCFIEhelm(t,G,sigma,k,eta);         % evaluate soln
figure; imagesc(g,g,log10(abs(u-f(t)))); caxis([-16 0]); colorbar; axis xy equal

