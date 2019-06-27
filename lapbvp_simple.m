G = smoothstar(.3,5); N = 150; G = curvquad(G,'ptr',N);      % set up curve
f = @(z) log(abs(z-(1+1.5i)));     % known Laplace soln in 2D (z is complex)
rhs = -2*f(G.x);                   % bdry data at nodes
A = nan(N,N); for i=1:N, for j=1:N, A(i,j) = -2*Dnystmatel(G,i,j); end, end
sigma = (eye(N) + A) \ rhs;        % dense solve
g = -1.4:0.01:1.4; [xx yy] = meshgrid(g); t = xx + 1i*yy;    % targets for plot
u = evalDLP(t,G,sigma);            % evaluate soln
figure; imagesc(g,g,log10(abs(u-f(t)))); caxis([-16 0]); colorbar; axis xy equal
