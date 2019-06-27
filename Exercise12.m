%Exercise 12

% Demo code for BIE solution of 2D Laplace BVP (PTR/panels). Barnett 6/7/14

G = smoothstar(.3,5);
%N = 200; G = curvquad(G,'panel',N,16); N = numel(G.x);  % ...or panel case

%f = @(z) real(exp(z));             % known Laplace soln in 2D (z is complex)

f = @(z) log(abs(z-(1+1.5i))); 

Ns = 20:20:200;
maxSolutionError=zeros(1,length(Ns));
k  = 1;

for N = Ns
    G = curvquad(G,'ptr',N);                       % PTR case
    rhs = -2*f(G.x);                   % bdry data at nodes
    A = nan(N,N); 
    for i=1:N
        for j=1:N 
            A(i,j) = -2*Dnystmatel(G,i,j); 
        end
    end

    sigma = (eye(N) + A) \ rhs;  % dense solve
    g = -5.4:0.01:-3.4; [xx yy] = meshgrid(g); t = xx + 1i*yy;    % targets for plot
    u = evalDLP(t,G,sigma);            % evaluate soln
    %u = reshape(lapDevalclose(t,G,sigma,'i'),size(t));  % accurate up to boundary
    
    maxSolutionError(k) = max(log10(abs(u-f(t))),[], 'all');
    fprintf("%d\n", maxSolutionError(k));
    k = k + 1;
    
%     figure; imagesc(g,g,log10(abs(u-f(t)))); caxis([-16 0]); colorbar;
%     hold on; plot([G.x;G.x(1)], '-'); plot(G.x, 'k.','markersize',5); % show geom
%     axis xy equal tight; title('interior Dirichlet 2D Laplace BVP log_{10} error');
%     set(gcf,'paperposition',[0 0 4 3]);
%     
    
    %print -depsc2 ../lapbvperr.eps
end
figure;
plot(Ns,maxSolutionError);
title('Maximum Solution Error');
