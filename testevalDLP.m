% test the DLP Laplace 2D evaluator evalDLP. Barnett 6/7/14

clear
G = smoothstar(.3,5);                 % curve shape
N = 150; G = curvquad(G, 'ptr', N);   % curve quadrature rule
g = -1.4:0.01:1.4; [xx yy] = meshgrid(g); t = xx + 1i*yy; % targets
sigma = ones(N,1);                                        % density
u = evalDLP_mine(t,G,sigma);               % use the evaluator you will write
figure; imagesc(g,g,log10(abs(u+1))); caxis([-16 0]); colorbar;
hold on; plot([G.x;G.x(1)], '-'); plot(G.x, 'k.','markersize',5); % show geom
axis xy equal tight; title('interior Laplace DLP log_{10} error');
set(gcf,'paperposition',[0 0 4 3]);
print -depsc2 ../testevalDLP.eps
