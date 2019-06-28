%Checking Green's Formula:

%The known harmonic function 
u = @(x,y) 4*(x.^2 - y.^2) + 8*x.*y;
ugrad = @(x,y) [8*x + 8*y; 8*x-8*y];



%The general solution to Laplace's equation , x,y E R2
I = @(x,y) 1/(2*pi)*log(1./sqrt(sum((x-y).^2)));
%Igrad = @(x,y) 1/(2*pi)*(x-y)./sum((x-y).^2);

N = 500;

%Our Region of integration G
G = curveCircle(3,N);

%target points over a whole grid, inside and outside of region G
[xx,yy] = meshgrid(-3:0.2:3.0);

%condense the points to a [2 row, n_targetPt col] matrix of (x,y) coordinates 
t = [reshape(xx, [1,length(xx)*length(xx)]); reshape(yy, [1,length(yy)*length(yy)])];

res = 0*length(t);

%Estimate the integral defined by Green's Formula for each target point
for j=1:1:length(t)
    
    dUdn = dot(G.n, ugrad(G.X(1,:),G.X(2,:))); 
    dIdn = dot(G.n, Igrad(t(:,j), G.X));
    integrand = sum((dUdn.*I(t(:,j),G.X) - u(G.X(1,:),G.X(2,:)).*dIdn).*(G.wts).*(G.spd),2); 
    res(j) = integrand; 
end


%Compare estimate with u(t) for every t
uOfT=u(t(1,:),t(2,:));
err = uOfT-res;

%plot error
figure;
imagesc(t(1,:),t(2,:),log10(abs(err)));  caxis([-16 0]);  colorbar; hold on;

%plot the region of integration 
plot(G.X(1,:), G.X(2,:), 'r-', 'markersize',16); ylim([-3,3]);
title('Greens Formula log_{10} Error');  colorbar; xlabel('X'); ylabel('Y');




