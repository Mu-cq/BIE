%Checking Green's Formula:

%The known harmonic function 
u = @(x,y) 4*(x.^2 - y.^2) + 8*x.*y;
ugrad = @(x,y) [8*x + 8*y; 8*x-8*y];

%The general solution to Laplace's equation , x,y E R2
I = @(x,y) 1/(2*pi)*log(1./sqrt(sum((x-y).^2)));
Igrad = @(x,y) 1/(2*pi)*(x-y)./sum((x-y).^2);

N = 500;

%Our Region of integration
G = curveCircle(3,N);

%target point

tx = [-3.0:0.2:3.0]; %only within the region 
ty = [-3.0:0.2:3.0];
t = [tx;ty]; 

%TO DO- get target points over a whole grid!!
%[xx,yy] = meshgrid(-3:0.2:3.0);
%t = xx;

x = [1.5; 1.1];  %inside the region R 
res = 0*length(t);

for j=1:1:length(t)
    dUdn = dot(G.n, ugrad(G.X(1,:),G.X(2,:))); 
    dIdn = dot(G.n, Igrad(t(:,j), G.X));
    integrand = sum((dUdn.*I(t(:,j),G.X) - u(G.X(1,:),G.X(2,:)).*dIdn).*(G.wts').*(G.spd),2); 
    res(j) = integrand; 
end



uOfT=u(t(1,:),t(2,:));
err = uOfT-res;

figure;
imagesc(tx,ty,log10(abs(err)));  caxis([-16 0]);  colorbar; hold on;
plot(G.X(1,:), G.X(2,:), 'r-', 'markersize',8); ylim([-3,3]);

