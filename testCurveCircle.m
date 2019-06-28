
t=0:0.2:2*pi;
%Example of what a successful parameterization looks like
%G=smoothstar(0.3,6);
%plot(real(G.Z(t)),imag(G.Z(t))); hold on;
%G = curvquad(G, 'panel', 100, 16);
%plot(real(G.x),imag(G.x),'o'); hold on;
%quiver(real(G.x),imag(G.x), real(G.nx),imag(G.nx))


%Test curve Circle
figure;
C = curveCircle(1,100);
c = C.Z(t);
plot(c(1,:),c(2,:)); hold on;
plot(C.X(1,:), C.X(2,:),'o'); hold on;
quiver(C.X(1,:),C.X(2,:),C.n(1,:),C.n(2,:));     %(dy/dx = (dy/dtheta)/(dx/dtheta)   

%if F is holomorphic, then the integral or F around a closed boundary is 0.

%strange: for even n, zero error after 2. For odd- fun oscillatory error 
Ns = 20:10:500; 
f =  @(z) sin(z); %Z is complex

err = 0*length(Ns);
integral = 0*length(Ns);
j = 1;
for n=Ns
    d = curveCircle(3,n); %define a quadrature scheme with n points
    z = d.X(1,:) + 1i*d.X(2,:);    %convert to complex
    integral(j) = sum(f(z).*d.wts.*d.spd);   %compute the integration estimation
    err(j) = abs(integral(j)); %(estimate - 0);
    j = j + 1;
end 

alpha_obs = (log(err(end))-log(err(1)))/(Ns(end)-Ns(1));
figure;
semilogy(Ns, err, 'r+-');
hold on;
plot(Ns, exp(-alpha_obs*Ns),'b+'); %convergence rate
xlabel('Ns'); title('Convergence Rate'); ylabel('Quadrature Error');
