%Checking Green's Formula:

%The known harmonic function
u = @(x,y) 4*(x.^2 - y.^2) + 8*x.*y;
ugrad = @(x,y) [8*x + 8*y, 8*x-8*y];

%The general solution to Laplace's equation , x,y E R2
I = @(x,y) 1/(2*pi)*log(1./sqrt(sum((x-y).^2)));
Igrad = @(x,y) 1/(2*pi)*(x-y)./sum((x-y).^2);

N = 500;

%Our Region of integration
G = curveCircle(3,N);

%target point

tx = [-4.0:0.2:4.0];
ty = [-4.0:0.2:4.0];
t = [tx;ty]; 

x = [1.5, 1.1];  %inside the region R 

result = 0;

%%
%integrate over source points
for i=1:1:length(G.s)

    s_x = G.X(1,i);
    s_y = G.X(2,i);

   if(and(s_x == x(1),s_y == x(2)))
       integrand = -G.cur(x)/(4*pi)*wt(i)*G.spd(i);
   else
    dUdn = dot(G.n(:,i), ugrad(s_x,s_y)); 
    dIdn = dot(G.n(:,i), Igrad(x, G.X(:,i)'));
    integrand = (dUdn*I(x,G.X(:,i)') - u(G.X(1,i),G.X(2,i))*dIdn)*G.wts(i)*G.spd(i); 
   
   end
   result = result + integrand;
end

fprintf("u(x) = %f\n", u(x(1),x(2)));
fprintf("Estimate(x) = %f\n", result);