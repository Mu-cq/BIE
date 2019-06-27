%GreensFormula Check Complex 
%Checking Green's Formula:

%The known harmonic function
u = @(x,y) 4*(x.^2 - y.^2) + 8*x.*y;
ugrad = @(x,y) [8*x + 8*y, 8*x-8*y];

%The general solution to Laplace's equation , x,y E R2
I = @(x,y) 1/(2*pi)*log(1./sqrt(sum((x-y).^2)));
Igrad = @(x,y) 1/(2*pi).*(x-y)./sum((x-y).^2);

N = 100;
a = 4;
omega=0.4;
%Our Region of integration 
%G = smoothstar(a,omega);
%Prepare the quadrature scheme
%G = curvquad(G, 'ptr', N);
%G = curvquad(G, 'panel', N, 10);

G = curveCircle(2,100);

%target point
x = [0.3, 0.1];  %inside the region R 

%wt_x = 10/length(G.s); %not sure about this  
%wt_y = 8/length(G.s);
%wt = sqrt(sum(wt_x^2 + wt_y^2));
%wt = sqrt((G.w/(2*pi)*8 + G.w/(2*pi)*10).^2); 
result = 0;

%integrate over source points
for i=1:1:length(G.s)
   
   %unparametrize x = Re(Z), y = Im(Z)
   %s_x = (1+a*cos(omega*G.s(i)))*cos(G.s(i));
   %s_y = (1+a*cos(omega*G.s(i)))*sin(G.s(i));
   
   %s_x = real(G.Z(i));
   %s_y = imag(G.Z(i));
    s_x = G.X(1,i);
    s_y = G.X(2,i);
   
   %fprintf("s_x: %f   ", s_x);
   %fprintf("s_y: %f   \n", s_y);
   
   if(and(s_x == x(1),s_y == x(2)))
       integrand = -G.cur(x)/(4*pi)*wt(i);
   else
    %n_x = real(G.nx(i));
    %n_y = imag(G.nx(i));
     n_x = G.dydx(1,i);
     n_y = G.dydx(2,i);
    dUdn = dot(G.dydx(:,i), ugrad(s_x,s_y));
    dIdn = dot(G.dydx(:,i), Igrad(x, G.X(:,i)'));
    integrand = (dUdn.*I(x,G.X(:,i)') - u(x(1),x(2)).*dIdn)*wt(i); 
    end       
  result = result + integrand;
end

fprintf("u(x) = %f\n", u(x(1),x(2)));
fprintf("Estimate(x) = %f\n", result);