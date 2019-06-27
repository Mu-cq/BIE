%Exercise 8
%n_y = 
%k = @(x,y) dot(n_y,(x-y))/(2*pi*abs(x-y).^2); %

fund_sol = @(x,y) (1/2*pi) * log(1/abs(x-y));
g = -1.4:0.01:1.4; [xx yy] = meshgrid(g);
Z = (1/2*pi)*log(1./abs(xx-yy));
plot3(xx,yy,Z)
