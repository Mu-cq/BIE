%%
%curveCircle a function defined on a circle of radius r
function C = curveCircle(r, N)

    C.Z = @(theta) [r.*cos(theta); r.*sin(theta)];
    C.Zp = @(theta) [-1*r.*sin(theta); r.*cos(theta)];	
    C.Zpp = @(theta) [-1*r.*cos(theta) ; -1*r.*sin(theta)] ;
   %setting up "PTR" Quadrature 

    pts = (2*pi*(1:N)'/N); %equally spaced pts on [0,2pi]
    wts = 2*pi/(N)*ones(N,1); 
    
    C.s = pts';
    C.X = C.Z(pts');
    C.wts = wts;
    
    C.Zp_eval = C.Zp(C.s);

    C.spd = sqrt(sum(C.Zp_eval.^2)); %r
    
    CN = C.Zp_eval./C.spd; 
    
    C.n = [CN(2,:); -1*CN(1,:)];

    C.curv = -1*C.Zpp(C.s).*C.n./C.spd.^2;  % curvatures
 

    
    
