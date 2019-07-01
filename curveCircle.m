%%
%parametric/polar equation for a circle
%and an associated quadrature scheme 
function C = curveCircle(r, N)

    C.Z = @(theta) [r.*cos(theta); r.*sin(theta)];
    C.Zp = @(theta) [-1*r.*sin(theta); r.*cos(theta)];	
    C.Zpp = @(theta) [-1*r.*cos(theta) ; -1*r.*sin(theta)] ;
 
   %setting up "PTR" Quadrature 

  %  pts = (2*pi*(1:N)'/N); %equally spaced pts on [0,2pi]
  %  wts = 2*pi/(N)*ones(N,1); 
  %  C.s = pts';
    
    %setting up a "gauss-legendre" quadrature - NOT finished 
    C.p = 16; %points per panel
    C.n_panels = ceil(N/C.p);    

    %divide [0,2pi] into n_panels 

    [pts_g, wts_g] = gauss(C.p); %intially on interval[-1,1]
    pts_g = (pts_g+1)/2; %shift to [0,1]
    wts_g = wts_g/2;
    
    %wts repeat for each panel
    wts = repmat(wts_g*2*pi/C.n_panels, [C.n_panels,1]);
    
    %pts for each panel 
    C.endpts = [0:C.n_panels]*2*pi/C.n_panels;

    pts = nan(N, 1);
    for i=1:C.n_panels
    	pts((i-1)*C.p + (1:C.p)) = C.endpts(i) + (C.endpts(i+1)-C.endpts(i))*pts_g;
    end    

    C.s = pts';
   
    C.X = C.Z(C.s); 
    C.wts = wts';
    
    C.Zp_eval = C.Zp(C.s);

    C.spd = sqrt(sum(C.Zp_eval.^2)); %r
    
    CN = C.Zp_eval./C.spd; 
    
    C.n = [CN(2,:); -1*CN(1,:)];

    C.curv = -1*C.Zpp(C.s).*C.n./C.spd.^2;  % curvatures

%    C.curr = @(x) -1*C.Zpp(x).*C.n./C.spd.^2;

    
    
