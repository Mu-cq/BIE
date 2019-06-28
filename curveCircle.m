%%
%parametric/polar equation for a circle
%and an associated quadrature scheme 
function C = curveCircle(r, N)

    C.Z = @(theta) [r.*cos(theta); r.*sin(theta)];
    C.Zp = @(theta) [-1*r.*sin(theta); r.*cos(theta)];	
    C.Zpp = @(theta) [-1*r.*cos(theta) ; -1*r.*sin(theta)] ;
 
   %setting up "PTR" Quadrature 

    pts = (2*pi*(1:N)'/N); %equally spaced pts on [0,2pi]
    wts = 2*pi/(N)*ones(N,1); 
  
%     %setting up a "gauss-legendre" quadrature - NOT finished 
%     C.p = 16; %points per panel
%     C.n_panels = ceil(N/C.p);    
% 
%     %divide [0,2pi] into n_panels 
% 
%     [pts_g, wts_g] = gauss(C.p); %intially on interval[-1,1]
%     pts_g = (pts_g+1)/2; %shift to [0,1]
%     wts_g = wts_g/2;  %why? n_pts doesn't change 
%     
%     %wts repeat for each panel
%     wts = repmat(wts_g', [1, C.n_panels]);
%     
%     %pts for each panel 
%     C.left_endpts = [0:C.n_panels-1]*2*pi/C.n_panels;
%     C.right_endpts = [1:C.n_panels]*2*pi/C.n_panels;
%     
%     %expand endspts to a vector containing per point displacement 
%     C.ptsDisplacement_left = reshape(repmat(C.left_endpts',[1, C.p])', 1, []);
%     C.ptsDisplacement_right = reshape(repmat(C.right_endpts',[1, C.p])', 1, []);  
%     
%     %repeat pts for each panel, but shift to start at left endpoint 
%     C.pts_repeat = reshape(repmat(pts_g, [1, C.n_panels]), 1 , []);
%     C.pts_shifted = C.pts_repeat + C.ptsDisplacement_left; 
% 
%     C.half_lengths = (C.ptsDisplacement_right-C.ptsDisplacement_left)/2;
%     C.midpts = (C.ptsDisplacement_left+C.ptsDisplacement_right)/2; 
% 
%     C.s = (C.half_lengths.*C.pts_shifted + C.midpts);
    
    C.s = pts';
    C.X = C.Z(C.s);
    C.wts = wts;
    
    C.Zp_eval = C.Zp(C.s);

    C.spd = sqrt(sum(C.Zp_eval.^2)); %r
    
    CN = C.Zp_eval./C.spd; 
    
    C.n = [CN(2,:); -1*CN(1,:)];

    C.curv = -1*C.Zpp(C.s).*C.n./C.spd.^2;  % curvatures

%    C.curr = @(x) -1*C.Zpp(x).*C.n./C.spd.^2;

    
    
