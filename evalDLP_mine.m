%%
function u = evalDLP_mine(t, G,sigma)

%EVALDP: a function that uses the quadrature rule in G to evaluate the
%interior Lapalace Problem 


%size(s) rows, n_target points columns

% %this only works when targets is a 1 dimensional array
% stackedNorms = repmat(G.nx, size(targets));
% quot = stackedNorms.*(targets-G.x)./(2*pi*(targets-G.x).^2);
% stackedWeights = repmat(G.w,size(targets));
% stackedSigma = repmat(arrayfun(sigma,G.s),size(targets));
% stackedSp = repmat(G.sp, size(targets));
% res = sum(stackedWeights.*(quot).*stackedSigma.*stackedSp,1); 
% % u = res;

%%
%for t = 2d (x,y)


N = numel(G.x);
u = 0*t; 
for j=1:N
    d = t-G.x(j);
    numerator = real(conj(G.nx(j)).*d); 
    denominator = (abs(d).^2);
    kernel = (numerator ./ denominator);
    prod = G.w(j)*sigma(j)*G.sp(j)*kernel;
    u = u+prod;
end
u = u/(2*pi);     

