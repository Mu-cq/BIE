%Exercise 2
%Exercise_2 
%Template code for convergence of PTR. Also serves as Matlab warm up.
% Barnett 6/5/14

%Regarding Exercise 3: does the target point x correspond to a? And 
%Does a approach the unit circle as it gets closer to 1??
          

%a = 0.1; %as a gets smaller, alpha gets smaller, 
          %and the convergence slows
          
          %f is vertically stretched / the parabola narrows 
          
for a=1:-0.2:0.3           

    f = @(s) a./(a^2+sin(s/2).^2);  % note elementwise ops
    s=0:1e-3:2*pi; figure; 
    plot(s,f(s),'-'); xlabel('s'); title('integrand');
    Ns = 10:10:100;   % convergence study.
    Is = 0*Ns;   % will get filled by you with each I_N looping over N's

    for i=1:numel(Ns), N = Ns(i);
      % ... insert quadrature here for Is(i) = I_N.
      w = 2*pi/N;
      j = 1:N;
      sj = 2*pi*j/N ; 
      Is(i) = sum(w*f(sj));
      fprintf('N=%d\t I_N = %.16g\n',N,Is(i))
    end

    %Taking convergence value of max N to be real value
    realValue = Is(numel(Ns));

    %Compute error for each Is(i)
    err = abs(Is-realValue);

    alpha = 2*imag(asin(1i*a));

    %By Davis Theorem, convergence  = 
    convergence = exp(-alpha*Ns) ; 

    %Convergence type?
    figure;
    semilogy(Ns, err, 'b'); hold on;
    plot(Ns, convergence, 'r-'); 
    xlabel('Ns'); title('Convergence Rate'); ylabel('Quadrature Error');
    v = axis; v(3:4)=[1e-16 1]; axis(v);

    
end
