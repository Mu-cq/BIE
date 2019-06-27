
function result = evaluateDLP2d(sigma,G,T)
    %t are 2 dimensional target points
    %G is a curve with embedded quadrature scheme
    %sigma is a density function

    n_srcs = length(G.X);

    result = 0*length(T);

    %loop over the target points

  for j =1:1:length(T)    %loop over the quadrature points
    res = 0;
    for i=1:1:n_srcs
           num = dot(G.n(:,i), T(:,j)-G.X(:,i));
           magdist = sum((T(:,j)-G.X(:,i)).^2);
           integrand = num/magdist;
           res = res + integrand*G.wts(i)*sigma(G.X(2,i))*G.spd(i); %spd(r)=r 
     end
     res = res/(2*pi);
     result(j) = res;
  end
    for k=1:1:length(T)
        fprintf("u(%f,%f) = %f\n", T(1,k), T(2,k), result(k));
    end
    end
