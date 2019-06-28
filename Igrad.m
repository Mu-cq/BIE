
function ans = Igrad(x,y)
	 if(and(x(1) == y(1), x(2) == y(2)))
	     ans = G.curr(x);
	     fprintf("HERE");
         else
	     ans =  1/(2*pi)*(x-y)./sum((x-y).^2);
	 end
end		
