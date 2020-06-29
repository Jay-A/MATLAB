function [Cfs] = Jacobi_Coeffs(n, alpha, beta)
	
	%if n<2
	%	n=2;
	%end
	
	Cfs = zeros(n+1, n+1);
	Cfs(1,n+1) = 1;
	
	switch n
		case 0
		case 1
			Cfs(2,n+1) = (1/2)*(2*(alpha+1)-(alpha+beta+2));			
			Cfs(2,n) = (1/2)*(alpha+beta+2);
		otherwise		
			Cfs(2,n+1) = (1/2)*(2*(alpha+1)-(alpha+beta+2));			
			Cfs(2,n) = (1/2)*(alpha+beta+2);
				
			c = @(n) n+alpha+beta; 
			for counter = 3:n+1	
				k=counter-1;
				temp_poly = zeros(1, n+1);
				temp_poly(1,n:n+1) = [c(2*k-1)*c(2*k-2)*c(2*k), c(2*k-1)*(alpha^2-beta^2)];

                tempMatlab = conv(temp_poly, Cfs(counter-1,:));
				temp = (1/(2*k*c(k)*c(2*k-2))) .* ...
                    (tempMatlab(1,n+1:2*n+1) - 2 .* (k-1+alpha) .* (k-1+beta) .* c(2*k) .* Cfs(counter-2,:));
				length(temp);
				size(Cfs);
				Cfs(counter, :) = temp(1,:);
			end
	end
	Cfs = Cfs(n+1,:);
			