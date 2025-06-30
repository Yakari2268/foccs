function [x,t] = fde_sol(f,q,t_span,ne,N,x0)
    %pre-allocate memory 
    x = zeros(N+1,ne);
    
    x(1,:) = x0(1);
    aj_n+1 = 1
    bj_n+1 = 1
    for m=1:N+1
        for n = 1:ne
            
            %predictor
            x(m+1,n) = x0(n) + (1/gamma(q(n)))*(summation_1);
            
            %corrector
            x(m+1,n) = x0(n) + (h^q(n)/gamma(q(n)+2))*f(x(m+1,n),t(m+1)) + (h^q/(gamma(q(n)+2))*summation_2;
        end
    end

end