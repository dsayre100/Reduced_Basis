function data = exact_soln

data = struct('f',@f,'exactu',@exactu);

% exact solution u
    function u =  exactu(p,s)        
        u = (8*pi^2)^(-s)*sin(2*pi*p(:,1)).*sin(2*pi*p(:,2));
    end

% right hand side
    function f =  f(p)        
        f = 2*(sin(2*pi*p(:,1)).*sin(2*pi*p(:,2)));
    end
end