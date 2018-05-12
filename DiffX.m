classdef DiffX <  matlab.System
    
    properties
       H1; H2; H3; H4;                 % Second derivative equasions 
    end
    
    
    methods(Access = public)
       
        function giveH(obj, H1_, H2_, H3_, H4_)
            obj.H1 = H1_;
            obj.H2 = H2_;
            obj.H3 = H3_;
            obj.H4 = H4_;
        end
        
        
       function [dX] = diffX(obj, t, X)
           % % returns the derivative of statevector X at time t
           dX = zeros(8,1);
           dX(1) = X(5);
           dX(2) = X(6);
           dX(3) = X(7);
           dX(4) = X(8);
           dX(5) = obj.H1(X(1), X(2), X(3), X(4), X(5), X(6), X(7), X(8));
           dX(6) = obj.H2(X(1), X(2), X(3), X(4), X(5), X(6), X(7), X(8));
           dX(7) = obj.H3(X(1), X(2), X(3), X(4), X(5), X(6), X(7), X(8));
           dX(8) = obj.H4(X(1), X(2), X(3), X(4), X(5), X(6), X(7), X(8));

    end    
 
        
    end
end
