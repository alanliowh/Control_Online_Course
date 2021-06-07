function [out]= f_switch(x0,x1,x)
    
    if x < x0
        out = 0;
    elseif ((x0<= x) && (x<= x1))
        out = (x-x0)/(x1-x0);
    elseif x>x1
        out = 1;
    end
end
        