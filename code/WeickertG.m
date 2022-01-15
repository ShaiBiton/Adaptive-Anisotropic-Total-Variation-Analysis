function [res] = WeickertG(x,K)
    C = 3.31488; m = 4;
    if (x <= 0)
        res = 1;
    else
        res = 1 - exp(-C./(x./K).^m);
    end
end

