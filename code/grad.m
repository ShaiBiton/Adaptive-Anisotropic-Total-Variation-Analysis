function dP = grad(P)
	fx = P(:,[2:end end])-P;
	fy = P([2:end end],:)-P;
    dP = cat(3,fx,fy);
end