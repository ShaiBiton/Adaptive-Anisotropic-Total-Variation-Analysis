function dP = gradA(A, P)
	fx = P(:,[2:end end])-P;
	fy = P([2:end end],:)-P;
    
    Afx = A{1,1}.*fx + A{1,2}.*fy;
    Afy = A{2,1}.*fx + A{2,2}.*fy;
    
    dP = cat(3,Afx,Afy);
end