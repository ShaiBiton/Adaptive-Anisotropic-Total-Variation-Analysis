function divP = divA(A, P)
    Px = P(:,:,1); Py = P(:,:,2);
    
    APx = A{1,1}.*Px + A{2,1}.*Py;
    APy = A{1,2}.*Px + A{2,2}.*Py;
    
	fx = APx-APx(:,[1 1:end-1]); fx(:,1) = APx(:,1); fx(:,end) = -APx(:,end-1);
	fy = APy-APy([1 1:end-1],:); fy(1,:) = APy(1,:); fy(end,:) = -APy(end-1,:);
    divP = fx+fy;
end