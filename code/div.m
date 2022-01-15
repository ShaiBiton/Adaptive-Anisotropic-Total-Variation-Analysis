function divP = div(P)
    Px = P(:,:,1); Py = P(:,:,2);
	fx = Px-Px(:,[1 1:end-1]); fx(:,1) = Px(:,1); fx(:,end) = -Px(:,end-1);
	fy = Py-Py([1 1:end-1],:); fy(1,:) = Py(1,:); fy(end,:) = -Py(end-1,:);
    divP = fx+fy;
end