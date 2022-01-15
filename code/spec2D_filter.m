function f_H = spec2D_filter( Phi, H, f_r, dt )
% private function by Guy Gilboa (Jan 2015)
% Construct filtered image from Phi, filter H and residual image f_r
% Example: f_H = specTV_filter( Phi, H, f_r, dt )

if (H(end) == 1  )
    f_H = f_r;  % include residual
else
    f_H = zeros(size(f_r));  % does not include residual
end

for i=1:length(H)
    f_H = f_H+H(i)*Phi(:,:,i)*dt;
end % for i

end

