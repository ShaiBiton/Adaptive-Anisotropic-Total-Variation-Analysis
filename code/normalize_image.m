function image = normalize_image(Im)
% Normalize image
image = double(Im);
image = image-min(image(:));
image = image/max(image(:));
image = image - mean(image(:));
end

