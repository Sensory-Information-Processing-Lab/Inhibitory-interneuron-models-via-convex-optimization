function masked_stimulus = gen_masked_stim(stimulus,radius,center,ph)
% stimulus: 6-D. contrast, orientation, freq,phase, the image.
% Generate masked stimulus at a certain position from the original
% full-sized stimulus. 
% radius: the range of radius to be generated.
% ph:
%   r_pixel: the radius of a pixel
% output:
%   index in the sequence: contrast, orientation, freqs, phase,
%   radius, x, y.

% Mengchen Zhu.

% The size of the image
img_sz = size(stimulus, 5);

%Add r dimension
masked_stimulus = repmat(stimulus, [1,1,1,1,1,1,length(radius)]); 
masked_stimulus = permute(masked_stimulus, [1 2 3 4 7 5 6]);

switch ph.sample
  case 'over'
    max_intersection = pi * ph.r_pixel^2;
    for x = 1: img_sz
        for y = 1: img_sz
            for z = 1:size(radius,2)
                if ((ph.r_pixel + radius(z)) < sqrt((img_sz/2-y)^2 + ...
                                                 (img_sz/2 -x)^2))
                    % If no intersection
                    masked_stimulus(:,:,:,:,z,x,y) = 0;
                elseif (((ph.r_pixel + radius(z)) >= sqrt((img_sz/2-y)^2 + ...
                                                       (img_sz/2 -x)^2)) ...
                        && ((radius(z)-ph.r_pixel) <=  sqrt((img_sz/2-y)^2 ...
                                                         + (img_sz/2 - x)^2)))
                    % If with intersection, calculate intersection
                    M = area_intersect_circle_analytical([img_sz/2, ...
                                        img_sz/2, radius(z); x, y, ...
                                        ph.r_pixel]);
                    w_area = M(1,2)/max_intersection; % weight of area
                    masked_stimulus(:,:,:,:,z,x,y) = masked_stimulus(:, ...
                                                                     :, :, :, z, x,y) * w_area;
                end    
            end
        end
    end
  case 'normal'
    for x = 1: img_sz
        for y = 1: img_sz
            for z = 1:size(radius,2)
                if ( sqrt((x- round(img_sz/2))^2 + (y - round(img_sz/2))^2) ...
                     > radius(z))
                    masked_stimulus(:,:,:,:,z, x,y) = 0;
                end    
            end
        end
    end
end

% Shift the masked stimulus to the center of the basis
trans_stim = permute(masked_stimulus,[6,7,1,2,3,4,5]);
%Translation transform
xform = [ 1  0  0
             0  1  0
 (center(2)-img_sz/2) (center(1)-img_sz/2)  1 ];
tform_translate = maketform('affine',xform);

trans_stim = imtransform(trans_stim, tform_translate,...
                        'XData', [1 img_sz],...
                        'YData', [1 img_sz]);

% contrast, orientation, freq,phase, radius,x,y.
masked_stimulus = permute(trans_stim,[3,4,5,6,7,1,2]);
