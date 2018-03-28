function [fromXYZ1] = bestfit(toXYZ, fromXYZ, mass, sum_of_mass)
%massx is a N by 1 atomic mass vector
%function [R,T,eRMSD] = kabsch(toXYZ:reference:t2,fromXYZ:current:t1)
% find R and T that best maps
% also find:
% eRMSD = std(R*fromXYZ + T - toXYZ)
% which will rotate (R: 3x3 matrix) and translate (T: 3 x 1 vector)
% coordinates: fromXYZ (3 x N matrix) to toXYZ (3 x N matrix) and returns root-mean-squared error (eRMSD)
% ||R*fromXYZ + T - toXYZ||^2

    ln1 = size(fromXYZ,1);
    mass_center1 = sum(fromXYZ .* mass) / sum_of_mass;

    t1 = fromXYZ-repmat(mass_center1,ln1,1);

    [u,~,w] = svd(toXYZ' * t1);
    R = w*[1 0 0 ; 0 1 0; 0 0 1]*u';
    
    fromXYZ1 = t1*R;