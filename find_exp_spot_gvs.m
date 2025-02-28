function [exp_spot_gv_list, exp_spot_details] = find_exp_spot_gvs(exp_imgs_bin,exp_imgs_intensity,P,pos)
% function [exp_spot_gv_list, exp_spot_details] = find_exp_spot_gvs(exp_imgs_bin,exp_imgs_intensity,P,pos)
% find experimental diffraction vectors from diffraction spots

rot_start = P.setup.rotation.start;
rot_step = P.setup.rotation.step;

spot_num = 0;
SPOTLIST =  struct('PixelList', [],'MeanIntensity', []);

exp_spot_gv_list = zeros(100000,7);
exp_spot_details(100000) = struct(SPOTLIST);

for i = 1:size(exp_imgs_bin,3)
    rot = rot_start + (i-1)*rot_step;
    [labeledImage, numBlobs] = bwlabel(exp_imgs_bin(:,:,i)); 
    spots = regionprops(labeledImage, exp_imgs_intensity(:,:,i),'Area','WeightedCentroid','PixelList','MeanIntensity');
    if  numBlobs >0
        
        rot_pos = euler2u(rot*degree,0,0)*pos';
        Gvs_i = spotpos2gvector(spots,P,rot_pos);
       
        exp_spot_gv_list(spot_num+1:spot_num + numBlobs,1) = [spot_num+1:spot_num+numBlobs];
        exp_spot_gv_list(spot_num+1:spot_num + numBlobs,2) = rot;
        exp_spot_gv_list(spot_num+1:spot_num + numBlobs,3:5) = Gvs_i;
        exp_spot_gv_list(spot_num+1:spot_num + numBlobs,6:7) = reshape([spots.WeightedCentroid],2,numBlobs)';
        for j = 1:numBlobs
            exp_spot_details(spot_num+j).PixelList = spots(j).PixelList;
            exp_spot_details(spot_num+j).MeanIntensity = spots(j).MeanIntensity;
        end
        
        spot_num = spot_num + numBlobs;
    end
end
if spot_num<100000
    exp_spot_gv_list(spot_num +1:end,:) = [];
    exp_spot_details(spot_num +1:end) = [];
end