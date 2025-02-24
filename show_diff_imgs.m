function diff_imgs = show_diff_imgs(spot_list,exp_spot_details,exp_imgs,parameters)
%peaksize = parameters.detector.peaksize;
detysize = parameters.detector.detysize;
detzsize = parameters.detector.detzsize;

% diff_imgs_mask = zeros(size(exp_imgs));
diff_imgs = zeros(size(exp_imgs));

for i = 1:size(spot_list,1)
    spot_id = spot_list(i,1);
    spot_pixel_list = [];
    spot_pixel_list = exp_spot_details(spot_id).PixelList;
    img_num = (spot_list(i,2)-parameters.setup.rotation.start)/parameters.setup.rotation.step +1;
    
    for j = 1:size(spot_pixel_list,1)
        diff_imgs(spot_pixel_list(j,2),spot_pixel_list(j,1),img_num) = exp_imgs(spot_pixel_list(j,2),spot_pixel_list(j,1),img_num);
    end
end
% imshow3D(diff_imgs)
