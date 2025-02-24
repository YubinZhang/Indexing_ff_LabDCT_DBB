function [exp_spot_details,exp_spot_gv_list_usigned,exp_spot_gv_list_signed_dual,exp_spot_gv_list_signed] = sign_exp_spot_gv_list(exp_spot_gv_list,exp_spot_details,grains_refined,exp_imgs_copy,parameters)
for i = 1:size(exp_spot_gv_list,1)
    exp_spot_details(i).grains = [];
end

for i = 1:size(grains_refined,2)
     if grains_refined(i).good_grain %not(grains_unique2(i).overlapped)         
            [~,idx] = unique(grains_refined(i).spot_list(:,1));
            spot_list = grains_refined(i).spot_list(idx,:);
            for j = 1:size(spot_list,1)
                if spot_list(j,18)<0.2
                    exp_spot_details(spot_list(j)).grains = [exp_spot_details(spot_list(j)).grains; i];
                end
            end
     end
end

ind = zeros(size(exp_spot_details));
for i = 1:size(exp_spot_details,2)
    if size(exp_spot_details(i).grains,1)==1 %uniquely indexed spots
        ind(i) = 1;
    elseif size(exp_spot_details(i).grains,1)>1 %overlapped spots
        ind(i) = 2;
    end
end

exp_spot_gv_list_usigned = exp_spot_gv_list(ind==0,:);
% imshow3D(show_diff_imgs(exp_spot_gv_list_usigined,exp_imgs_bg_corr_bin,parameters),[])
figure(1),imshow(sum(show_diff_imgs(exp_spot_gv_list_usigned,exp_spot_details,exp_imgs_copy,parameters),3))
%% more than one grain
exp_spot_gv_list_signed_dual = exp_spot_gv_list(ind==2,:);
%imshow3D(show_diff_imgs(exp_spot_gv_list_sigined_dual,exp_imgs_bg_corr_bin,parameters),[])
figure(2),imshow(sum(show_diff_imgs(exp_spot_gv_list_signed_dual,exp_spot_details,exp_imgs_copy,parameters),3))
%% all indexed
exp_spot_gv_list_signed = exp_spot_gv_list(ind>0,:);
%imshow3D(show_diff_imgs(exp_spot_gv_list_sigined,exp_imgs_bg_corr_bin,parameters),[])
figure(3),imshow(sum(show_diff_imgs(exp_spot_gv_list_signed,exp_spot_details,exp_imgs_copy,parameters),3))
