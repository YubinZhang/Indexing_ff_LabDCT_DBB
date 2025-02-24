tic
input_main
check_input
load parameter_2024_02_21.mat % load best fitted parameters
parameters.beam.E_max = 110;
%% load dictionary orientations
% orientation dictionary resolution
Ori_res = 2.5;

% Set the resolution in degrees
resolution = Ori_res*degree; % Adjust for desired density

CS = crystalSymmetry('cubic');
SS = specimenSymmetry('orthorhombic');

% Generate a uniform orientation grid in the fundamental zone
cry_orien = equispacedSO3Grid(CS, 'resolution', resolution);

%%
grains_index = [];

%% load experimental data
img_y = detysize/binning;
img_z = detzsize/binning;
exp_imgs = zeros(img_y,img_z,image_num);

for i = 1:image_num
     im_i = imread(['FF_LabDCT_pure_iron\Iron_sample_tomo-B_Export' num2str(i,'%0.4d') '.tiff']);    
%       im_i = imread(['Al_ff2/proj_' num2str(i,'%0.4d') '.tif']);  %far-field single crystal Al  
%       im_i = flipud(im_i);
    exp_imgs(:,:,i) = im_i;
end

%% background process
% correct source intensity variation
for i = 1:image_num
    expi = sum(sum(exp_imgs(1:10,1:10,i)))/100;
    meani(i) = expi;
end
%figure,plot(meani/mean(meani))

% intensity rescaling
exp_imgs_scaled = zeros(img_y,img_z,image_num);
for i = 1:image_num
    exp_imgs_scaled(:,:,i) = exp_imgs(:,:,i)/meani(i)*mean(meani);
end

% remove background noise - median along the third dimension
bg_scale = median(exp_imgs_scaled,3);
exp_imgs_bg_scale_corr = exp_imgs_scaled - repmat(bg_scale,1,1,image_num);

% median filter to remove high frequency noise
exp_imgs_bg_scale_corr_medfiltered = zeros(img_y,img_z,image_num);
for i = 1:image_num
    exp_imgs_bg_scale_corr_medfiltered(:,:,i) = medfilt2(exp_imgs_bg_scale_corr(:,:,i));
end 

%% indexing with different thresholding
% defined maximum angular deviation caused by grain center of mass 
delta = 0.5/100/degree;

for threshold = [200 100 60 40]
    exp_imgs_bg_corr_bin = exp_imgs_bg_scale_corr_medfiltered>threshold;
    
    % remove small spots less than 3 pixels
    for i = 1:image_num
        exp_imgs_bg_corr_bin(:,:,i) = bwareaopen(exp_imgs_bg_corr_bin(:,:,i),3);
    end
    
    % determine experimental g-vectors 
    [exp_spot_gv_list, exp_spot_details] = find_exp_spot_gvs(exp_imgs_bg_corr_bin,exp_imgs_bg_scale_corr,parameters,Pos);

    % indexing using DBB
    for cri_completeness = [0.8,0.7,0.6,0.5,0.4]
        upper_bound_angle = Ori_res*sqrt(3)/2+delta; %delta is deviation angle caused by grain center of mass deviation
        refined_angle = 1;

        %%% update input and output
        grains_index = update_grains_spot_list(grains_index,exp_spot_gv_list,parameters,B,Ahkl,checking_angle,'nearst');
        [exp_spot_details,exp_spot_gv_list_usigned,exp_spot_gv_list_signed_dual,exp_spot_gv_list_signed] = sign_exp_spot_gv_list(exp_spot_gv_list,exp_spot_details,grains_index,exp_imgs_bg_scale_corr,parameters);

        %%% first matching - find candidate orientations
        grains = first_match(exp_spot_gv_list_usigned,cry_orien,parameters,Pos,B,Ahkl_indexing,upper_bound_angle,cri_completeness);

        %%% refine grain orientations and remove duplicates (should these be merged?)
        grains_refined = refine_candidate_grains(exp_spot_gv_list_usigned,grains,parameters,B,Ahkl_indexing,Ahkl_output,upper_bound_angle,refined_angle);

        %%% find large grains
        grains_refined_updated = update_grains_spot_list(grains_refined,exp_spot_gv_list,parameters,B,Ahkl,checking_angle,'nearst');
        [~,indx] = sort([grains_refined_updated(:).num_matched_gv]);
        grains_refined_updated = grains_refined_updated(indx);
        grains_refined_updated = grains_refined_updated(end:-1:1);
        if threshold >= 100
            grains_refined_updated = grains_refined_updated([grains_refined_updated(:).num_matched_gv]>120); %can also use completeness value
        else
            grains_refined_updated = grains_refined_updated([grains_refined_updated(:).num_matched_gv]>150); %can also use completeness value
        end

        %%% append the indexed grains to grains_index list
        if isempty(grains_index)
            grains_index = grains_refined_updated;
        else
            s = size(grains_index,2);
            for i = 1:size(grains_refined_updated,2)
                grains_index(i+s) = grains_refined_updated(i);
            end
        end
    end
end
% remove duplicate 
grains_index = unique_grains(grains_index,1);
grains_index_3hkls = update_grains_spot_list(grains_index,exp_spot_gv_list,parameters,B,Ahkl_indexing,checking_angle,'nearst');

% save results
save Indexing_result_ori2p5_2025_01_23.mat '*'
toc
%% reindexing with lower completeness value
if reindex
    grains_reindex = [];

    grains_index_updated = update_grains_spot_list(grains_index,exp_spot_gv_list,parameters,B,Ahkl,checking_angle,'nearst');
    [exp_spot_details,exp_spot_gv_list_usigned,exp_spot_gv_list_signed_dual,exp_spot_gv_list_signed] = sign_exp_spot_gv_list(exp_spot_gv_list,exp_spot_details,grains_index_updated,exp_imgs_bg_scale_corr,parameters);
    % indexing parameters
    upper_bound_angle = Ori_res*sqrt(3)/2+delta;
    cri_completeness = 0.35;
    refined_angle = 1;
    
    % first matching - find candidate orientations
    grains = first_match(exp_spot_gv_list_usigned,cry_orien,parameters,Pos,B,Ahkl_indexing,upper_bound_angle,cri_completeness);
    
    %%% refine grain orientations and remove duplicates (should these be merged?)
    grains_refined = refine_candidate_grains(exp_spot_gv_list_usigned,grains,parameters,B,Ahkl_indexing,Ahkl_output,upper_bound_angle,refined_angle);
    
    %%% find large grains
    grains_refined_updated = update_grains_spot_list(grains_refined,exp_spot_gv_list,parameters,B,Ahkl,checking_angle,'nearst');
    [~,indx] = sort([grains_refined_updated(:).num_matched_gv]);
    grains_refined_updated = grains_refined_updated(indx);
    grains_refined_updated = grains_refined_updated(end:-1:1);
    grains_refined_updated = grains_refined_updated([grains_refined_updated(:).num_matched_gv]>150); %check if these grains are all good grains.
    %%% append the indexed grains to grains_index list
    if isempty(grains_reindex)
        grains_reindex = grains_refined_updated;
    else
        s = size(grains_reindex,2);
        for i = 1:size(grains_refined_updated,2)
            grains_reindex(i+s) = grains_refined_updated(i);
        end
    end
end

grains_reindex = update_grains_spot_list(grains_reindex,exp_spot_gv_list,parameters,B,Ahkl,checking_angle,'nearst');

%% fit setup geometry - typically involves multiple iterations, 
%%% starting with the large grains indexed with high confidence and relaxed criteria
if fitdetector
    %% fit detector center
    det_center = fit_detector_center(grains_index,parameters,B)
%     parameters.detector.dety0 = det_center(1);
%     parameters.detector.detz0 = det_center(2);
    
    %% fit detector distance
    det_dist = fit_detector_dist(grains_index,parameters,B)
%       parameters.setup.Lsd = det_dist;
    
    %% fit source distance (this should not be fitted in this way, maybe a better way?)
    sou_pos = fit_source_dist(grains_index,parameters,B)
%      parameters.setup.Lss = sou_pos;
    
    %% fit detector tilt_new
    det_tilt1 = fit_det_tilt_new(grains_index,parameters)

    tilt_x = det_tilt1(1);
    tilt_y = det_tilt1(2);
    tilt_z = det_tilt1(3);
    Rx = [1 0 0; 0 cos(tilt_x) -sin(tilt_x); 0 sin(tilt_x) cos(tilt_x)];
    Ry = [cos(tilt_y) 0 sin(tilt_y); 0 1 0; -sin(tilt_y) 0 cos(tilt_y)];
    Rz = [cos(tilt_z) -sin(tilt_z) 0; sin(tilt_z) cos(tilt_z) 0; 0 0 1];

    parameters.detector.tilt_x = tilt_x;
    parameters.detector.tilt_y = tilt_y;
    parameters.detector.tilt_z = tilt_z;
    parameters.detector.tilt = Rz*Ry*Rx;
    
end

%% Check results
if check_result
    for i = 1

    grains_updated = update_grains_spot_list(grains_index(i),exp_spot_gv_list,parameters,B,Ahkl,checking_angle,'nearst');
    grains_updated.spot_list(:,19) =  grains_updated.spot_list(:,6) - grains_updated.spot_list(:,14);
    grains_updated.spot_list(:,20) =  grains_updated.spot_list(:,7) - grains_updated.spot_list(:,15);
    figure(1),imshow(sum(show_diff_imgs(grains_updated.spot_list,exp_spot_details,exp_imgs_bg_scale_corr,parameters),3), [0 500])
    hold on
    plot(grains_updated.spot_list(:,6),grains_updated.spot_list(:,7),'ro')
    plot(grains_updated.spot_list(:,14),grains_updated.spot_list(:,15),'b+')
    quiver(grains_updated.spot_list(:,14),grains_updated.spot_list(:,15),grains_updated.spot_list(:,19),grains_updated.spot_list(:,20))
    %ES_tensor(:,:,i) = fit_elastic_strain(grains_updated,B,parameters);
    hold off
    end
end

%% intensity vs energy vs polorization/Lorenzen factor...
spot_ind = grains_updated.spot_list(:,1);
hkl = grains_updated.spot_list(:,8:10);
hkln = sum(hkl.^2,2);
hkl = sort(abs(hkl),2);
ind1 = find(hkl(:,2)==3);
hkln(ind1) =[];
grain1_meanintensity = [exp_spot_details(spot_ind).MeanIntensity];
grain1_meanintensity(ind1)= [];
energy1 = grains_updated.spot_list(:,16);
energy1(ind1) = [];
tth = grains_updated.spot_list(:,17);
tth(ind1) = [];
PL = (1+cosd(tth).^2)./sind(tth)/2;
hklnu = unique(sort(hkln));
%f_asf = sind(tth/2)./(1.2398./energy1);
colorset = {'.','b.','r.','g.','k.','m.','c.','.'} ;
lengeds = {'{110}','{002}','{112}','{013}','{222}','{123}','{114}','{024}'};
figure,     hold on
t = tiledlayout(1,1);
ax1 = axes(t);
xlabel('Lorentz-Polarization factor','FontSize', 25)
ylabel('Spot Mean Intensity','FontSize', 25)

ax1.XAxisLocation = 'bottom';
ax1.YAxisLocation = 'left';
ax1.XColor = 'k';
ax1.YColor = 'k';
hold on
ax2 = axes(t);
ax2.XAxisLocation = 'top';
ax1.Box = 'on';
ax2.Box = 'on';
hold on
for i = 1:length(hklnu)-1
    indx = (hkln == hklnu(i));
    plot(ax1,PL(indx),grain1_meanintensity(indx),colorset{i})
    hold on
    plot(ax2,energy1(indx),grain1_meanintensity(indx),colorset{i})
    hold on
end
legend(lengeds)
xlabel('Energy (keV)','FontSize', 25)
hold off

%% refit the grains using all spots associated with
% grains_index_updated1 = update_grains_spot_list(grains_index,exp_spot_gv_list,parameters,B,Ahkl,checking_angle,'nearst');
grains_index_updated = grains_index;
for i = 1:length(grains_index)    
    [grains_index_updated(i).refined_ori_matrix,grains_index_updated(i).refined_pos] = refinegrain(grains_index(i).spot_list,parameters,grains_index(i).refined_ori_matrix,B,grains_index(i).refined_pos);
end
grains_index_updated = update_grains_spot_list(grains_index_updated,exp_spot_gv_list,parameters,B,Ahkl,checking_angle,'nearst');

%% find spot overlap rate for a given grain
for l = 1:length(ind)
    j = 1;
    for i = 1:length(exp_spot_details)
        if ismember(ind(l),exp_spot_details(i).grains)
            if length(exp_spot_details(i).grains)>1
                graini{j}=exp_spot_details(i).grains;
                j = j+1;
            end
        end
    end
    j-1
end