function grains_refined = refine_candidate_grains(exp_spot_gv_list,grains,parameters,B,Ahkl_indexing,Ahkl_output,upper_bound_angle,refined_angle)


parfor i = 1:size(grains,2)
    %initial refinement using first guess position
    if isempty(grains(i).num_matched_gv) || grains(i).num_matched_gv ==0
        continue
    end
    U = grains(i).ori_matrix;
    pos = grains(i).pos;
    [refined_ori_matrix_init,refined_pos_init] = refinegrain(grains(i).spot_list,parameters,U,B,pos);
    
    %update experimental g-vector using the refined_position
%     [exp_spot_gv_listi, ~] = find_exp_spot_gvs(exp_imgs_bg_corr_bin,exp_imgs_bg_corr,parameters,refined_pos_init);% new function update_exp_spot_gvs()
    exp_spot_gv_listi = update_exp_spot_gvs(exp_spot_gv_list,parameters,refined_pos_init);% new function update_exp_spot_gvs()
    
    %update matched g-vector list
    [match_gvs_init,~,~,~] = find_matched_gvs(exp_spot_gv_listi,parameters,refined_pos_init,refined_ori_matrix_init,B,Ahkl_indexing,upper_bound_angle,'nearst');
    disp( ['grain#:' num2str(i) ', size of init matched gvs (2p5) = ' num2str(size(match_gvs_init,1))])
    
    %further refinement of the orientation and position
    [refined_ori_matrix,refined_pos] = refinegrain(match_gvs_init,parameters,refined_ori_matrix_init,B,refined_pos_init);
    
    %%update experimental g-vector using newly refined position 
    exp_spot_gv_listi = update_exp_spot_gvs(exp_spot_gv_list,parameters,refined_pos);

    %update matched g-vector list using more critical angular tolerance 
    [match_gvs,num_match_gv,num_multi_match,Spot_list] = find_matched_gvs(exp_spot_gv_listi,parameters,refined_pos,refined_ori_matrix,B,Ahkl_output,refined_angle,'nearst');
    disp( ['grain#:' num2str(i) ', size of matched gvs(1) = ' num2str(size(match_gvs,1))])
    
    %further refinement of the orientation and position (typically no need)
    flag = 0;
    while flag
        %second refinement using first refined position
        
        [match_gvs,num_match_gv,num_multi_match,Spot_list] = find_matched_gvs(exp_spot_gv_listi,parameters,refined_pos,refined_ori_matrix,B,Ahkl_output,refined_angle,'nearst');
        disp( ['grain#:' num2str(i) ', size of matched gvs(1) = ' num2str(size(match_gvs,1))])
        
        [refined_ori_matrix_new,refined_pos_new] = refinegrain(match_gvs,parameters,refined_ori_matrix,B,refined_pos);
        exp_spot_gv_listi_new = update_exp_spot_gvs(exp_spot_gv_list,parameters,refined_pos_new);
        [match_gvs_new,num_match_gv_new,~,~] = find_matched_gvs(exp_spot_gv_listi_new,parameters,refined_pos_new,refined_ori_matrix_new,B,Ahkl_indexing,refined_angle,'closest');
        disp( ['grain#:' num2str(i) ', size of matched gvs(1) aftere further refinement = ' num2str(size(match_gvs,1))])
        
        if size(match_gvs_new,1)>size(match_gvs,1)
            match_gvs = match_gvs_new;
            num_match_gv = num_match_gv_new;
            refined_ori_matrix = refined_ori_matrix_new;
            refined_pos = refined_pos_new;
            %             Spot_list = Spot_list_new;
            flag = 1;
        else
            flag = 0;
        end
    end
    
    %save results
    if num_match_gv/size(Spot_list,1)>0.35 && size(unique(match_gvs(:,8:10),'rows'),1)>size(Ahkl_indexing,1)/2
        grains(i).good_grain = 1;
        grains(i).num_matched_gv = num_match_gv;
        grains(i).num_multi_matched = num_multi_match;
        grains(i).spot_list = match_gvs;
        grains(i).refined_ori_matrix = refined_ori_matrix;
        grains(i).refined_pos = refined_pos;
        mis = misori2(U',refined_ori_matrix');
        grains(i).mis = mis;
        grains(i).dict = sqrt(sum((pos-refined_pos).^2));
        
    else % seems no need to define bad grains...
        grains(i).good_grain = 0;
        grains(i).num_matched_gv = num_match_gv;
        grains(i).num_multi_matched = num_multi_match;
        grains(i).spot_list = match_gvs;
        grains(i).refined_ori_matrix = refined_ori_matrix;
        grains(i).refined_pos = refined_pos;
        mis = misori2(U',refined_ori_matrix');
        grains(i).mis = mis;
        grains(i).dict = sqrt(sum((pos-refined_pos).^2));
    end
end
%remove duplicates
grains_refined = unique_grains(grains,1);