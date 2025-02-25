function grains = update_grains_spot_list(grains_old,exp_spot_gv_list,parameters,B,Ahkl,cri_angle,option)

grains = grains_old;

parfor i = 1:size(grains_old,2)
    U = grains_old(i).refined_ori_matrix;
    pos = grains_old(i).refined_pos;
    exp_spot_gv_listi = update_exp_spot_gvs(exp_spot_gv_list,parameters,pos);% new function update_exp_spot_gvs()

    [match_gvs,num_match_gv,~,~] = find_matched_gvs(exp_spot_gv_listi,parameters,pos,U,B,Ahkl,cri_angle,option);
%     disp( ['grain#:' num2str(i) ', size of matched gvs (' num2str(cri_angle) ') = ' num2str(size(match_gvs,1))])

    grains(i).spot_list = match_gvs;
    grains(i).num_matched_gv = num_match_gv;

    Spot_list = forward_simulation_ff_LabDCT_v3(parameters,pos,U,B,Ahkl);
    %grains(i).expected = size(Spot_list,1);
    grains(i).completeness = grains(i).num_matched_gv/size(Spot_list,1);%grains(i).expected;
    disp( ['grain#:' num2str(i) ', size of matched gvs' '(' num2str(cri_angle) ') = ' num2str(size(match_gvs,1)) '; completeness = ' num2str(grains(i).completeness)])

end