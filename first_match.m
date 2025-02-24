function grains = first_match(exp_spot_gv_list,cry_orien,parameters,Pos,B,Ahkl_indexing,upper_bound_angle,cri_completeness)
% function grains = first_match(exp_spot_gv_list,cry_orien,parameters,Pos,B,Ahkl_indexing,upper_bound_angle,cri_completeness)

GRAINSTR = struct(...
    'num_matched_gv',  0,...
    'num_multi_matched', 0,...
    'spot_list',       [],...
    'ori_matrix',      [0 0 0; 0 0 0; 0 0 0],...
    'pos',      [0,0,0],...
    'refined_ori_matrix', zeros(3),...
    'refined_pos', [0 0 0],...
    'overlapped', 0,...
    'good_grain', 0,...
    'dict',       0,...
    'strain',   [0 0 0; 0 0 0; 0 0 0]);
grains(size(cry_orien,1)) = struct(GRAINSTR);


parfor ori_i = 1:length(cry_orien)
    U = cry_orien(ori_i).matrix;
    [match_gvs,num_match_gv,num_multi_match,Spot_list] = find_matched_gvs(exp_spot_gv_list,parameters,Pos,U,B,Ahkl_indexing,upper_bound_angle,'all');
    
    if num_match_gv/size(Spot_list,1)>cri_completeness && size(unique(match_gvs(:,8:10),'rows'),1)>size(Ahkl_indexing,1)/2
        disp('found a grain');
        ori_i
        if num_multi_match > 0.45*num_match_gv
            disp('possible orientation overlap')
            grains(ori_i).overlapped = 1;
        else
            grains(ori_i).overlapped = 0;
        end
        grains(ori_i).num_matched_gv = num_match_gv;
        grains(ori_i).num_multi_matched = num_multi_match;
        grains(ori_i).spot_list = match_gvs;
        grains(ori_i).ori_matrix = U;
        grains(ori_i).pos = Pos;
    end
end

ind = zeros(size(grains));
for i = 1:size(grains,2)
    if isempty(grains(i).num_matched_gv) || grains(i).num_matched_gv ==0
        ind(i) = 1;
    end
end
grains(ind==1)=[];