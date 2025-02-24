function exp_spot_gv_list_update = update_exp_spot_gvs(exp_spot_gv_list,P,pos)

exp_spot_gv_list_update = exp_spot_gv_list;

for i = 1:size(exp_spot_gv_list,1)
    rot = exp_spot_gv_list(i,2);
    rot_pos = euler2u(rot*degree,0,0)*pos';
    spots.WeightedCentroid(1) = exp_spot_gv_list(i,6);
    spots.WeightedCentroid(2) = exp_spot_gv_list(i,7);      
    Gvs_i = spotpos2gvector(spots,P,rot_pos);
    exp_spot_gv_list_update(i,3:5) = Gvs_i;
end
