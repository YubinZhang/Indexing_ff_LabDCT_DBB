function [match_gvs,num_match_gv,num_multi_match,Spot_list] = find_matched_gvs(exp_spot_gv_list,parameters,Pos,U,B,Ahkl,criangle,option)

num_match = 0;
num_multi_match = 0;
num_match_gv = 0;
match_gvs = [];
rot_start = parameters.setup.rotation.start;
rot_step = parameters.setup.rotation.step;
rot_end = parameters.setup.rotation.end;

Spot_list = forward_simulation_ff_LabDCT_v3(parameters,Pos,U,B,Ahkl);
for rot_i = rot_start:rot_step:rot_end
    theo_spot_rot_i = Spot_list(Spot_list(:,2)==rot_i,:);
    if ~isempty(theo_spot_rot_i)
        exp_spot_rot_i = exp_spot_gv_list(exp_spot_gv_list(:,2)==rot_i,:);
        
        for spot_i = 1:size(theo_spot_rot_i,1)
            angle = acosd(dot(exp_spot_rot_i(:,3:5)',repmat(normr(theo_spot_rot_i(spot_i,7:9))',1,size(exp_spot_rot_i,1))));
            ind = find(angle<criangle);% add a option and a criterion to select only the one with min angle deviation

            if strcmp(option, 'closest')
                distance = sqrt(sum((exp_spot_rot_i(:,6:7)-repmat(theo_spot_rot_i(spot_i,12:13),size(exp_spot_rot_i,1),1)).^2,2));
                [~,b] = min(distance);

                if min(distance)<5 && angle(b)<criangle 
                    
                    match_gvs(num_match+1,1:7) = exp_spot_rot_i(b,1:7);
                    match_gvs(num_match+1,8:10) = theo_spot_rot_i(spot_i,3:5);
                    match_gvs(num_match+1,11:13) = theo_spot_rot_i(spot_i,7:9);
                    match_gvs(num_match+1,14:15) = theo_spot_rot_i(spot_i,12:13);
                    match_gvs(num_match+1,16:17) = theo_spot_rot_i(spot_i,10:11);
                    match_gvs(num_match+1,18) = angle(b);
                    num_match = num_match + 1;
                    num_match_gv = num_match_gv +1;
                end
            else
                if ~isempty(ind)
                    if length(ind)>1
                        if strcmp(option, 'all')
                            match_gvs(num_match+1:num_match + length(ind),1:7) = exp_spot_rot_i(ind,1:7);
                            match_gvs(num_match+1:num_match + length(ind),8:10) = repmat(theo_spot_rot_i(spot_i,3:5),length(ind),1);
                            match_gvs(num_match+1:num_match + length(ind),11:13) = repmat(theo_spot_rot_i(spot_i,7:9),length(ind),1);
                            match_gvs(num_match+1:num_match + length(ind),14:15) = repmat(theo_spot_rot_i(spot_i,12:13),length(ind),1);
                            match_gvs(num_match+1:num_match + length(ind),16:17) = repmat(theo_spot_rot_i(spot_i,10:11),length(ind),1);%theo_spot_rot_i(spot_i,10:11);
                            match_gvs(num_match+1:num_match + length(ind),18) = angle(ind);

%                             match_gvs(num_match+1:num_match + length(ind),16) = angle(ind);
                            num_match = num_match + length(ind);
                            num_match_gv = num_match_gv +1;
                            %if length(ind) > 1
                            num_multi_match = num_multi_match +1;
                        elseif strcmp(option,'nearst')
                            [~,b] = min(angle);
                            %if sqrt(sum((exp_spot_rot_i(b,6:7)-theo_spot_rot_i(spot_i,12:13)).^2))<5
                            match_gvs(num_match+1,1:7) = exp_spot_rot_i(b,1:7);
                            match_gvs(num_match+1,8:10) = theo_spot_rot_i(spot_i,3:5);
                            match_gvs(num_match+1,11:13) = theo_spot_rot_i(spot_i,7:9);
                            match_gvs(num_match+1,14:15) = theo_spot_rot_i(spot_i,12:13);
                            match_gvs(num_match+1,16:17) = theo_spot_rot_i(spot_i,10:11);
                            match_gvs(num_match+1,18) = angle(b);

%                             match_gvs(num_match+1,16) = angle(b);
                            num_match = num_match + 1;
                            num_match_gv = num_match_gv +1;
                            %if length(ind) > 1
                            %num_multi_match = num_multi_match +1;
                            %end
                        end
                    elseif length(ind)==1
                        %                     if strcmp(option, 'all')
                        match_gvs(num_match+1,1:7) = exp_spot_rot_i(ind,1:7);
                        match_gvs(num_match+1,8:10) = theo_spot_rot_i(spot_i,3:5);
                        match_gvs(num_match+1,11:13) = theo_spot_rot_i(spot_i,7:9);
                        match_gvs(num_match+1,14:15) = theo_spot_rot_i(spot_i,12:13);
                        match_gvs(num_match+1,16:17) = theo_spot_rot_i(spot_i,10:11);
                        match_gvs(num_match+1,18) = angle(ind);
% match_gvs(num_match+1,16) = angle(ind);
                        num_match = num_match + 1;
                        num_match_gv = num_match_gv +1;
                        %                     elseif strcmp(option,'nearst')
                        %                         %if sqrt(sum((exp_spot_rot_i(ind,6:7)-theo_spot_rot_i(spot_i,12:13)).^2))<5
                        %                             match_gvs(num_match+1,1:7) = exp_spot_rot_i(ind,1:7);
                        %                             match_gvs(num_match+1,8:10) = theo_spot_rot_i(spot_i,3:5);
                        %                             match_gvs(num_match+1,11:13) = theo_spot_rot_i(spot_i,7:9);
                        %                             match_gvs(num_match+1,14:15) = theo_spot_rot_i(spot_i,12:13);
                        %                             match_gvs(num_match+1,16) = angle(ind);
                        %                             num_match = num_match + 1;
                        %                             num_match_gv = num_match_gv +1;
                        %                         %end
                        %                     end
                    end
                end
            end
        end
    end
end
end

