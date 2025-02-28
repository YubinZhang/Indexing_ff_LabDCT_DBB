function det_center = fit_detector_center(grains,parameters,B)

x0 = [parameters.detector.dety0,parameters.detector.detz0];

options = optimoptions('fmincon','Display','iter','Algorithm','sqp-legacy');%'sqp'

LB = x0 - [10 10];%[490 490];
UB = x0 + [10 10];%[520 520];
%this works fine
[det_center,fval,~,~]  = fmincon(@(x) fit_det_center(x,grains,parameters,B),x0,[],[],[],[],LB,UB,[],options);


function angle_sum = fit_det_center(x,grains,parameters,B)
parameters.detector.dety0 = x(1);
parameters.detector.detz0 = x(2);


Gv_angle = zeros(size(grains,2),1300);
for i = 1:size(grains,2)
    if grains(i).good_grain
        spot_list = grains(i).spot_list;
        pos = grains(i).refined_pos;
        U = grains(i).refined_ori_matrix;
        refined_Gvs = [];
        Gv = [];
        for j = 1:size(spot_list,1)            
            rot = spot_list(j,2);
            Omega = euler2u(rot*pi/180,0,0);
            spots.WeightedCentroid = spot_list(j,6:7);%experimental spot
            pos_rot = Omega*pos';
            refined_Gvs(j,:) = spotpos2gvector(spots,parameters,pos_rot);
            
            Gv(j,:) = Omega*U*B*spot_list(j,8:10)';
            % Compute angle difference (ensure numerical stability)
            dot_product = dot(normr(Gv(j, :)), normr(refined_Gvs(j, :)));
            dot_product = max(min(dot_product, 1), -1); % Avoid NaN due to floating-point errors
            Gv_angle(i,j) = acos(dot_product);

            %Gv_angle(i,j) = acos(dot(normr(Gv(j,:)),normr(refined_Gvs(j,:))));
        end
    end
end
angle_sum = sum(abs(Gv_angle(:).^2));

