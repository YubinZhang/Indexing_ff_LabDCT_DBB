function sou_pos = fit_source_dist(grains,parameters,B)

x0 = [-100 0 0];
ConstraintFunction = @constraintfun;
options = optimoptions('fmincon','Display','iter','Algorithm','sqp-legacy');%'sqp'
% options.InitialPopulationMatrix = x0(:);
% % options = optimoptions(@ga,'UseVectorized',true);
% options = optimoptions(@fmincon,'PlotFcn',{@gaplotbestf});
LB = [-120 -1 -1];
UB = [-90 1 1];
%this works fine
[sou_pos,fval,~,~]  = fmincon(@(x) fit_sou_dist(x,grains,parameters,B),x0,[],[],[],[],LB,UB,[],options);


function angle_sum = fit_sou_dist(x,grains,parameters,B)
parameters.setup.Lss = x;

Gv_angle = zeros(size(grains,2),400);
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
            Gv_angle(i,j) = acos(dot(normr(Gv(j,:)),normr(refined_Gvs(j,:))));
        end
    end
end
angle_sum = sum(abs(Gv_angle(:)));

function [c, ceq]= constraintfun(x)
c = [];%sum(c(:));
ceq = x(1)+100;

