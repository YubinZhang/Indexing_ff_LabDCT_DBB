function [refined_ori,refined_pos,refined_error] = refinegrain(grain_spot_list,parameters,U,B,pos)

x0 = [U(:);pos']';
ConstraintFunction = @(x) constraintfun(x);
options = optimoptions('fmincon','UseParallel',true,'Display','off');%'iter','Algorithm','sqp-legacy');%'sqp''off');%
LB = [-1 -1 -1 -1 -1 -1 -1 -1 -1 -0.5 -0.5 -0.5];
UB = [1 1 1 1 1 1 1 1 1 0.5 0.5 0.5];
%this works fine
[refined_grain,fval,~,~]  = fmincon(@(x) fitgrains(x,grain_spot_list,parameters,B),x0,[],[],[],[],LB,UB,ConstraintFunction,options);

% obj = @(x)fitgrains(x,grain_spot_list,parameters,B);
% problem = createOptimProblem('fmincon','objective',obj,'x0',x0,'lb',LB,'ub',UB,'nonlcon',ConstraintFunction,'options',options);
% gs = GlobalSearch;
% [refined_grain,fval,flg,og] = run(gs,problem)
refined_ori = reshape(refined_grain(1:9),3,3);
refined_pos = refined_grain(10:12);
refined_error = fval;

function angle_sum = fitgrains(x,spot_list,parameters,B)
U = reshape(x(1:9),3,3);
pos = x(10:12);

refined_Gvs = zeros(size(spot_list,1),3);

for i = 1:size(spot_list,1)
    rot = spot_list(i,2);
    Omega = euler2u(rot*pi/180,0,0);
    spots.WeightedCentroid = spot_list(i,6:7);%experimental spot
    pos_rot = Omega*pos';
    refined_Gvs(i,:) = spotpos2gvector(spots,parameters,pos_rot);
    
    Gv(i,:) = Omega*U*B*spot_list(i,8:10)';
    % Compute angle difference (ensure numerical stability)
    dot_product = dot(normr(Gv(i, :)), normr(refined_Gvs(i, :)));
    dot_product = max(min(dot_product, 1), -1); % Avoid NaN due to floating-point errors
    Gv_angle(i) = acos(dot_product);
    %Gv_angle(i) = acos(dot(normr(Gv(i,:)),normr(refined_Gvs(i,:))));
end
angle_sum = sum(abs(Gv_angle));

function [c, ceq]= constraintfun(x)
U = reshape(x(1:9),3,3);
c = [];%sum(c(:));
ceq = U*U'-eye(3);
