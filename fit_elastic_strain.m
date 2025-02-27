function e = fit_elastic_strain(grains,B,parameters)
UBI_init = grains.refined_ori_matrix * B; % Initial UBI matrix
%%%
options = optimoptions('fmincon','Display','off');%'iter','Algorithm','sqp-legacy');%'sqp''off');%
LB = -ones(3)*B;
UB = ones(3)*B;

[UBI,fval,~,~]  = fmincon(@(x) min_vec_angle(x,grains,parameters),UBI_init,[],[],[],[],LB,UB,@(y) nodilat(y,B),options);

% [UBI,fval, exitflag, output] = fmincon(@(x) min_vec_angle(x,grains,parameters),UB,[],[],[],[],[],[],@(y)nodilat(y,B));
% UBI = ga(@(x) sum(abs(acos(dot(Qvs_index,normc(reshape(x,3,3)*hkl))).^2)),9,[],[],[],[],[],[],@(y)nodilat(y,B));
% UBI = reshape(UBI,3,3);
% UBI = fmincon(@(x) sum(abs(acos(dot(Qvs_index,normc(x*hkl)).^2))),UB,[],[],[],[],[],[],@(y) constraint(y,B));

[~,R] = qr(UBI');
D=diag(sign(diag(R)));
B_exp=D*R;
% U_exp=Q*D;

% A_exp=eye(3)/B_exp;
% T=A_exp/Amatrix;
% eps_measured=1/2*(T+T')-eye(3);

% e = (B*inv(B_exp)+(B*inv(B_exp))')/2-eye(3);
e = (B_exp\B+(B_exp\B)')/2-eye(3);

function angle_sum = min_vec_angle(x,grains,parameters)

refined_Gvs = zeros(size(grains.spot_list,1),3);
Gv = zeros(size(grains.spot_list,1),3);
Gv_angle = zeros(size(grains.spot_list,1),1);
pos = grains.refined_pos;

for i = 1:size(grains.spot_list,1)
    rot = grains.spot_list(i,2);
    Omega = euler2u(rot*pi/180,0,0);
    spots.WeightedCentroid = grains.spot_list(i,6:7);%experimental spot
    pos_rot = Omega*pos';
    refined_Gvs(i,:) = spotpos2gvector(spots,parameters,pos_rot);
    
    Gv(i,:) = Omega*x*grains.spot_list(i,8:10)';
    %Gv_angle(i) = acos(dot(normr(Gv(i,:)),normr(refined_Gvs(i,:))));
    
    % Ensure numerical stability for acos by clamping values to [-1,1]
    dot_product = dot(normr(Gv(i, :)), normr(refined_Gvs(i, :)));
    dot_product = max(min(dot_product, 1), -1); % Avoid NaN due to floating point errors

    % Compute angular difference
    Gv_angle(i) = acos(dot_product);
end
angle_sum = sum(abs(Gv_angle).^2);


function [c ceq]=nodilat(x,B)
%NODILAT sets the condition for zero dilatational strain
%   Function to enable the zero dilatational strain criterion when
%   minimizing the angle between simulated and measured G-vectors
c=[];
x = reshape(x,3,3);
[Q,R]=qr(x);
D=diag(sign(diag(R)));
B_exp=D*R;
% A_exp=eye(3)/B_exp;
% T=A_exp/(inv(B));
% eps_measured=1/2*(T+T')-eye(3);
eps_measured = (B_exp\B+(B_exp\B)')/2-eye(3);


ceq=sum(diag(eps_measured));

