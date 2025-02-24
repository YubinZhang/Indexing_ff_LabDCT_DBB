function det_tilt = fit_det_tilt_new(grains,parameters)
%x0 = eye(3);

x0 = [0 0 0];
ConstraintFunction = @constraintfun;
options = optimoptions('fmincon','Display','iter','Algorithm','sqp-legacy');%'sqp'
% options.InitialPopulationMatrix = x0(:);
% % options = optimoptions(@ga,'UseVectorized',true);
% options = optimoptions(@fmincon,'PlotFcn',{@gaplotbestf});
LB = [-1 -1 -1 ];%-ones(3);
UB = [1 1 1 ];%ones(3);
%this works fine
[det_tilt,fval,~,~]  = fmincon(@(x) fit_spot_dist(x,grains,parameters),x0,[],[],[],[],LB,UB,[],options);


function dist = fit_spot_dist(x,grains,parameters)
tilt_x = x(1);
tilt_y = x(2);
tilt_z = x(3);
Rx = [1 0 0; 0 cos(tilt_x) -sin(tilt_x); 0 sin(tilt_x) cos(tilt_x)];
Ry = [cos(tilt_y) 0 sin(tilt_y); 0 1 0; -sin(tilt_y) 0 cos(tilt_y)];
Rz = [cos(tilt_z) -sin(tilt_z) 0; sin(tilt_z) cos(tilt_z) 0; 0 0 1];

parameters.detector.tilt_x = tilt_x;
parameters.detector.tilt_y = tilt_y;
parameters.detector.tilt_z = tilt_z;
parameters.detector.tilt = Rz*Ry*Rx;
Rtilt = parameters.detector.tilt;


Lsam2sou = parameters.setup.Lss;
Lsam2det = parameters.setup.Lsd;
L = parameters.setup.Lsd - parameters.setup.Lss(1);
detysize = parameters.detector.detysize;
detzsize = parameters.detector.detzsize;
dety0 = parameters.detector.dety0;
detz0 = parameters.detector.detz0;
pixelysize = parameters.detector.pixelysize;
pixelzsize = parameters.detector.pixelzsize;


% Gv_angle = zeros(size(grains,2),400);
dist = 0;
for i = 1:size(grains,2)
    if grains(i).good_grain
%         grains_updated = update_grains_spot_list(grains(i),exp_spot_gv_list,parameters,B,Ahkl,check_angle,'nearst');
% 
%         dist = dist + sum(sqrt((grains_updated.spot_list(:,14)-grains_updated.spot_list(:,6)).^2+(grains_updated.spot_list(:,15)-grains_updated.spot_list(:,7).^2)));
%         % 
        spot_list = grains(i).spot_list;
        %pos = x(4:6);
        pos = grains(i).refined_pos;
        U = grains(i).refined_ori_matrix;
        
        for j = 1:size(spot_list,1)            
            rot = spot_list(j,2);
            Omega = euler2u(rot*pi/180,0,0);

            %spots.WeightedCentroid = spot_list(j,6:7);%experimental spot
            
            voxelpos = Omega*pos';

            %%%%%%
            alpha = atan(sqrt((voxelpos(2)-Lsam2sou(2))^2+(voxelpos(3)-Lsam2sou(3))^2)/(voxelpos(1)-Lsam2sou(1)));
        
            
            %diffraction center
            N_det = parameters.detector.tilt*[-1 0 0]';
            center_vector = normc(cross(normc(cross(voxelpos,voxelpos-Lsam2sou')), N_det));
            alpha1 = acos(dot(center_vector,normc(voxelpos-Lsam2sou')));
            center_length = L*sin(alpha)/sin(alpha1);
            center = center_vector*center_length;
            K_in = voxelpos - Lsam2sou';
    
            Gw = spot_list(j,11:13)';
    
            beta = acos(dot(repmat(normc(K_in),1,size(Gw,2)),normc(Gw),1));
            theta = beta-pi/2;% diffraction angle
            sintth = sin(2*theta);
    
            diff_dir = normc(cross(normc(cross(Gw,repmat(normc(K_in),1,size(Gw,2)))),repmat(N_det,1,size(Gw,2))));
            phix= acos(dot(repmat(normc(K_in),1,size(theta,2)),diff_dir,1));
            phiy = phix-2*theta;
            difflength = norm([Lsam2det 0 0]' + center -voxelpos)./sin(phiy).*sintth;
            diff_vector = difflength.*diff_dir + repmat(center,1,size(theta,2));

            
            diff_vector(1,:) = 0;
            diff_vector_rot = Rtilt'*diff_vector;
            
            
            dety = dety0+diff_vector_rot(2,:)/pixelysize;
            detz = detz0+diff_vector_rot(3,:)/pixelzsize;

%             dety = dety0+diff_vector(2,:)/pixelysize;
%             detz = detz0+diff_vector(3,:)/pixelzsize;
        
            dist = dist + norm(spot_list(j,6:7)-[dety detz]);

            %%%%
            %refined_Gvs(j,:) = spotpos2gvector(spots,parameters,pos_rot);
            
            %Gv(j,:) = Omega*U*B*spot_list(j,8:10)';
            %Gv_angle(i,j) = acos(dot(normr(Gv(j,:)),normr(refined_Gvs(j,:))));
        end
    end
%     x,dist
end
% angle_sum = sum(abs(Gv_angle(:)));

function [c, ceq]= constraintfun(x)
% tilt_x = x(1);
% tilt_y = x(2);
% tilt_z = x(3);
% Rx = [1 0 0; 0 cos(tilt_x) -sin(tilt_x); 0 sin(tilt_x) cos(tilt_x)];
% Ry = [cos(tilt_y) 0 sin(tilt_y); 0 1 0; -sin(tilt_y) 0 cos(tilt_y)];
% Rz = [cos(tilt_z) -sin(tilt_z) 0; sin(tilt_z) cos(tilt_z) 0; 0 0 1];
% 
% U = Rz*Ry*Rx;
U = reshape(x(4:12),3,3);
c = [];%sum(c(:));
ceq = U*U'-eye(3);

