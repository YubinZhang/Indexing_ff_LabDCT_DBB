function [Spot_list,diff_imgs] = forward_simulation_ff_LabDCT_v3(P,Pos,U,B,Ahkl)
% function Spot_list = forward_simulation_ff_LabDCT_v3(P,Pos,U,D,S,B,Ahkl)

rot_start = P.setup.rotation.start;
rot_step = P.setup.rotation.step;
rot_end = P.setup.rotation.end;
Lsam2sou = P.setup.Lss;
Lsam2det = P.setup.Lsd;
L = P.setup.Lsd - P.setup.Lss(1);

detysize = P.detector.detysize;
detzsize = P.detector.detzsize;
dety0 = P.detector.dety0;
detz0 = P.detector.detz0;
pixelysize = P.detector.pixelysize;
pixelzsize = P.detector.pixelzsize;
beamstop = P.detector.beamstop;
peaksize = P.detector.peaksize;
peakshape = P.detector.peakshape;
Rtilt = P.detector.tilt;

E_max = P.beam.E_max;
E_min = P.beam.E_min;
K1 = P.beam.K1;
K2 = P.beam.K2;
I0 = P.beam.I0;

criangle = (atan(sqrt((detysize/2)^2+(detzsize/2)^2)*(pixelysize+pixelzsize)/2/Lsam2det)+pi/2)/degree;

image_num = (rot_end -rot_start)/rot_step+1;

Spot_list = zeros(500,15);
spot_num = 0;

% Epsilon = [-0.000223304751099440,1.95432435808773e-05,-4.69928147855185e-05;1.95432435808773e-05,-0.000251783426028829,8.28988683262049e-06;-4.69928147855185e-05,8.28988683262049e-06,0.000475223659616386];
Epsilon = zeros(3,3);
for rot = rot_start:rot_step:rot_end
    Omega = euler2u(rot*pi/180,0,0);
    voxelpos = Omega*Pos';
    alpha = atan(sqrt((voxelpos(2)-Lsam2sou(2))^2+(voxelpos(3)-Lsam2sou(3))^2)/(voxelpos(1)-Lsam2sou(1)));

    %diffraction center
    N_det = Rtilt*[-1 0 0]';
    center_vector = normc(cross(normc(cross(voxelpos,voxelpos-Lsam2sou')), N_det));
    alpha1 = acos(dot(center_vector,normc(voxelpos-Lsam2sou')));
    center_length = L*sin(alpha)/sin(alpha1);
    center = center_vector*center_length;
    %center = [L, voxelpos(2)*L/(-Lsam2sou(1)+voxelpos(1)), voxelpos(3)*L/(-Lsam2sou(1)+voxelpos(1))];
%     alpha = atan(sqrt(voxelpos(2)^2+voxelpos(3)^2)/(-Lsam2sou(1)+voxelpos(1)));
    K_in = voxelpos - Lsam2sou';
%     D = eye(3)+Epsilon;%this is not the correct way to include strain% strain comes into A matrix, and then propagates to B matrix
    Gw = Omega*U*B*Ahkl(:,1:3)';%Omega*D^-1*U*B*Ahkl(:,1:3)'; 
    Ahkl1 = Ahkl(:,1:5)';
%    Ahkl1 = Ahkl1(:,Gw(1,:)<0);%Gw(x) negtive - This is too simplified.some negtive g(3) could lead to diffraction if the grain sits at the bottom of the gauge volume
%     Gw = Gw(:,Gw(1,:)<0);
    
    beta = acos(dot(repmat(normc(K_in),1,size(Gw,2)),normc(Gw),1));
    
    Gw = Gw(:,beta>=90*pi/180);%more general 
    Ahkl1 = Ahkl1(:,beta>=90*pi/180);
    
    beta = beta(beta>=90*pi/180);
    Gw = Gw(:,beta<criangle*pi/180);
    Ahkl1 = Ahkl1(:,beta<criangle*pi/180);
    theta = beta(beta<criangle*pi/180)-pi/2;% diffraction angle
    sintth = sin(2*theta);
    d = 1./vecnorm(Gw)*2*pi;
    lambdahkl = 2*d.*sin(theta);
    Energy = 12.398./lambdahkl;
    Energy(Energy<E_min) = 1000;%special cases
    
    Gw = Gw(:,Energy<E_max);
    Ahkl1 = Ahkl1(:,Energy<E_max);
    sintth = sintth(Energy<E_max);
    theta = theta(Energy<E_max);
    lambdahkl = lambdahkl(Energy<E_max);
    Energy = Energy(Energy<E_max);%special cases
    
    if size(theta,2)>0
        diff_dir = normc(cross(normc(cross(Gw,repmat(normc(K_in),1,size(Gw,2)))),repmat(N_det,1,size(Gw,2))));
        %K_out = Gw + repmat((K_in)',1,size(theta,2));
        phix= acos(dot(repmat(normc(K_in),1,size(theta,2)),diff_dir,1));
        phiy = phix-2*theta;
        difflength = norm([Lsam2det 0 0]' + center -voxelpos)./sin(phiy).*sintth;

        diff_vector = difflength.*diff_dir + repmat(center,1,size(theta,2));

        diff_vector(1,:) = 0;
        diff_vector_rot = Rtilt'*diff_vector;


        dety = dety0+diff_vector_rot(2,:)/pixelysize;
        detz = detz0+diff_vector_rot(3,:)/pixelzsize;
        
        spot_pos = [diff_vector_rot(1,:); dety; detz]';
        
        spot_num = spot_num + size(theta,2);
        
        Spot_list(spot_num-size(theta,2)+1:spot_num,1) = spot_num-size(theta,2)+1:spot_num;
        Spot_list(spot_num-size(theta,2)+1:spot_num,2) = rot;
        Spot_list(spot_num-size(theta,2)+1:spot_num,3:5) = Ahkl1(1:3,:)';
        Spot_list(spot_num-size(theta,2)+1:spot_num,6) = Ahkl1(5,:);
        Spot_list(spot_num-size(theta,2)+1:spot_num,7:9) = Gw';
        Spot_list(spot_num-size(theta,2)+1:spot_num,10) = Energy;%w*180/pi;
        Spot_list(spot_num-size(theta,2)+1:spot_num,11)= 2*theta*180/pi;
        Spot_list(spot_num-size(theta,2)+1:spot_num,12) = spot_pos(:,2);%dety;
        Spot_list(spot_num-size(theta,2)+1:spot_num,13) = spot_pos(:,3);%detz;
        Spot_list(spot_num-size(theta,2)+1:spot_num,14) = difflength;
        Iwh = (lambdahkl/(12.398/E_max)-1)./(lambdahkl).^2;
        %Diffracted intensity
        Spot_list(spot_num-size(theta,2)+1:spot_num,15) = K1*K2*I0*Ahkl1(5,:).*Iwh;%to be updated with real intensity
    end
end
Spot_list(Spot_list(:,12)<1,:)=[];
Spot_list(Spot_list(:,12)>=detysize-1,:)=[];
Spot_list(Spot_list(:,13)<1,:)=[];
Spot_list(Spot_list(:,13)>=detzsize-1,:)=[];

indx_beamstop = find(Spot_list(:,12) > beamstop(3) & Spot_list(:,12) < beamstop(4) ...
    & Spot_list(:,13) > beamstop(1) & Spot_list(:,13) < beamstop(2));

Spot_list(indx_beamstop,:)=[];

if nargout ==2
    diff_imgs = zeros(detysize,detzsize,image_num,'uint16');
    
%     TO update with different spot size and shape (defined by
%     parameter.indexing.spot
    if peakshape == 0
        peakwsig = 0;
        pixellimit = 0;
    end
    
    num_spot = size(Spot_list,1);
    for spot_i = 1:num_spot
        int=Spot_list(spot_i,15);
        img_i = (Spot_list(spot_i,2)-rot_start)/rot_step+1;
        dety=round(Spot_list(spot_i,12))+1;
        detz=round(Spot_list(spot_i,13))+1;
        if (-peaksize*pixellimit >= dety) || (dety >= detysize+peaksize*pixellimit) || (-peaksize*pixellimit >= detz) || (detz >= detzsize+peaksize*pixellimit)

        else
            if peakshape == 0                
                for k1=dety-peaksize:dety+peaksize
                    for k2=detz-peaksize:detz+peaksize
                        if (0 < k1) && (k1 <= detysize) && (0 < k2) && (k2 <= detzsize)
                            diff_imgs(k2,k1,img_i)=diff_imgs(k2,k1,img_i)+int;
                        end
                    end
                end
            end
        end
    end
    diff_imgs(beamstop(1):beamstop(2),beamstop(3):beamstop(4),:)=0;        %remove beamstop area
    % for i = 1:image_num
    %     frame = diff_imgs(:,:,i);        %remove beamstop area
    %     frame = fliplr(flipud(frame));
    %     diff_imgs(:,:,i) = frame;
    % end
end