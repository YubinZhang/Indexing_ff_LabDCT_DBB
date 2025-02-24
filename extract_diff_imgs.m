function [diff_imgs,spot_list] = extract_diff_imgs(spot_list,exp_imgs,parameters,U,D,S,B,pos)

peaksize = parameters.detector.peaksize;
detysize = parameters.detector.detysize;
detzsize = parameters.detector.detzsize;
img_nums = (parameters.setup.rotation.end-parameters.setup.rotation.start)/parameters.setup.rotation.step +1;
exp_imgs_flipped = zeros(size(exp_imgs));
diff_imgs_flipped = zeros(size(exp_imgs));
diff_imgs = zeros(size(exp_imgs));

dety0 = parameters.detector.dety0;
detz0 = parameters.detector.detz0;
pixelysize = parameters.detector.pixelysize;
pixelzsize = parameters.detector.pixelzsize;

Lsam2sou = parameters.setup.Lss;
Lsam2det = parameters.setup.Lsd; 
L = Lsam2sou + Lsam2det;

for i = 1:img_nums
    exp_imgs_flipped(:,:,i) = flipud(fliplr(exp_imgs(:,:,i))); 
end 

for i = 1:size(spot_list,1)
    img_num = (spot_list(i,2)-parameters.setup.rotation.start)/parameters.setup.rotation.step +1;
    dety = round(spot_list(i,12))+1;
    detz = round(spot_list(i,13))+1;
    
    % calculated spot center
    rot = spot_list(i,2);
    Omega = euler2u(rot*pi/180,0,0);
    voxelpos = Omega*pos';
    %diffraction center
    center = [L, voxelpos(2)*L/(Lsam2sou+voxelpos(1)), voxelpos(3)*L/(Lsam2sou+voxelpos(1))];
    alpha = atan(sqrt(voxelpos(2)^2+voxelpos(3)^2)/(Lsam2sou+voxelpos(1)));
    K_in = [Lsam2sou+voxelpos(1) voxelpos(2) voxelpos(3)];
    
    Gw = Omega*D^-1*S*U*B*spot_list(i,3:5)';
    
    beta = acos(dot(normr(K_in),normc(Gw)));
    theta = beta-pi/2;% diffraction angle
    sintth = sin(2*theta);
    
    v1 = [0 Gw(2) Gw(3)];
    phix= acos(dot(normr(K_in),normc(v1)));
    phiy = phix-2*theta;
    L2 = (Lsam2det-voxelpos(1))/cos(alpha);%norm(center-grainpos)
    difflength = L2*sintth/sin(phiy);% may need revision if the detector is not align perfectly perpendicular to the beam
    Gwnorm = normc(Gw);
        
    dety_theory = round(dety0+(center(2)+ difflength*Gwnorm(2))/pixelysize)+1;
    detz_theory = round(detz0+(center(3)+ difflength*Gwnorm(3))/pixelzsize)+1;

    if dety-peaksize>0 && dety+peaksize<detysize && detz-peaksize>0 && detz+peaksize<detzsize
        spot_segments = exp_imgs_flipped(detz-peaksize:detz+peaksize,dety-peaksize:dety+peaksize,img_num);
        [labeledImage, numBlobs] = bwlabel(spot_segments>0);
        spots = regionprops(labeledImage, 'Area','Centroid');
        
        if numBlobs==1 %&& sqrt((spots.Centroid(1)- ((size(labeledImage,1)-1)/2))^2 + (spots.Centroid(2)-(size(labeledImage,1)-1)/2)^2)<3%size(spots) == 1
            diff_imgs_flipped(detz-peaksize:detz+peaksize,dety-peaksize:dety+peaksize,img_num) = spot_segments;
            spot_list(i,16) = sqrt((spots.Centroid(1)- ((size(labeledImage,1)-1)/2))^2 + (spots.Centroid(2)-(size(labeledImage,1)-1)/2)^2);
        elseif numBlobs>1
            dist = zeros(10,1);
            
            for num_i = 1:numBlobs
                dist(num_i) = (spots(num_i).Centroid(1)- ((size(labeledImage,1)-1)/2))^2 + (spots(num_i).Centroid(2)-(size(labeledImage,1)-1)/2)^2;
            end
            [~,b] = min(dist(1:numBlobs));
            spot_id = labeledImage(round(spots(b).Centroid(2)),round(spots(b).Centroid(1)));
            spot_segments = (labeledImage == spot_id);
            %spot = regionprops(spot_segments, 'Area','Centroid');
            %if dist(b)<3 %&& sqrt((dety-dety_theory)^2+(detz-detz_theory)^2)<2
                diff_imgs_flipped(detz-peaksize:detz+peaksize,dety-peaksize:dety+peaksize,img_num) = spot_segments;
            %end 
            spot_list(i,16) = dist(b);
        end
    end
end

for i = 1:img_nums
    diff_imgs(:,:,i) = fliplr(flipud(diff_imgs_flipped(:,:,i)));
end

% imshow3D(diff_imgs)
