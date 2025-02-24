%
% Parameters describing experimental setup 
%
% Experimental setup
Lsam2det = 100%99.2125		% sample-detector distance (mm)
Lsam2sou = [-100 0 0]%109.8628 %109.2740
rot_start = 180;
rot_step = -3;
rot_end = -180;

image_num = (rot_end-rot_start)/rot_step+1;
L = Lsam2det-Lsam2sou;

%detector 
binning = 2
dety0 = 512.4836*2+1%1024		% beamcenter, y in pixel coordinatees%512.5787 
detz0 = 512.0230*2+1%1024		% beamcenter, z in pixel coordinatees%511.9820
pixelysize = 0.034167%3588 		% Pixel size y (mm)
pixelzsize = 0.034167%3588	% Pixel size z (mm)
detysize = 2048            % detector y size (pixels)
detzsize = 2048             % detector z size (pixels)
beamstop = [840 1220 840 1220]%[700 1300 700 1300]%

tilt   = eye(3)           % detector tilt counterclockwise around lab x axis in rad 
% tilt = [0.999999996504214,8.15483617871188e-05,-1.84779904920508e-05;-8.15506191331780e-05,0.999999989208788,-0.000122196241111896;1.84680253893717e-05,0.000122197747576289,0.999999992363321];
% [0.999999992973457,0.000114728771289443,-2.98394744957016e-05;-0.000114732917493999,0.999999983759660,-0.000138985745717046;2.98235283472655e-05,0.000138989168310423,0.999999989896284]
% polfactor = 0               % Do not apply polarisation factor
polfactor = 1               % Apply polarisation factor
% Lorentz = 0                 % Do not apply Lorentz factor
Lorentz = 1                 % Apply Lorentz factor
beampol = 1.0               % Beam polarisation (1 = totally horizontally polarised)
                            %                   (0 = nonpolarised beam)
I0 = 1e13                   % Beam flux (Ph/s/mm2)
thetamax = 6.0              % Maximum theta angle for reflection generation

% Generate and calc SF or read list from file
readhkl = 0                 % generate reflections
%readhkl = 1                % read reflections from file

%Do we make the frames
save_indexing_pattern = 1;          % Do not make tif frames

direc = 'Indexing_results'
prefix = 'Iron_ff_LabDCT_' % prefix of  frame names default is 'frame'
mkdir 'Indexing_results'           % save frames in this directory 

% Define peak shape if diffraction images are to be formed 
peakshape = 0;          % Represent peak as a spike 3x3 pixels  
peaksize = 5;

E_max = 150;
E_min = 15;

lambda_min = 12.39818746/150;
Al = 0;
Fe = 1;
if Fe
    input_fe
    % cell = [2.8665 2.8665 2.8665 90 90 90];
    B = FormB(cell);
    V = cellvolume(cell);
    emass =9.1093826e-31;
    echarge = 1.60217653e-19;
    pi4eps0 = 1.11265e-10;
    c = 299792458;
    K1 =  (echarge^2/(pi4eps0*emass*c*c)*1000)^2; % Unit is mm
    grainsize = 20;
    grainvolume = grainsize*grainsize*grainsize*pi/6;
    K2 = lambda_min*lambda_min*lambda_min*grainvolume*10^12/V^2;
    load('Ahkl_Fe_new.mat');
    % load('Ahkl_Al.mat');
    AAhkl = Ahkl(:,1).^2+Ahkl(:,2).^2+Ahkl(:,3).^2;
    ind1 = find(AAhkl==8);
    ind2 = find(AAhkl==16);
    Ahkl([ind1;ind2],:) = [];
    Ahkl_indexing = Ahkl(sum(Ahkl(:,1:3).^2,2)<=6,:);
    Ahkl_output = Ahkl(sum(Ahkl(:,1:3).^2,2)<=6,:);
elseif Al
    input_al
    B = FormB(cell);
    V = cellvolume(cell);
    emass =9.1093826e-31;
    echarge = 1.60217653e-19;
    pi4eps0 = 1.11265e-10;
    c = 299792458;
    K1 =  (echarge^2/(pi4eps0*emass*c*c)*1000)^2; % Unit is mm
    grainsize = 20;
    grainvolume = grainsize*grainsize*grainsize*pi/6;
    K2 = lambda_min*lambda_min*lambda_min*grainvolume*10^12/V^2;
    %load('Ahkl_Fe_new.mat');
    load('Ahkl_Al.mat');
    AAhkl = Ahkl(:,1).^2+Ahkl(:,2).^2+Ahkl(:,3).^2;
    Ahkl_indexing = Ahkl(sum(Ahkl(:,1:3).^2,2)<=11,:);
    Ahkl_output = Ahkl(sum(Ahkl(:,1:3).^2,2)<=11,:);
end
%sample
Pos = [0 0 0];
% PosList = [0 0 0; 0 0.150 0; 0.150 0 0;0 0 0.150;...
%            0.150 0.150 0;0 0.150 0.150; 0.150 0 0.150;0.150 0.150 0.150;...
%            0 -0.150 0; -0.150 0 0;0 0 -0.150;...
%            -0.150 -0.150 0;0 -0.150 -0.150; -0.150 0 -0.150;...
%            0.150 -0.150 0;0 0.150 -0.150; 0.150 0 -0.150;...
%            -0.150 0.150 0;0 -0.150 0.150; -0.150 0 0.150;...
%            0.150 -0.150 0.150;0.150 0.150 -0.150;-0.150 0.150 0.150;
%            -0.150 -0.150 0.150;-0.150 0.150 -0.150;0.150 -0.150 -0.150;-0.150 -0.150 -0.150];
% PosList = unique(PosList,'rows');
fitdetector = 0;
reindex = 0;
check_result = 0;
checking_angle  = 0.15;