% clear all
% clear lambda
% load('C:\Users\adli\Documents\LabDCT_datasets\2019_10_03_Ref_Iron_2019-10-03_170704\Results.mat')

cell =[2.8665,2.8665,2.8665,90,90,90]; %crystallographic parameters
B=FormB(cell);
Amatrix=2*pi*inv(B);

S=[1 0 0;0 -1 0;0 0 1]; %microscope system is left handed

% sample to source/detector distances, only used for fitting detector
Lsam2sou=100;
Lsam2det=100;

% expected strain, only used for statistical purposes (finding error)
epsilon=[0,0,0;0,0,0;0,0,0];
% epsilon=[0,1e-3,0;1e-3,0,0;0,0,0];
% epsilon=[0,0,1e-3;0,0,0;1e-3,0,0];
% epsilon=[0,0,0;0,0,1e-3;0,1e-3,0];
% epsilon=[1e-3,0,0;0,-0.5e-3,0;0,0,-0.5e-3];

% flags for what to fit
fit_detector=0;
fit_detector_rot=0;
fit_COM=1;
grainCOM=zeros(size(sampos));


options = optimoptions('fmincon','Display','off'); % removes clutter from command window

if fit_detector==1
    Spotpos_init=SpotInfo(:,9:11);
    Spotpos_simu=[Lsam2det*ones(size(SpotInfo(:,1))),SpotInfo(:,16:17)];
    
%     funaff=@(x0,y0,z0,a1,a2,a3,x,y,z)[x;y;z]+[cos(a3) sin(a3) 0; -sin(a3) cos(a3) 0; 0 0 1]*...
%         [cos(a2) 0 -sin(a2); 0 1 0; sin(a2) 0 cos(a2)]*[1 0 0; 0 cos(a1) sin(a1); 0 -sin(a1) cos(a1)]*[x0;y0;z0];
    if fit_detector_rot==1
        disp("fitting detector position with respect to y,z, and rotation about x")
        % function defining the affine transformation
        funaff=@(x0,y0,z0,a1,y,z)ones(size([x0,y0,z0])).*[0,y,z]+[x0,y0,z0]*[1 0 0; 0 cos(a1) sin(a1); 0 -sin(a1) cos(a1)];
        % objective function, difference between simulated and real spot COM
        funSpotpos=@(x)sum(sum(abs(funaff(Spotpos_init(:,1),Spotpos_init(:,2),Spotpos_init(:,3),x(1),x(2),x(3))-Spotpos_simu),2));
        Spotpos_shift=fmincon(funSpotpos,[0,0,0],[],[],[],[],[],[],[],options);
        Spotpos_fit=funaff(Spotpos_init(:,1),Spotpos_init(:,2),Spotpos_init(:,3),Spotpos_shift(1),Spotpos_shift(2),Spotpos_shift(3));
    else
        disp("fitting detector position with respect to y,z")
        funaff=@(x0,y0,z0,y,z)ones(size([x0,y0,z0])).*[0,y,z]+[x0,y0,z0];
        funSpotpos=@(x)sum(sum(abs(funaff(Spotpos_init(:,1),Spotpos_init(:,2),Spotpos_init(:,3),x(1),x(2))-Spotpos_simu),2));
        Spotpos_shift=fmincon(funSpotpos,[0,0],[],[],[],[],[],[],[],options);
        Spotpos_fit=funaff(Spotpos_init(:,1),Spotpos_init(:,2),Spotpos_init(:,3),Spotpos_shift(1),Spotpos_shift(2));
    end
    
%     funSpotpos=@(x)sum(abs(sqrt(funaff(Spotpos_init(:,1),Spotpos_init(:,2),Spotpos_init(:,3),x(1),x(2),x(3))-Spotpos_simu)),'All');
%     funSpotpos=@(x)sum(abs(sqrt((funaff(Spotpos_init(:,1),Spotpos_init(:,2),Spotpos_init(:,3),x(1),x(2))-Spotpos_simu).^2)),'All');
%     funSpotpos=@(x)sum(abs(sqrt(sum((funaff(Spotpos_init(:,1),Spotpos_init(:,2),Spotpos_init(:,3),x(1),x(2))-Spotpos_simu),2).^2)),'All');
%     funSpotpos=@(x)sum(sum(abs(funaff(Spotpos_init(:,1),Spotpos_init(:,2),Spotpos_init(:,3),x(1),x(2),x(3))-Spotpos_simu),2));
%     funSpotpos=@(x)sum(sum(abs(funaff(Spotpos_init(:,1),Spotpos_init(:,2),Spotpos_init(:,3),x(1),x(2))-Spotpos_simu),2));
    
    
%     
    
    Spotpos_diff_init_fit=abs(Spotpos_init-Spotpos_fit);
    Spotpos_diff_fit_simu=abs(Spotpos_fit-Spotpos_simu);
    Spotpos_diff_init_simu=abs(Spotpos_init-Spotpos_simu);
else
    
end

if fit_COM==1
    disp("fitting grain center of mass positions")
for grainno=1:grains
% for grainno = [5 10 17 35 42]
% for grainno=[4 5 9 42]
    % associate spot data with grains
    GrainInfo{grainno}=SpotInfo(find(SpotInfo(:,2)==grainno),:);
    G{grainno}=GrainInfo{grainno}(:,3:5);
    H{grainno}=GrainInfo{grainno}(:,6:8);
    lambda{grainno}=GrainInfo{grainno}(:,15);
    Grainpos{grainno}=GrainInfo{grainno}(:,12:14);
    if fit_detector==1
        Spotpos{grainno}=Spotpos_fit(find(SpotInfo(:,2)==grainno),:);
    else
        Spotpos{grainno}=GrainInfo{grainno}(:,9:11);
    end
    rot=GrainInfo{grainno}(:,1);
    % function to calculate G vector from spot and grain position
    Gfun=@(p,x)2*pi./(12.398./(lambda{grainno})).*(p-x)./sqrt(sum((p-x).^2,2))-2*pi./(12.398./(lambda{grainno})).*(x+[Lsam2sou 0 0])./sqrt(sum((x+[Lsam2sou 0 0]).^2,2));
    
    %objective function, sum of angles between simualated and measured G
    %vectors as a function of grain center of mass
    funCOM=@(x)abs(sum(acos(dot(((iOmegafun3(rot,Gfun(Spotpos{grainno},Omegafun4(rot,x,S)),eye(3))) ...
        ./sqrt(sum(((iOmegafun3(rot,Gfun(Spotpos{grainno},Omegafun4(rot,x,S)),eye(3)))).^2,2)))', ...
        (((S*U_simu{grainno}*B*H{grainno}')')./sqrt(sum(((S*U_simu{grainno}*B*H{grainno}')').^2,2)))'))));
    
%     funCOM=@(x)abs(sum(acos(dot((2*pi./(12.398./(lambda{grainno})).*(Spotpos{grainno}-Omegafun3(Omega_tot,x).*ones(size(Spotpos{grainno})) ... %K_out
%     ./sqrt(sum((Spotpos{grainno}-Omegafun3(Omega_tot,x).*ones(size(Spotpos{grainno}))).^2,2))) ... %abs(K_out)
%     -2*pi./(12.398./(lambda{grainno})).*((Omegafun3(Omega_tot,x)+[Lsam2sou 0 0])./sqrt(sum((Omegafun3(Omega_tot,x)+[Lsam2sou 0 0]).^2,2)))) ... %K_in/abs(K_in)
%     ./sqrt(sum((2*pi./(12.398./(lambda{grainno})).*(Spotpos{grainno}-Omegafun3(Omega_tot,x).*ones(size(Spotpos{grainno})) ... %abs(Gvec)
%     ./sqrt(sum((Spotpos{grainno}-Omegafun3(Omega_tot,x).*ones(size(Spotpos{grainno}))).^2,2)))-2*pi./(12.398./(lambda{grainno})).* ... 
%     ((Omegafun3(Omega_tot,x)+[Lsam2sou 0 0])./sqrt(sum((Omegafun3(Omega_tot,x)+[Lsam2sou 0 0]).^2,2)))),2).^2) ...
%     ,(S*U_simu{grainno}*B*H{grainno}'./sqrt(sum((S*U_simu{grainno}*B*H{grainno}').^2,2)))')))); % simulated Gvec
    disp("grain " + num2str(grainno))
    grainCOM(grainno,:)=fmincon(funCOM,sampos(grainno,:),[],[],[],[],[],[],[],options);
    gCOMdiff=sampos-grainCOM;
end
else
    grainCOM=sampos;
end

disp("fitting strain and orientation")
for grainno=1:grains
% for grainno=[5 10 17 35 42]
% for grainno=[4 5 9 42]
    % associate spot data with grains
    GrainInfo{grainno}=SpotInfo(find(SpotInfo(:,2)==grainno),:);
    G{grainno}=GrainInfo{grainno}(:,3:5);
    H{grainno}=GrainInfo{grainno}(:,6:8);
    lambda{grainno}=GrainInfo{grainno}(:,15);
    Grainpos{grainno}=GrainInfo{grainno}(:,12:14);
    if fit_detector==1
        Spotpos{grainno}=Spotpos_fit(find(SpotInfo(:,2)==grainno),:);
    else
        Spotpos{grainno}=GrainInfo{grainno}(:,9:11);
    end
    rot=GrainInfo{grainno}(:,1);
%     if fit_COM==1
%         Kout=2*pi./(12.398./(lambda{grainno})).*(Spotpos{grainno}-Omegafun4(rot,grainCOM(grainno,:),S));
%         Kin=2*pi./(12.398./(lambda{grainno})).*(Omegafun4(rot,grainCOM(grainno,:),S)+[Lsam2sou 0 0]);
%     else
%         Kout=2*pi./(12.398./(lambda{grainno})).*(Spotpos{grainno}-Grainpos{grainno})./sqrt(sum((Spotpos{grainno}-Grainpos{grainno}).^2,2));
%         Kin=2*pi./(12.398./(lambda{grainno})).*(Grainpos{grainno}+[Lsam2sou 0 0])./sqrt(sum((Grainpos{grainno}+[Lsam2sou 0 0]).^2,2));
%     end
%     Gvec=Kout-Kin;
%     Gvec_w=zeros(size(Gvec));
%     for i = 1:length(rot)
%         omega=rot(i)*pi/180;
%         Omega = [cos(omega) -sin(omega) 0;sin(omega) cos(omega) 0;0 0 1];
%         Gvec_w(i,:)=(Omega\Gvec(i,:)')';
%     end
%     funeps=@(x)abs(sum(acos(dot(((Kout-Kin)./sqrt(sum((Kout-Kin).^2,2)))',S*x*H{grainno}'./sqrt(sum((S*x*H{grainno}').^2,2))))));
%     funeps=@(x)abs(sum(acos(dot(((Gvec)./sqrt(sum((Gvec).^2,2))),Omegafun4(Omega_tot,x*H{grainno}'./sqrt(sum((S*x*H{grainno}').^2,2)))))));
    
    % function to calculate G vector from spot and grain position
    Gfun=@(p,x)2*pi./(12.398./(lambda{grainno})).*(p-x)./sqrt(sum((p-x).^2,2))-2*pi./(12.398./(lambda{grainno})).*(x+[Lsam2sou 0 0])./sqrt(sum((x+[Lsam2sou 0 0]).^2,2));
    
    % function to rotate G vector
    Gvec_w=iOmegafun3(rot,Gfun(Spotpos{grainno},Omegafun4(rot,grainCOM(grainno,:),S)),eye(3));
    
%     funeps=@(x)abs(sum(acos(dot(Gvec_w'./sqrt(sum(G{grainno}'.^2,1)),S*x*eye(3)*H{grainno}'./sqrt(sum((S*x*eye(3)*H{grainno}').^2))))));
    
    % objective function, sum of angles between simualated and measured G
    % vectors as a function of U*B
    funeps=@(x)abs(sum(acos(dot(Gvec_w'./sqrt(sum(Gvec_w'.^2,1)),S*x*eye(3)*H{grainno}'./sqrt(sum((S*x*eye(3)*H{grainno}').^2))))));
    
    
    disp("grain " + num2str(grainno))
    [Q,R]=qr(fmincon(funeps,U_simu{grainno}*B,[],[],[],[],[],[],@(x)nodilat(x,B),options));
    % result of fitting QR-decomposed into B and U with positive diagonals
    D=diag(sign(diag(R)));
    B_exp{grainno}=D*R;
    U_exp{grainno}=Q*D;
    
    A_exp=eye(3)/B_exp{grainno}*2*pi;
    T=A_exp/Amatrix;
    eps_measured{grainno}=1/2*(T+T')-eye(3);
    eps_rotated{grainno}=inv(U_exp{grainno})*eps_measured{grainno}*U_exp{grainno};
    epsilon_rotated{grainno}=inv(U_simu{grainno})*epsilon*U_simu{grainno};
    rotated_diff{grainno}=eps_rotated{grainno}-epsilon_rotated{grainno};
end
statistics_epsilon

% figure
% hold on
% plot(Spotpos_simu(:,2),Spotpos_simu(:,3),'r*')
% plot(Spotpos_init(:,2),Spotpos_init(:,3),'b*')
% for i=1:length(Spotpos_init(:,1))
%     plot([Spotpos_init(i,2) Spotpos_simu(i,2)], [Spotpos_init(i,3) Spotpos_simu(i,3)])
% end
% 
% figure
% hold on
% plot(Spotpos_simu(:,2),Spotpos_simu(:,3),'r*')
% plot(Spotpos_fit(:,2),Spotpos_fit(:,3),'b*')
% for i=1:length(Spotpos_fit(:,1))
%     plot([Spotpos_fit(i,2) Spotpos_simu(i,2)], [Spotpos_fit(i,3) Spotpos_simu(i,3)])
% end