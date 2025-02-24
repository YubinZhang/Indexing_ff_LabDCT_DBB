function Gvs = spotpos2gvector(spots,parameters,pos)

dety0 = parameters.detector.dety0;
detz0 = parameters.detector.detz0;
pixelysize = parameters.detector.pixelysize;
pixelzsize = parameters.detector.pixelzsize;
Rtilt = parameters.detector.tilt;

Lsam2sou = parameters.setup.Lss;
Lsam2det = parameters.setup.Lsd;

K_in =  [pos(1) pos(2) pos(3)] - Lsam2sou;

% if size(spots,1)>1
%     disp('spot overlap');
% end
for i = 1:size(spots,1)
    dety = spots(i).WeightedCentroid(1); % flip is considered already during spot finding
    detz = spots(i).WeightedCentroid(2);
    spot_pos = [Lsam2det 0 0]'+ Rtilt*[0 (dety-dety0)*pixelysize (detz-detz0)*pixelzsize]';
    K_out = spot_pos - pos;
    %K_out = [Lsam2det-pos(1) (dety-dety0)*pixelysize-pos(2) (detz-detz0)*pixelzsize-pos(3)];
    Gvs(i,:) = normr(normr(K_out')-normr(K_in));
end
