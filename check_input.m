%check_input
if exist('Lsam2det')
    parameters.setup.Lsd = Lsam2det;
else
    parameters.setup.Lsd = 100;
end

if exist('Lsam2sou')
    parameters.setup.Lss = Lsam2sou;
else
    parameters.setup.Lsd = [-100 0 0];
end

if exist('rot_start')
    parameters.setup.rotation.start = rot_start;
else
    parameters.setup.rotation.start = -180;
end

if exist('rot_step')
    parameters.setup.rotation.step = rot_step;
else
    parameters.setup.rotation.step = 2;
end

if exist('rot_end')
    parameters.setup.rotation.end = rot_end;
else
    parameters.setup.rotation.end = 0;
end

if exist('binning')
    parameters.detector.binning = binning;
else
     parameters.detector.binning = 4;
end

if exist('detysize')
    parameters.detector.detysize = detysize/parameters.detector.binning;
else
     parameters.detector.detysize = 2048/parameters.detector.binning;
end

if exist('detzsize')
    parameters.detector.detzsize = detzsize/parameters.detector.binning;
else
     parameters.detector.detzsize = 2048/parameters.detector.binning;
end

if exist('pixelysize')
    parameters.detector.pixelysize = pixelysize*parameters.detector.binning;
else
     parameters.detector.pixelysize = 0.033588*parameters.detector.binning;
end

if exist('pixelzsize')
    parameters.detector.pixelzsize = pixelzsize*parameters.detector.binning;
else
     parameters.detector.pixelzsize = 0.033588*parameters.detector.binning;
end

if exist('dety0')
    parameters.detector.dety0 = dety0/parameters.detector.binning;
else
     parameters.detector.dety0 = 1024/parameters.detector.binning;
end

if exist('detz0')
    parameters.detector.detz0 = detz0/parameters.detector.binning;
else
     parameters.detector.detz0 = 1024/parameters.detector.binning;
end

if exist('tilt')
    parameters.detector.tilt = tilt;
else
     parameters.detector.tilt = eye(3);
end

if exist('tilt_y')
    parameters.detector.tilt_y = tilt_y;
else
     parameters.detector.tilt_y = 0;
end

if exist('tilt_z')
    parameters.detector.tilt_z = tilt_z;
else
     parameters.detector.tilt_z = 0;
end

if exist('polfactor')
    parameters.beam.polfactor = polfactor;
else
    parameters.beam.polfactor = 1;
end

if exist('Lorentz')
    parameters.beam.Lorentz = Lorentz;
else
    parameters.beam.Lorentz = 1;
end

if exist('I0')
    parameters.beam.I0 = I0;
else
    parameters.beam.I0 = 1e13;
end

if exist('K1')
    parameters.beam.K1 = K1;
end

if exist('K2')
    parameters.beam.K2 = K2;
end

if exist('E_max')
    parameters.beam.E_max = E_max;
else 
    parameters.beam.E_max = 120;
end

if exist('E_min')
    parameters.beam.E_min = E_min;
else 
    parameters.beam.E_min = 15;
end

if exist('beamstop')
    parameters.detector.beamstop = round(beamstop/binning);
else 
    parameters.detector.beamstop = [700 1300 700 1300]/binning;
end

if exist('peakshape')
    parameters.detector.peakshape = peakshape;
else 
    parameters.detector.peakshape = 0;
end

if exist('peaksize')
    parameters.detector.peaksize = peaksize;
else 
    parameters.detector.peaksize = 5;
end
