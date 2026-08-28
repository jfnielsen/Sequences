% writeB0.m
%
% Create the file 'b0.seq', a 3D spoiled GRE B0 mapping sequence.

% RF/gradient delay (sec). 
% Conservative choice that should work across all GE scanners.
psd_rf_wait = 200e-6;  
addpath(genpath("pulseq/"))
addpath(genpath("B0shimming/"))
% System/design parameters.
% Here we extend rfRingdownTime by psd_rf_wait to ensure that 
% the subsequent wait pulse (delay block) doesn't overlap with
% the 'true' RF ringdown time (54 us).
sys = mr.opts('maxGrad', 30, 'gradUnit','mT/m', ...
              'maxSlew', 100, 'slewUnit', 'T/m/s', ...
              'gradRasterTime', 4e-6, ...
              'rfRasterTime', 4e-6, ...
              'rfDeadTime', 100e-6, ...
              'rfRingdownTime', 60e-6 + psd_rf_wait, ...
              'adcDeadTime', 40e-6, ...
              'adcRasterTime', 2e-6, ...
              'blockDurationRaster', 4e-6, ...
              'B0', 3.0);

% Acquisition parameters
fov = [240e-3 240e-3 240e-3];   % FOV (m)
Nx = 80; Ny = Nx; Nz = 80;      % Matrix size
dwell = 12e-6;                  % ADC sample time (s)
alpha = 3;                      % flip angle (degrees)
fatChemShift = 3.5e-6;          % 3.5 ppm
fatOffresFreq = sys.gamma*sys.B0*fatChemShift;  % Hz
TE = 1/fatOffresFreq*[1 2];     % fat and water in phase for both echoes
TR = 8e-3*[1 1];                % constant TR
nCyclesSpoil = 2;               % number of spoiler cycles
Tpre = 1.0e-3;                  % prephasing trapezoid duration
rfSpoilingInc = 117;            % RF spoiling increment

% Create a new sequence object
seq = mr.Sequence(sys);           

% Create non-selective pulse
[rf] = mr.makeBlockPulse(alpha/180*pi, sys, 'Duration', 0.2e-3, 'use', 'excitation');

% Define other gradients and ADC events
% Cut the redaout gradient into two parts for optimal spoiler timing
deltak = 1./fov;
Tread = Nx*dwell;

gyPre = mr.makeTrapezoid('y', sys, ...
    'Area', Ny*deltak(2)/2, ...   % maximum PE1 gradient, max positive amplitude
    'Duration', Tpre);
gzPre = mr.makeTrapezoid('z', sys, ...
    'Area', Nz*deltak(3)/2, ...   % maximum PE2 gradient, max positive amplitude
    'Duration', Tpre);

gx = mr.makeTrapezoid('x', sys, ...  
    'Amplitude', Nx*deltak(1)/Tread, ...
    'FlatTime', Tread);
gxPre = mr.makeTrapezoid('x', sys, ...
    'Area', -gx.area/2, ...
    'Duration', Tpre);

adc = mr.makeAdc(Nx, sys, ...
    'Duration', Tread,...
    'Delay', gx.riseTime);

gzSpoil = mr.makeTrapezoid('z', sys, ...
    'Area', Nx*deltak(1)*nCyclesSpoil);
gxSpoil = mr.makeTrapezoid('x', sys, ...
    'Area', Nx*deltak(1)*nCyclesSpoil);
%gxSpoil = mr.makeTrapezoid('x', sys, ...
%    'Area', Nx*deltak(1)*nCyclesSpoil);

% y/z PE steps
pe1Steps = ((0:Ny-1)-Ny/2)/Ny*2;
pe2Steps = ((0:Nz-1)-Nz/2)/Nz*2;

% Calculate timing
TEmin = rf.shape_dur/2 + rf.ringdownTime + mr.calcDuration(gxPre) ...
      + adc.delay + Nx/2*dwell;
delayTE = ceil((TE-TEmin)/seq.gradRasterTime)*seq.gradRasterTime;
TRmin = mr.calcDuration(rf) + delayTE + mr.calcDuration(gxPre) ...
      + mr.calcDuration(gx) + mr.calcDuration(gxSpoil);
delayTR = ceil((TR-TRmin)/seq.gradRasterTime)*seq.gradRasterTime;

% Loop over phase encodes and define sequence blocks
% iZ < 0: Dummy shots to reach steady state
% iZ = 0: ADC is turned on and used for receive gain calibration on GE scanners
% iZ > 0: Image acquisition

nDummyZLoops = 1;
rf_phase = 0;
rf_inc = 0;

textprogressbar('Generating sequence object: ');
for iZ = -nDummyZLoops:Nz
    textprogressbar(max(iZ,1)/Nz*100);

    isDummyTR = iZ < 0;

    for iY = 1:Ny
        % Turn on y and z prephasing lobes, except during dummy scans and
        % receive gain calibration (auto prescan)
        yStep = (iZ > 0) * pe1Steps(iY);
        zStep = (iZ > 0) * pe2Steps(max(1,iZ));

        if yStep == 0
            yStep = eps;
        end
        if zStep == 0
            zStep = eps;
        end

        for c = 1:length(TE)

            % RF spoiling
            rf.phaseOffset = rf_phase/180*pi;
            adc.phaseOffset = rf_phase/180*pi;
            rf_inc = mod(rf_inc+rfSpoilingInc, 360.0);
            rf_phase = mod(rf_phase+rf_inc, 360.0);
            
            % Excitation
            % Mark start of segment (block group) by adding label.
            % Subsequent blocks in block group are NOT labelled.
            seq.addBlock(rf, mr.makeLabel('SET', 'TRID', 2-isDummyTR));
            
            % Encoding
            seq.addBlock(mr.makeDelay(delayTE(c)));
            seq.addBlock(gxPre, ...
                mr.scaleGrad(gyPre, yStep), ...
                mr.scaleGrad(gzPre, zStep));
            if isDummyTR
                seq.addBlock(gx);
            else
                seq.addBlock(gx, adc);
            end

            % rephasing/spoiling
            seq.addBlock(gxSpoil, ...
                mr.scaleGrad(gyPre, -yStep), ...
                mr.scaleGrad(gzPre, -zStep));
            seq.addBlock(mr.makeDelay(delayTR(c)));
        end
    end
end
textprogressbar('');
fprintf('Sequence ready\n');

% Check sequence timing
[ok, error_report]=seq.checkTiming;
if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end

% Output for execution
seq.setDefinition('FOV', fov);
seq.setDefinition('Name', 'b0');
seq.write([seq_name '.seq']);

% Plot sequence
Noffset = length(TE)*Ny*(nDummyZLoops+1);
seq.plot('timerange',[Noffset Noffset+4]*TR(1), 'timedisp', 'ms');
