function makeBSscan
% Fully sampled stack-of-spirals Bloch-Siegert B1 mapping scan
%
% Outputs:
%    readoutDur    spiral leaf duration (msec)
%    sampWin       [1 ndat]    sampling window (indeces) to be used in reconstruction
%
% Usage:
%  >> addpath ../..
%  >> seq = getparams();
%  >> main(seq);

%% Sequence parameters

FOV = [24 24 20];      % cm
N = [60 60 50];        % matrix size
res = FOV./N;          % voxel size

tr = 150;              % ms (to reduce SAR)
nCyclesSpoil = 2;

bs.amp = 0.05;            % Amplitude of Fermi pulse (Gauss)
bs.freq = 4000;           % offset frequency (Hz)

rf.flip = 10;             % low flip angle for spin-density weighting
rf.slThick = 0.8*FOV(3);  % avoid wrap-around
rf.tbw = 8;
rf.dur = 1;               % ms
rf.type = 'st';           % small-tip design
rf.ftype = 'min';         % minimum-phase design is ok for 3D imaging

fatsat.flip    = 90;
fatsat.slThick = 1000;       % dummy value (determines slice-select gradient, but we won't use it; just needs to be large to reduce dead time before+after rf pulse)
fatsat.tbw     = 2.5;        % time-bandwidth product
fatsat.dur     = 5;          % pulse duration (ms)
fatsat.freqOffset = -440;    % Hz

spiral.F0 = FOV(1);
spiral.F1 = 0;
spiral.F2 = 0;
spiral.nLeafs = 4;
spiral.dir = 1;   % spiral-out


%% System hardware specs
% It is recommended that maxRF match the physical scanner limits, 
% to ensure accurate B1 scaling.
% maxGrad and maxSlew can be < physical limit
sys = toppe.systemspecs('maxSlew', 9, 'slewUnit', 'Gauss/cm/ms', ...
    'maxGrad', 5, 'gradUnit', 'Gauss/cm', ...
    'myrfdel', 152, ... % (us) best if multiple of 4us
    'daqdel', 152, ...  % (us) best if multiple of 4us
    'timessi', 200);    % us


%% Write modules.txt
mods.fatsat          = 'fatsat.mod';
mods.tipdown         = 'tipdown-norep.mod';  
mods.bs              = 'bs.mod';
mods.tipdownRephaser = 'tipdown-rephaser.mod';
mods.partEncode      = 'partEncode.mod'; % z partition-encoding gradient (at full strength)
mods.readout         = 'readout.mod';    % 2D spiral leafs
mods.spoiler         = 'spoiler.mod';

fid = fopen('modules.txt', 'wt');
fprintf(fid, 'Total number of unique cores\n');
fprintf(fid, '%d\n', length(fieldnames(mods)));
fprintf(fid, 'fname  duration(us)    hasRF?  hasDAQ?\n');
fprintf(fid, '%s\t0\t1\t0\n'   , mods.fatsat          );
fprintf(fid, '%s\t0\t1\t0\n'   , mods.tipdown         );
fprintf(fid, '%s\t0\t1\t0\n'   , mods.bs              );
fprintf(fid, '%s\t0\t0\t0\n'   , mods.tipdownRephaser ); 
fprintf(fid, '%s\t0\t0\t0\n'   , mods.partEncode      ); 
fprintf(fid, '%s\t0\t0\t1\n'   , mods.readout         );
fprintf(fid, '%s\t0\t0\t0\n'   , mods.spoiler         );
fclose(fid);


%% Create .mod files 

% fat sat pulse
b1 = toppe.utils.rf.makeslr(fatsat.flip, fatsat.slThick, fatsat.tbw, fatsat.dur, 1e-6, sys, ...
    'type', 'ex', ...    % fatsat pulse is a 90 so is of type 'ex', not 'st' (small-tip)
    'writeModFile', false);
b1 = toppe.makeGElength(b1);
toppe.writemod(sys, 'rf', b1, 'ofname', mods.fatsat, 'desc', 'fat sat pulse');

% tipdown-norep.mod and tipdown-rephaser.mod
% We do this since the BS pulse goes before the rephaser.
[~,~,freq] = toppe.utils.rf.makeslr(rf.flip, rf.slThick, rf.tbw, rf.dur, 0.001, sys, ...
    'type', rf.type, ... 
    'ftype', rf.ftype, ...
    'ofname', 'tipdown.mod', ...
    'forBlochSiegert', true);
system('rm tipdown.mod');

% Bloch-Siegert pulse
toppe.utils.rf.makebs(bs.amp, 'system', sys, 'ofname', mods.bs);

% spiral leafs (balanced, since we don't want net spoiling area to rotate)
if false
%[gx, gy, gz, paramsint16, paramsfloat] = toppe.utils.spiral.makeVDSreadout(...
[gx, gy, gz, paramsint16, paramsfloat] = toppe.utils.spiral.makeVDSreadout(...
    FOV(1), spiral.F0, spiral.F1, spiral.F2, N(1), spiral.nLeafs, ...
    sys.maxSlew/sqrt(2), spiral.dir, sys, ...
    'gmax', sys.maxGrad, 'bal', true);
end

% partition-encoding (z) gradient 
gamma = 4257.6;       % Hz/G
raster = 4e-3;          % gradient/ADC sample duration (ms)
zres = FOV(3)/N(3);     % cm
kmax = 1/(2*zres);       % cycles/cm
areaPartEncode = 2*kmax/gamma;    % G/cm*s
gzPartEncode= -toppe.utils.trapwave2(areaPartEncode/2, sys.maxGrad, sys.maxSlew, raster);
gzPartEncode = toppe.makeGElength(gzPartEncode(:));
toppe.writemod(sys, 'gz', gzPartEncode, 'ofname', mods.partEncode, ...
    'desc', 'partition-encoding gradient');

% spoiler
gspoiler = toppe.utils.makecrusher(nCyclesSpoil, res(1), sys, 0, 0.5*sys.maxSlew/sqrt(2), sys.maxGrad/sqrt(2));
gspoiler = toppe.makeGElength(gspoiler);
toppe.writemod(sys, 'ofname', mods.spoiler, 'gx', gspoiler, 'gz', gspoiler);


%% Get min TR so we can calculate delay to achieve desired TR
toppe.write2loop('setup', sys);
toppe.write2loop(mods.fatsat, sys);
toppe.write2loop(mods.spoiler, sys);
toppe.write2loop(mods.tipdown, sys);
toppe.write2loop(mods.bs, sys);
toppe.write2loop(mods.tipdownRephaser, sys);
toppe.write2loop(mods.partEncode, sys);
toppe.write2loop(mods.readout);
toppe.write2loop(mods.partEncode, sys);
toppe.write2loop(mods.spoiler, sys);
toppe.write2loop('finish', sys);
mintr = toppe.getTRtime(1,9,sys)*1e3;  % ms


%% Write scanloop.txt
nz = N(3);

% Initial setup
rfphs = 0;
rf_spoil_seed = 117;
rf_spoil_seed_cnt = 0;
toppe.write2loop('setup', sys);

for iim = 1:2
    fprintf('.');

    bsf = bs.freq* (-1)^(iim-1);   % toggle frequency 

    for iz = -0:nz             % iz < 1 are discarded acquisitions (to reach steady state)
        for ileaf = 1:spiral.nLeafs
            if iz > 0
                a_gz = ((iz-1+0.5)-nz/2)/(nz/2);   % z phase-encode amplitude, scaled to (-1,1) range
            else
                a_gz = 0; 
            end

            % fat sat
            toppe.write2loop(mods.fatsat, sys, ...
                'RFphase', rfphs, ...
                'Gamplitude', [0 0 0]', ...
                'RFoffset', fatsat.freqOffset);
            toppe.write2loop(mods.spoiler, sys);

            % rf excitation
            toppe.write2loop(mods.tipdown, sys, ...
                'RFphase', rfphs);

            % Bloch-Siegert pulse
            toppe.write2loop(mods.bs, sys, ...
                'RFoffset', bsf, 'RFphase', rfphs);

            % Slice-select gradient rephaser
            toppe.write2loop(mods.tipdownRephaser, sys);

            % Partition-encode gradient
            toppe.write2loop(mods.partEncode, sys, ...
                'Gamplitude', [0 0 a_gz]');

            % readout. Data is stored in 'slice', 'echo', and 'view' indeces.
            slice = max(iz,1);
            echo = iim;
            view = ileaf;
            toppe.write2loop('readout.mod', sys, ...
                'DAQphase', rfphs, ...
                'slice', slice, 'echo', echo, 'view', view, ...
                'waveform', ileaf, ...
                'dabmode', 'on');

            % Partition-encode rephaser
            toppe.write2loop(mods.partEncode, sys, ...
                'Gamplitude', [0 0 -a_gz]');

            % spoiler
            toppe.write2loop(mods.spoiler, sys, ...
                'textra', max(tr - mintr, 0));

            % update rf phase (RF spoiling)
            rfphs = rfphs + (rf_spoil_seed/180 * pi)*rf_spoil_seed_cnt ;  % radians
            rf_spoil_seed_cnt = rf_spoil_seed_cnt + 1;
        end
    end
end
fprintf('\n');
toppe.write2loop('finish', sys);

return;

%% create tar file
system('tar czf ~/tmp/scan,mwf,bs,spiral.tgz ../../getparams.m main.m *.mod modules.txt scanloop.txt sampWin.mat');

return;
