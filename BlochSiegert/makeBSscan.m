function [readoutDur, sampWin] = main(seq)
% Stack-of-spirals Bloch-Siegert B1 mapping scan
%
% Outputs:
%    readoutDur    spiral leaf duration (msec)
%    sampWin       [1 ndat]    sampling window (indeces) to be used in reconstruction
%
% Usage:
%  >> addpath ../..
%  >> seq = getparams();
%  >> main(seq);

%% Write modules.txt
mods.fatsat          = 'fatsat.mod';
mods.tipdown         = 'tipdown-norep.mod';  
mods.bs              = 'bs.mod';
mods.tipdownRephaser = 'tipdown-rephaser.mod';
mods.readout         = 'readout.mod';
mods.spoil           = 'spoiler.mod';
%mods.ex              = 'ex.mod';

loopFile = 'scanloop.txt';
moduleListFile = 'modules.txt';

modFileText = ['' ...
'Total number of unique cores\n' ...
sprintf('%d\n', length(fieldnames(mods))) ...
'fname	duration(us)	hasRF?	hasDAQ?\n' ...
sprintf('%s\t0\t1\t0\n'   , mods.fatsat          ), ...
sprintf('%s\t0\t1\t0\n'   , mods.tipdown         ), ...
sprintf('%s\t0\t1\t0\n'   , mods.bs              ), ...
sprintf('%s\t1000\t0\t0\n', mods.tipdownRephaser ), ...  
sprintf('%s\t0\t0\t1\n'   , mods.readout         ), ...
sprintf('%s\t0\t0\t0\n'   , mods.spoil           )];
%sprintf('%s\t0\t1\t0'     , mods.ex              )];
fid = fopen(moduleListFile, 'wt');
fprintf(fid, modFileText);
fclose(fid);

%% Create .mod files 

% fat saturation pulse
nCyclesSpoil = 1e-6;      % small enough so that no gradient trapezoid is added before RF pulse
toppe.utils.rf.makeslr(seq.fatsat.flip, seq.fatsat.slThick, seq.fatsat.tbw, seq.fatsat.dur, nCyclesSpoil, ...
    'type', 'ex', ...    % fatsat pulse is a 90 so is of type 'ex', not 'st' (small-tip)
    'ofname', mods.fatsat, 'system', seq.sys);

% tipdown-norep.mod and tipdown-rephaser.mod
if seq.rf.isHard
	% add spoiler before + after Bloch-Siegert RF pulse (bs.mod)
	dur = seq.rf.dur*1e-3;   % sec
	freq = 0;
	rftmp = toppe.utils.rf.makehardpulse(seq.rf.flip, dur);
	gtmp = toppe.utils.makecrusher(1, seq.stfr.res(1), 0, seq.sys.maxSlew/sqrt(2), seq.sys.maxGrad/sqrt(2));
	rf.tipdown = toppe.utils.makeGElength([  rftmp(:); 0*gtmp(:)]);
	g.tipdown  = toppe.utils.makeGElength([0*rftmp(:); gtmp(:)]);
	g.tipdownRephaser  = toppe.utils.makeGElength([-gtmp(:)]);
	toppe.writemod('rf', rf.tipdown, 'gz', g.tipdown, 'ofname', mods.tipdown, 'system', seq.sys);
	toppe.writemod('gz', g.tipdownRephaser, 'ofname', mods.tipdownRephaser, 'system', seq.sys);
else
	[~,~,freq] = toppe.utils.rf.makeslr(seq.bs.flip, seq.rf.slThick, seq.rf.tbw, seq.rf.dur, 0.001, ...
		'sliceOffset', seq.sliceOffset, 'type', 'st', ... 
		'ofname', 'tipdown.mod', 'forBlochSiegert', true, 'system', seq.sys);
	system('rm tipdown.mod');
end

% bs.mod
toppe.utils.rf.makebs(seq.bs.amp, 'system', seq.sys, 'ofname', mods.bs);

% stack-of-spirals readout
[g,roInfo] = toppe.utils.spiral.makesosreadout(seq.bs.fov, seq.bs.matrix, ...
	seq.bs.nLeafs, seq.sys.maxSlew, 'system', seq.sys, 'ofname', mods.readout);
readoutDur = seq.sys.raster*1e3*length(roInfo.sampWin);  % msec
sampWin = roInfo.sampWin;
save sampWin

% spoiler
gspoil = toppe.utils.makecrusher(seq.bs.nCyclesSpoil, seq.bs.res(1), 0, 0.5*seq.sys.maxSlew/sqrt(2), seq.sys.maxGrad/sqrt(2));
toppe.writemod('ofname', mods.spoil, 'gx', gspoil, 'gz', gspoil, 'system', seq.sys);

% excitation for double-angle b1 measurement
if 0
	toppe.utils.rf.makeslr(90, seq.rf.slThick, seq.rf.tbw, 1.5*seq.rf.dur, 0.001, ...
		'spoilDerate', 0.8, ...
		'type', 'ex', ...
		'ofname', mods.ex, 'system', seq.sys);
end

%% Write scanloop.txt

% get min TR so we can calculate textra (to achieve desired TR)
toppe.write2loop('setup');
toppe.write2loop(mods.fatsat);
toppe.write2loop(mods.spoil);
toppe.write2loop(mods.tipdown);
toppe.write2loop(mods.bs);
toppe.write2loop(mods.tipdownRephaser);
toppe.write2loop(mods.readout);
toppe.write2loop(mods.spoil);
toppe.write2loop('finish');
mintr = toppe.getTRtime(1,7)*1e3;  % ms

nz = seq.bs.matrix(3);

% Initial setup
rfphs = 0;
rf_spoil_seed = 117;
rf_spoil_seed_cnt = 1;
toppe.write2loop('setup', 'version', 3);

% Bloch-Siegert scan
for iim = 1:2
	fprintf('.');

	bsfreq = 4000 * (-1)^(iim-1);   % Fermi pulse frequency offset (Hz). See toppe.utils.rf.makebs

   for iz = -0:nz             % iz < 1 are discarded acquisitions (to reach steady state)
		for ileaf = 1:seq.bs.nLeafs
         if iz > 0
        		a_gz = ((iz-1+0.5)-nz/2)/(nz/2);   % z phase-encode amplitude, scaled to (-1,1) range
         else
         	a_gz = 0; 
         end

			% fat sat
  	 		toppe.write2loop(mods.fatsat, 'RFphase', rfphs, 'Gamplitude', [0 0 0]', 'RFoffset', seq.fatsat.freqOffset);
   		toppe.write2loop(mods.spoil, 'rotmat', seq.rotmat);

	   	% rf excitation
	  		toppe.write2loop(mods.tipdown, 'RFphase', rfphs, ...
				'rotmat', seq.rotmat, ...
				'RFoffset', round(freq));

			% Bloch-Siegert pulse
  			toppe.write2loop(mods.bs, 'RFoffset', bsfreq, 'RFphase', rfphs);

			% Slice-select gradient rephaser
  			toppe.write2loop(mods.tipdownRephaser, ...
				'rotmat', seq.rotmat);

   		% readout. Data is stored in 'slice', 'echo', and 'view' indeces.
			phi = 2*pi*(ileaf-1)/seq.bs.nLeafs;              % leaf rotation angle (radians)
			slice = max(iz,1);
			echo = iim;
			view = ileaf;
   		toppe.write2loop('readout.mod', 'DAQphase', rfphs, ...
				'slice', slice, 'echo', echo, 'view', view, ...
				'rotmat', seq.rotmat, ...        % prescribed scan plane rotation (logical frame)
				'rot',    phi, ...               % in-plane rotation (in logical frame)
				'dabmode', 'on', 'Gamplitude', [1 1 a_gz]');

			% spoiler
   		toppe.write2loop(mods.spoil, 'textra', max(seq.bs.tr - mintr, 0), ...
				'rotmat', seq.rotmat);

		   % update rf phase (RF spoiling)
			rfphs = rfphs + (rf_spoil_seed/180 * pi)*rf_spoil_seed_cnt ;  % radians
			rf_spoil_seed_cnt = rf_spoil_seed_cnt + 1;
      end
	end
end

% Double-angle b1 mapping
if 0
for iz = 1:-1 %1:nz
	for ileaf = 1:seq.bs.nLeafs
		for iim = 3:4
			if iim == 3
				a_rf = 0.5;
			else
				a_rf = 1;
			end

			a_gz = ((iz-1+0.5)-nz/2)/(nz/2);   % z phase-encode amplitude, scaled to (-1,1) range

			% fat sat
  			toppe.write2loop(mods.fatsat, ...
				'RFphase', rfphs, ...
				'Gamplitude', [0 0 0]', ...
				'RFoffset', seq.fatsat.freqOffset);
  			toppe.write2loop(mods.spoil, 'rotmat', seq.rotmat);

		   % update rf phase (RF spoiling)
			rfphs = rfphs + (rf_spoil_seed/180 * pi)*rf_spoil_seed_cnt ;  % radians
			rf_spoil_seed_cnt = rf_spoil_seed_cnt + 1;

  		 	% rf excitation
  			toppe.write2loop(mods.ex, ...
				'RFphase', rfphs, ...
				'RFamplitude', a_rf, ...
				'rotmat', seq.rotmat, ...
				'RFoffset', round(freq));

  		 	% readout. Data is stored in 'slice', 'echo', and 'view' indeces.
			phi = 2*pi*(ileaf-1)/seq.bs.nLeafs;              % leaf rotation angle (radians)
			slice = max(iz,1);
			echo = iim;
			view = ileaf;
  		 	toppe.write2loop('readout.mod', ...
				'DAQphase', rfphs, ...
				'Gamplitude', [1 1 a_gz]', ...
				'slice', slice, 'echo', echo, 'view', view, ...
				'rotmat', seq.rotmat, ...        % prescribed scan plane rotation (logical frame)
				'rot',    phi, ...               % in-plane rotation (in logical frame)
				'dabmode', 'on');

			% spoiler
  		 	toppe.write2loop(mods.spoil, ...
				'rotmat', seq.rotmat, ... 
				'textra', 5000);

		   % update rf phase (RF spoiling)
			rfphs = rfphs + (rf_spoil_seed/180 * pi)*rf_spoil_seed_cnt ;  % radians
			rf_spoil_seed_cnt = rf_spoil_seed_cnt + 1;
		end
	end
end
end

% finish
fprintf('\n');
toppe.write2loop('finish');

%% create tar file
system('tar czf ~/tmp/scan,mwf,bs,spiral.tgz ../../getparams.m main.m *.mod modules.txt scanloop.txt sampWin.mat');

return;
