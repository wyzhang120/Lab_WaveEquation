%% unpack from saved file
imgDir='C:\DFiles\Geophysics\Project\Figs_Crosswell';
fmod1 = 'vp22_pad_smooth.mat';
mod1 = load(fullfile(imgDir, fmod1));
if mod1.fastz; vel = mod1.vp'; else; vel = mod1.vp; end
[nx, nz] = size(vel); dx = mod1.dx; dz = mod1.dz;
x = (0:nx-1)*dx; z = (0:nz-1)*dz;
% source and receiver geometry
sx = mod1.xsrc; sz = mod1.zsrc; gx = mod1.xrec; gz = mod1.zrec; 
g=1:numel(gx); 
%% DF params
nbc=40; nt=2001; dt=0.0005; t=(0:nt-1)*dt; isFS=false; 
freq=300; s=ricker(freq,dt);
