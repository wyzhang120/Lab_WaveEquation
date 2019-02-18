%% unpack from saved file
imgDir='C:\DFiles\Geophysics\Project\Figs_Crosswell';
% path to save figure
figDir = fullfile(pwd, 'Figures');

fmod1 = 'vp22_pad_smooth.mat';
mod1 = load(fullfile(imgDir, fmod1));
if mod1.fastz; vel = mod1.vp'; else; vel = mod1.vp; end
[nz, nx] = size(vel); dx = mod1.dx; dz = mod1.dz;
x = (0:nx-1)*dx; z = (0:nz-1)*dz;
% source and receiver geometry
idsrc=53; idxTrace = idsrc;  
sx = mod1.xsrc(idsrc); sz = mod1.zsrc(idsrc); gx = mod1.xrec; gz = mod1.zrec; 
g=1:numel(gx); 
%% FD params
velMax = max(max(vel)); velMin = min(min(vel)); 
freq = 300; % central frequnecy of Ricker
fmax = 1.5*freq;
% spatial constraint of FD 2D
%    dh <= vmin/(n fmax), 
%    where n is the FD order, fmax is max source frequncy
% tempal constraint
%    dt <= dh / ( h sqrt(2) vmax)
%    where h is constant related to FD order, 8th order h = 2161/1680

nbc=40; nt=double(mod1.nt); dt=mod1.dt; t=(0:nt-1)*dt; isFS=false; 
s=ricker(mod1.fc,dt);
tic; seis=a2d_mod_abc28(vel,nbc,dx,nt,dt,s,sx,sz,gx,gz,isFS); toc;
%% plotting
% set path to save figures
if ~exist('figDir', 'var')
    figDir = fullfile(pwd, 'Figures');
    if ~exist(figDir, 'dir')
        mkdir(figDir);
    end
end

figure(1);ax1=subplot(221);
imagesc(x,z,vel);colorbar;caxis([3000 6000]);colormap(ax1, flipud(jet));
xlabel('X (m)'); ylabel('Z (m)'); title('velocity');
figure(1);ax2=subplot(222);
imagesc(gz,t,seis);colormap(ax2, gray);ylabel('Time (t)');ylim([0 0.15]);
xlabel('X(m)'); title(sprintf('seismic gather src (z=%.2f m)', sz));
% select trace and time window

figure(1);subplot(223);
plot(t,seis(:,idxTrace));xlabel('Time (t)'); 
title(sprintf('seismic trace (sz=%.2f m, gz = %.2f m)', sz, gz(idxTrace)));
idxMute = round(0.066/dt);
trace=seis(:,idxTrace); trace(idxMute:end)=0;
tic; wp=a2d_wavepath_abc28(seis(:,idxTrace),vel,nbc,dx,nt,dt,s,sx,sz,gx(idxTrace),gz(idxTrace)); toc;
figure(1);ax4=subplot(224);
imagesc(x,z,wp);colormap(ax4, gray); caxis([-5e-4 5e-4]);
% use histogram to set scale; figure(2);histogram(reshape(wp, [numel(wp) 1]));
xlabel('X (m)'); ylabel('Z (m)'); title(sprintf('wave path (sz=%.2f m, gz = %.2f m)', sz, gz(idxTrace)));
hold on;
plot(sx, sz, '*r', 'LineWidth', 1, 'MarkerSize', 10);
plot(gx(idxTrace),gz(idxTrace), '<g', 'MarkerFaceColor', 'g', 'MarkerSize', 5);
hold off;
savefig(fullfile(figDir, sprintf('xwell_sz%.fm_gz%.fm.fig', sz, gz(idxTrace))));