clear;
%% unpack from saved file
% imgDir = 'C:\DFiles\Geophysics\Project\Figs_Crosswell';
imgDir = 'E:/Geophysics/Project/Crosswell/FWI_2arr';
rayDir = 'E:\Geophysics\Inversion\Shuster_labs\Raytracinglab';
addpath(rayDir);
% path to save figure
figDir = fullfile(pwd, 'Figures');

fmod1 = 'vp22_pad_smooth.mat';
mod1 = load(fullfile(imgDir, fmod1));
if mod1.fastz 
    vel = mod1.vpGaussBlur';
    vel0 = mod1.vp';
else 
    vel = mod1.vpGaussBlur;
    vel0 = mod1.vp;
end
sln = 1./vel;
refl_ss=(sln-1./vel0)./sln;
[nz, nx] = size(vel); dx = mod1.dx; dz = mod1.dz;
x = (0:nx-1)*dx; z = (0:nz-1)*dz;
% source and receiver geometry
idsrc=53; idxTrace = idsrc;  
sx = mod1.xsrc(idsrc); sz = mod1.zsrc(idsrc); gx = mod1.xrec; gz = mod1.zrec; 
g=1:numel(gx);
%% FD params
nbc=40; nt=double(mod1.nt); dt=mod1.dt; t=(0:nt-1)*dt; isFS=false; 
s=ricker(mod1.fc,dt);
csg0 =CSG(sx, sz, gx, gz);
csg0 = csg0.getSeis(vel,nbc,dx,nt,dt,mod1.fc,isFS);
csg0 = csg0.getFirstArrival();
csg0.plotCSG('FA', 'y');
csg0.plotTrace(idxTrace);
t1=0.0608; t2=0.066;
wp1 = csg0.getWP(idxTrace, 'tStart', t1, 'tEnd', t2, 'wavetype', 'refl', 'reflCoeff', refl_ss);
hd_wp1=wp1.plotWP();
% wp0.savefig(hdl);
% tic; seis=a2d_mod_abc28(vel,nbc,dx,nt,dt,s,sx,sz,gx,gz,isFS); toc;
% %% plotting
% % set path to save figures
% if ~exist('figDir', 'var')
%     figDir = fullfile(pwd, 'Figures');
%     if ~exist(figDir, 'dir')
%         mkdir(figDir);
%     end
% end
% 
% figure(1);ax1=subplot(221);
% imagesc(x,z,vel);colorbar;caxis([3000 6000]);colormap(ax1, flipud(jet));
% xlabel('X (m)'); ylabel('Z (m)'); title('velocity');
% figure(1);ax2=subplot(222);
% imagesc(gz,t,seis);colormap(ax2, gray);ylabel('Time (t)');ylim([0 0.15]);
% xlabel('X(m)'); titllle(sprintf('seismic gather src (z=%.2f m)', sz));
% % select trace and time window
% figure(1);subplot(223);
% plot(t,seis(:,idxTrace));xlabel('Time (t)'); 
% title(sprintf('seismic trace (sz=%.2f m, gz = %.2f m)', sz, gz(idxTrace)));
% idxMute = round(0.066/dt);
% trace=seis(:,idxTrace); trace(idxMute:end)=0;
% tic; wp=a2d_wavepath_abc28(trace,vel,nbc,dx,nt,dt,s,sx,sz,gx(idxTrace),gz(idxTrace)); toc;
% figure(1);ax4=subplot(224);
% imagesc(x,z,wp);colormap(ax4, gray); caxis([-5e-4 5e-4]);
% % use histogram to set scale; figure(2);histogram(reshape(wp, [numel(wp) 1]));
% xlabel('X (m)'); ylabel('Z (m)'); 
% title(sprintf('wave path (sz=%.2f m, gz = %.2f m)', sz, gz(idxTrace)));
% hold on;
% plot(sx, sz, '*r', 'LineWidth', 1, 'MarkerSize', 10);
% plot(gx(idxTrace),gz(idxTrace), '<g', 'MarkerFaceColor', 'g', 'MarkerSize', 5);
% hold off;
% savefig(fullfile(figDir, sprintf('xwell_sz%.fm_gz%.fm.fig', sz, gz(idxTrace))));