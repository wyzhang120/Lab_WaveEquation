function [img,illum]=a2d_rtm_abc28_snapshot(seis,v,nbc,dx,nt,dt,s,sx,sz,gx,gz)
%  Copyright (C) 2010 Center for Subsurface Imaging and Fluid Modeling (CSIM),
%  King Abdullah University of Science and Technology, All rights reserved.
%
%  author:   Xin Wang
%  email:    xin.wang@kaust.edu.sa
%  date:     Sep 26, 2012
%  purpose:  2DTDFD reverse time migration with accuracy of 2-8
%            use the absorbing boundary condition
%
%  IN   seis(:,:) -- seismogram,    v(:,:) -- velocity
%       nbc       -- grid number of boundary
%       dx        -- grid intervel, nt          -- number of sample
%       dt        -- time interval, s(:)        -- wavelet
%       sx,sz     -- src position,  gx(:),gz(:) -- rec position
%  OUT  img(:,:)           -- Output reverse time migration image
%       illum(:,:)         --        illumination compensation

[nz,nx]=size(v); ng=numel(gx); img=zeros(nz,nx); 
x = (0:nx-1)*dx; z = (0:nz-1)*dx;

c1 = -205.0/72.0; c2 = 8.0/5.0; c3 = -1.0/5.0; c4 = 8.0/315.0; c5 = -1.0/560.0;

% setup ABC and temperary variables
v=padvel(v,nbc);
abc=AbcCoef2D(v,nbc,dx);
alpha=(v*dt/dx).^2; kappa=abc*dt;
temp1=2+2*c1*alpha-kappa; temp2=1-kappa;
beta_dt = (v*dt).^2;
s=expand_source(s,nt);

[isx,isz,igx,igz]=adjust_sr(sx,sz,gx,gz,dx,nbc);

p1=zeros(size(v)); p0=zeros(size(v)); q=zeros(size(v)); bp=zeros(size(v));

bc_top=zeros(5,nx,nt);
bc_bottom=zeros(5,nx,nt);
bc_left=zeros(nz,5,nt);
bc_right=zeros(nz,5,nt);

figure(3); subplot(234);imagesc(x,z,bp(nbc+1:nbc+nz,nbc+1:nbc+nx));
figure_title=['Reconstruct Forward Field, t=',num2str((nt-1)*dt),' s'];
title(figure_title);ylabel('Z (m)');xlabel('X (M)');caxis([-0.25 0.25]);
figure(3); subplot(235);imagesc(x,z,q(nbc+1:nbc+nz,nbc+1:nbc+nx));
figure_title=['Backward Propagate Field, t=',num2str((nt-1)*dt),' s'];
title(figure_title);ylabel('Z (m)');xlabel('X (M)');caxis([-0.25 0.25]);
figure(3); subplot(236);imagesc(x,z,img);
figure_title=['RTM image, t=',num2str((nt-1)*dt),' s'];
title(figure_title);ylabel('Z (m)');xlabel('X (M)');caxis([-0.25 0.25]);
% Time Looping
for it=1:nt
    p=temp1.*p1-temp2.*p0+alpha.*...
        (c2*(circshift(p1,[0,1,0])+circshift(p1,[0,-1,0])+circshift(p1,[1,0,0])+circshift(p1,[-1,0,0]))...
        +c3*(circshift(p1,[0,2,0])+circshift(p1,[0,-2,0])+circshift(p1,[2,0,0])+circshift(p1,[-2,0,0]))...
        +c4*(circshift(p1,[0,3,0])+circshift(p1,[0,-3,0])+circshift(p1,[3,0,0])+circshift(p1,[-3,0,0]))...
        +c5*(circshift(p1,[0,4,0])+circshift(p1,[0,-4,0])+circshift(p1,[4,0,0])+circshift(p1,[-4,0,0])));
    p(isz,isx) = p(isz,isx) + beta_dt(isz,isx) * s(it);
    % dipole source
    %p(isz-2,isx) = p(isz-2,isx) - beta_dt(isz-2,isx) *wavelet(it);
    if mod(it,100)==1
        figure(3);subplot(233);colormap(gray);imagesc(x,z,p(nbc+1:nbc+nz,nbc+1:nbc+nx));
        figure_title=['Forward Propagate Field, t=',num2str((it-1)*dt),' s'];
        title(figure_title);ylabel('Z (m)');xlabel('X (M)');caxis([-0.25 0.25]); pause(0.2);
    end
    % save BC
    [bc_top(:,:,it),bc_bottom(:,:,it),bc_left(:,:,it),bc_right(:,:,it)]=save_boundary(p,nz,nx,nbc);
    p0=p1;
    p1=p;
end
% save final wavefield
bc_p_nt_1=p0;
bc_p_nt=p1;

bp1=bc_p_nt_1;
bp0=bc_p_nt;

illum=bp1(nbc+1:nbc+nz,nbc+1:nbc+nx).^2+bp0(nbc+1:nbc+nz,nbc+1:nbc+nx).^2;

q0=zeros(size(v)); q1=zeros(size(v));

% Time Loop
for it=nt-2:-1:1
    % subtrace source
    bp0(isz,isx)=bp0(isz,isx)-s(it+2)*beta_dt(isz,isx);
    % dipole source
    %p0(isz-2,isx)=p0(isz-2,isx)-wavelet(it+2)*beta_dt(isz-2,isx);
    bp=temp1.*bp1-temp2.*bp0+alpha.*...
        (c2*(circshift(bp1,[0,1,0])+circshift(bp1,[0,-1,0])+circshift(bp1,[1,0,0])+circshift(bp1,[-1,0,0]))...
        +c3*(circshift(bp1,[0,2,0])+circshift(bp1,[0,-2,0])+circshift(bp1,[2,0,0])+circshift(bp1,[-2,0,0]))...
        +c4*(circshift(bp1,[0,3,0])+circshift(bp1,[0,-3,0])+circshift(bp1,[3,0,0])+circshift(bp1,[-3,0,0]))...
        +c5*(circshift(bp1,[0,4,0])+circshift(bp1,[0,-4,0])+circshift(bp1,[4,0,0])+circshift(bp1,[-4,0,0])));
    bp=load_boundary(bp,bc_top(:,:,it),bc_bottom(:,:,it),bc_left(:,:,it),bc_right(:,:,it),nz,nx,nbc);
    
    q=temp1.*q1-temp2.*q0+alpha.*...
        (c2*(circshift(q1,[0,1,0])+circshift(q1,[0,-1,0])+circshift(q1,[1,0,0])+circshift(q1,[-1,0,0]))...
        +c3*(circshift(q1,[0,2,0])+circshift(q1,[0,-2,0])+circshift(q1,[2,0,0])+circshift(q1,[-2,0,0]))...
        +c4*(circshift(q1,[0,3,0])+circshift(q1,[0,-3,0])+circshift(q1,[3,0,0])+circshift(q1,[-3,0,0]))...
        +c5*(circshift(q1,[0,4,0])+circshift(q1,[0,-4,0])+circshift(q1,[4,0,0])+circshift(q1,[-4,0,0])));
    % Add seismogram
    for ig=1:ng
        q(igz(ig),igx(ig))=q(igz(ig),igx(ig))+beta_dt(igz(ig),igx(ig))*seis(it,ig);
    end
    img=image_condition(img,bp1,q0,nz,nx,nbc);
    if mod(it,100)==1
        figure(3); subplot(234);imagesc(x,z,bp(nbc+1:nbc+nz,nbc+1:nbc+nx));
        figure_title=['Reconstruct Forward Field, t=',num2str((it-1)*dt),' s'];
        title(figure_title);ylabel('Z (m)');xlabel('X (M)');caxis([-0.25 0.25]);
        figure(3); subplot(235);imagesc(x,z,q(nbc+1:nbc+nz,nbc+1:nbc+nx));
        figure_title=['Backward Propagate Field, t=',num2str((it-1)*dt),' s'];
        title(figure_title);ylabel('Z (m)');xlabel('X (M)');caxis([-1 1]);
        figure(3); subplot(236);imagesc(x,z,img);
        figure_title=['RTM image, t=',num2str((it-1)*dt),' s'];
        title(figure_title);ylabel('Z (m)');xlabel('X (M)');caxis([-5 5]);pause(0.2);
    end
    illum=illum+bp(nbc+1:nbc+nz,nbc+1:nbc+nx).^2;
    
    % wf refresh
    bp0=bp1; bp1=bp;
    q0=q1; q1=q;
end
end