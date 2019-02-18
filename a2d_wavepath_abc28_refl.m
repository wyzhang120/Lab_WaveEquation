function wp=a2d_wavepath_abc28_refl(seis,v,refl,nbc,dx,nt,dt,s,sx,sz,gx,gz)
%  Copyright (C) 2010 Center for Subsurface Imaging and Fluid Modeling (CSIM),
%  King Abdullah University of Science and Technology, All rights reserved.
%
%  author:   Xin Wang
%  email:    xin.wang@kaust.edu.sa
%  date:     Sep 26, 2012
%  purpose:  2DTDFD solver to calculate wavepath with a born modeling ,
%            accuracy of 2-8, use the absorbing boundary condition
%
%  IN   seis(:,:) -- seismogram,    v(:,:) -- velocity,  refl(:,:) -- reflectivity
%       nbc       -- grid number of boundary
%       dx        -- grid intervel, nt          -- number of sample
%       dt        -- time interval, s(:)        -- wavelet
%       sx,sz     -- src position,  gx(:),gz(:) -- rec position
%  OUT  wp(:,:)   -- Output wavepath image

[nz,nx]=size(v); ng=numel(gx); wp=zeros(nz,nx);

c1 = -205.0/72.0; c2 = 8.0/5.0; c3 = -1.0/5.0; c4 = 8.0/315.0; c5 = -1.0/560.0;

% setup ABC and temperary variables
v=padvel(v,nbc); refl=padvel(refl,nbc);
abc=AbcCoef2D(v,nbc,dx);
alpha=(v*dt/dx).^2; kappa=abc*dt;
temp1=2+2*c1*alpha-kappa; temp2=1-kappa;
beta_dt = (v*dt).^2;
s=expand_source(s,nt);

[isx,isz,igx,igz]=adjust_sr(sx,sz,gx,gz,dx,nbc);

p0=zeros(size(v)); p1=zeros(size(v));
q0=zeros(size(v)); q1=zeros(size(v));

pwf=zeros(nz+2,nx+2,nt);
qwf=zeros(nz+2,nx+2,nt);

% pbc_top=zeros(5,nx,nt);
% pbc_bottom=zeros(5,nx,nt);
% pbc_left=zeros(nz,5,nt);
% pbc_right=zeros(nz,5,nt);
% pbc_nt=zeros(size(v));
% pbc_nt_1=zeros(size(v));
% 
% qbc_top=zeros(5,nx,nt);
% qbc_bottom=zeros(5,nx,nt);
% qbc_left=zeros(nz,5,nt);
% qbc_right=zeros(nz,5,nt);
% qbc_nt=zeros(size(v));
% qbc_nt_1=zeros(size(v));

for it=1:nt
    p=temp1.*p1-temp2.*p0+alpha.*...
        (c2*(circshift(p1,[0,1,0])+circshift(p1,[0,-1,0])+circshift(p1,[1,0,0])+circshift(p1,[-1,0,0]))...
        +c3*(circshift(p1,[0,2,0])+circshift(p1,[0,-2,0])+circshift(p1,[2,0,0])+circshift(p1,[-2,0,0]))...
        +c4*(circshift(p1,[0,3,0])+circshift(p1,[0,-3,0])+circshift(p1,[3,0,0])+circshift(p1,[-3,0,0]))...
        +c5*(circshift(p1,[0,4,0])+circshift(p1,[0,-4,0])+circshift(p1,[4,0,0])+circshift(p1,[-4,0,0])));
    p(isz,isx) = p(isz,isx) + beta_dt(isz,isx) * s(it);
    % Save Wave field
    pwf(:,:,it)=p(nbc:nbc+nz+1,nbc:nbc+nx+1);
    % save BC
    %[pbc_top(:,:,it),pbc_bottom(:,:,it),pbc_left(:,:,it),pbc_right(:,:,it)]=save_boundary(p,nz,nx,nbc);
    q=temp1.*q1-temp2.*q0+alpha.*...
        (c2*(circshift(q1,[0,1,0])+circshift(q1,[0,-1,0])+circshift(q1,[1,0,0])+circshift(q1,[-1,0,0]))...
        +c3*(circshift(q1,[0,2,0])+circshift(q1,[0,-2,0])+circshift(q1,[2,0,0])+circshift(q1,[-2,0,0]))...
        +c4*(circshift(q1,[0,3,0])+circshift(q1,[0,-3,0])+circshift(q1,[3,0,0])+circshift(q1,[-3,0,0]))...
        +c5*(circshift(q1,[0,4,0])+circshift(q1,[0,-4,0])+circshift(q1,[4,0,0])+circshift(q1,[-4,0,0])));
    q=pertubation(p1,q,refl,beta_dt);
    % Save Wave Field
    qwf(:,:,it)=q(nbc:nbc+nz+1,nbc:nbc+nx+1);
    % Save BC
    %[qbc_top(:,:,it),qbc_bottom(:,:,it),qbc_left(:,:,it),qbc_right(:,:,it)]=save_boundary(q,nz,nx,nbc);
    p0=p1; p1=p;
    q0=q1; q1=q;
end

% Time Loop
p1(nbc:nbc+nz+1,nbc:nbc+nx+1)=pwf(:,:,nt);
q1(nbc:nbc+nz+1,nbc:nbc+nx+1)=qwf(:,:,nt);

bp1=zeros(size(v));bp0=zeros(size(v));
bq1=zeros(size(v));bq0=zeros(size(v));

for it=nt-2:-1:1
    p(nbc:nbc+nz+1,nbc:nbc+nx+1)=pwf(:,:,it);
    q(nbc:nbc+nz+1,nbc:nbc+nx+1)=qwf(:,:,it);
    bq=temp1.*bq1-temp2.*bq0+alpha.*...
        (c2*(circshift(bq1,[0,1,0])+circshift(bq1,[0,-1,0])+circshift(bq1,[1,0,0])+circshift(bq1,[-1,0,0]))...
        +c3*(circshift(bq1,[0,2,0])+circshift(bq1,[0,-2,0])+circshift(bq1,[2,0,0])+circshift(bq1,[-2,0,0]))...
        +c4*(circshift(bq1,[0,3,0])+circshift(bq1,[0,-3,0])+circshift(bq1,[3,0,0])+circshift(bq1,[-3,0,0]))...
        +c5*(circshift(bq1,[0,4,0])+circshift(bq1,[0,-4,0])+circshift(bq1,[4,0,0])+circshift(bq1,[-4,0,0])));
    % Add seismogram
    for ig=1:ng
        bq(igz(ig),igx(ig))=bq(igz(ig),igx(ig))+beta_dt(igz(ig),igx(ig))*seis(it,ig);
    end
    bp=temp1.*bp1-temp2.*bp0+alpha.*...
        (c2*(circshift(bp1,[0,1,0])+circshift(bp1,[0,-1,0])+circshift(bp1,[1,0,0])+circshift(bp1,[-1,0,0]))...
        +c3*(circshift(bp1,[0,2,0])+circshift(bp1,[0,-2,0])+circshift(bp1,[2,0,0])+circshift(bp1,[-2,0,0]))...
        +c4*(circshift(bp1,[0,3,0])+circshift(bp1,[0,-3,0])+circshift(bp1,[3,0,0])+circshift(bp1,[-3,0,0]))...
        +c5*(circshift(bp1,[0,4,0])+circshift(bp1,[0,-4,0])+circshift(bp1,[4,0,0])+circshift(bp1,[-4,0,0])));
    bp=pertubation(bq1,bp,refl,beta_dt);
    wp=image_condition(wp,p1,bp0,nz,nx,nbc);
    wp=image_condition(wp,q1,bq0,nz,nx,nbc);
    % wf refresh
    p1=p; bp0=bp1; bp1=bp;
    q1=q; bq0=bq1; bq1=bq;
end

end