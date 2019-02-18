function [isx,isz,igx,igz]=adjust_sr(sx,sz,gx,gz,dx,nbc)
%  Copyright (C) 2010 Center for Subsurface Imaging and Fluid Modeling (CSIM),
%  King Abdullah University of Science and Technology, All rights reserved.
%
%  author:   Xin Wang
%  email:    xin.wang@kaust.edu.sa
%  date:     Sep 26, 2012
%  purpose:  adjust the source and receiver location to avoid putting on
%            the free surface
%
%  IN   nbc      -- grid number of boundary, 
%       dt       -- time interval, dx -- grid intervel
%       sx,sz    -- src position,  gx(:),gz(:) -- rec position
%  OUT  isx,isz  -- Output source location on the grid
%       igx(:),igz(:)      receiver location on the grid

% set and adjust the free surface position
isx=round(sx/dx)+1+nbc;isz=round(sz/dx)+1+nbc;
igx=round(gx/dx)+1+nbc;igz=round(gz/dx)+1+nbc;
if abs(sz) <0.5, isz=isz+1; end
igz=igz+(abs(gz)<0.5)*1;
end