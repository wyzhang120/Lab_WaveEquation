function d_out=highpass(d_in,nit1,nit2,nit3)
%  Copyright (C) 2010 Center for Subsurface Imaging and Fluid Modeling (CSIM),
%  King Abdullah University of Science and Technology, All rights reserved.
%
%  author:   Xin Wang
%  email:    xin.wang@kaust.edu.sa
%  date:     Sep 26, 2012
%  purpose:  highpass filter for 1/2/3D
%
%  IN   d_in(:,:,:) -- input data
%       nit1        -- number of iterations on the first direction
%       nit2        -- number of iterations on the second direction
%       nit3        -- number of iterations on the third direction
%  OUT  d_out(:,:,:)   -- Output filtered data

nd=ndims(d_in); d_in0=d_in;
if nd==3
    if nargin==2
        if nit1==0
            d_out=d_in;
        else
            nit2=nit1; nit3=nit1;
            d_in=padvel(d_in,1);
            for it=1:nit1
                d_temp=0.25*(circshift(d_in,[-1,0,0])+circshift(d_in,[1,0,0]))+0.5*d_in;
                d_in=d_temp;
            end
            for it=1:nit2
                d_temp=0.25*(circshift(d_in,[0,-1,0])+circshift(d_in,[0,1,0]))+0.5*d_in;
                d_in=d_temp;
            end
            for it=1:nit3
                d_temp=0.25*(circshift(d_in,[0,0,-1])+circshift(d_in,[0,0,1]))+0.5*d_in;
                d_in=d_temp;
            end
            d_out=d_in0-d_in(2:end-1,2:end-1);
        end
    elseif nargin==4
        d_in=padvel(d_in,1);
        for it=1:nit1
            d_temp=0.25*(circshift(d_in,[-1,0,0])+circshift(d_in,[1,0,0]))+0.5*d_in;
            d_in=d_temp;
        end
        for it=1:nit2
            d_temp=0.25*(circshift(d_in,[0,-1,0])+circshift(d_in,[0,1,0]))+0.5*d_in;
            d_in=d_temp;
        end
        for it=1:nit3
            d_temp=0.25*(circshift(d_in,[0,0,-1])+circshift(d_in,[0,0,1]))+0.5*d_in;
            d_in=d_temp;
        end
        d_out=d_in0-d_in(2:end-1,2:end-1);
    else
        error('Not Enough input arguments')
    end
elseif nd==2
    [n1,n2]=size(d_in);
    if min(n1,n2)==1 % 1D
        if nit1==0
            d_out=d_in;
        else
            d_in=d_in(:); d_in=[d_in(1);d_in;d_in(end)];
            for it=1:nit1
                d_temp=0.25*(circshift(d_in,[-1,0,0])+circshift(d_in,[1,0,0]))+0.5*d_in;
                d_in=d_temp;
            end
            d_out=d_in0-d_in(2:end-1,2:end-1);
            d_out=reshape(d_out,n1,n2);
        end
    else
        if nargin==2
            if (nit1==0)
                d_out=d_in;
            else
                nit2=nit1;
                d_in=padvel(d_in,1);
                for it=1:nit1
                    d_temp=0.25*(circshift(d_in,[-1,0,0])+circshift(d_in,[1,0,0]))+0.5*d_in;
                    d_in=d_temp;
                end
                for it=1:nit2
                    d_temp=0.25*(circshift(d_in,[0,-1,0])+circshift(d_in,[0,1,0]))+0.5*d_in;
                    d_in=d_temp;
                end
                d_out=d_in0-d_in(2:end-1,2:end-1);
            end
        else
            d_in=padvel(d_in,1);
            for it=1:nit1
                d_temp=0.25*(circshift(d_in,[-1,0,0])+circshift(d_in,[1,0,0]))+0.5*d_in;
                d_in=d_temp;
            end
            for it=1:nit2
                d_temp=0.25*(circshift(d_in,[0,-1,0])+circshift(d_in,[0,1,0]))+0.5*d_in;
                d_in=d_temp;
            end
            d_out=d_in0-d_in(2:end-1,2:end-1);
        end
    end
end

end