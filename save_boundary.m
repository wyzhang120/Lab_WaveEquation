function [top,bottom,left,right]=save_boundary(p,nz,nx,nbc)
top=p(nbc:nbc+4,nbc+1:nbc+nx);
bottom=p(nz+nbc-3:nz+nbc+1,nbc+1:nbc+nx);
left=p(nbc+1:nbc+nz,nbc:nbc+4);
right=p(nbc+1:nbc+nz,nx+nbc-3:nx+nbc+1);
end
