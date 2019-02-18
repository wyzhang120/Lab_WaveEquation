function p=load_boundary(p,bc_top,bc_bottom,bc_left,bc_right,nz,nx,nbc)
p(nbc:nbc+4,nbc+1:nbc+nx)=bc_top;
p(nz+nbc-3:nz+nbc+1,nbc+1:nbc+nx)=bc_bottom;
p(nbc+1:nbc+nz,nbc:nbc+4)=bc_left;
p(nbc+1:nbc+nz,nx+nbc-3:nx+nbc+1)=bc_right;
end