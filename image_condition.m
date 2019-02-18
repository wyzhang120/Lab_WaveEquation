function img=image_condition(img,p1,q0,nz,nx,nbc)
p11=p1(nbc:nbc+nz+1,nbc:nbc+nx+1);
q00=q0(nbc:nbc+nz+1,nbc:nbc+nx+1);
img_temp=(circshift(p11,[-1,0,0])+circshift(p11,[1,0,0])+circshift(p11,[0,-1,0])...
    +circshift(p11,[0,1,0])-4*p11).*q00;
img=img+img_temp(2:end-1,2:end-1);
end