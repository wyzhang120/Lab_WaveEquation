function q=pertubation(p1,q,refl,beta_dt)
q=q+(circshift(p1,[1,0,0])+circshift(p1,[-1,0,0])+circshift(p1,[0,1,0])+circshift(p1,[0,-1,0])...
    -4.0*p1).*refl.*beta_dt;
end