function v=padvel(v0,nbc)
v=[repmat(v0(:,1),1,nbc), v0, repmat(v0(:,end),1,nbc)];
v=[repmat(v(1,:),nbc,1); v; repmat(v(end,:),nbc,1)];
end