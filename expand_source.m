function s=expand_source(s0,nt)
nt0=numel(s0);
if nt0<nt
    s=zeros(nt,1);s(1:nt0)=s0;
else
    s=s0;
end
end