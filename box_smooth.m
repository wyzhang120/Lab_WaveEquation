function s=box_smooth(s,h1,h2,h3,h4)
check_vel(s,nargin);
s0=s;s(:)=0;
nd=ndims(s);
if nd==2
    if nargin == 3
        nit=1;
    else
        nit=h3;
    end
    hh1=(h1-1)/2; hh2=(h2-1)/2;
    [n1,n2]=size(s0);
    for it=1:nit
        %display(['Iteration ',num2str(it)]);
        % s0=box_padding(s0,h1,h2);
        s0=padarray(s0,[hh1,hh2],'replicate');
        % s0=mirror_padding(s0,h1,h2);
        for i2=hh2+1:hh2+n2
            for i1=hh1+1:hh1+n1
                i11=i1-hh1; i12=i1+hh1;
                i21=i2-hh2; i22=i2+hh2;
                ss=s0(i11:i12,i21:i22);
                s(i1,i2)=sum(sum(ss))/(h1*h2);
            end
        end
        s0=s(hh1+1:hh1+n1,hh2+1:hh2+n2);
    end
    s=s0;
elseif nd==3
    if nargin==4
        nit=1;
    else
        nit=h4;
    end
    hh1=(h1-1)/2; hh2=(h2-1)/2;hh3=(h3-1)/2;
 
    [n1,n2,n3]=size(s0);
    csh=zeros(h1*h2*h3,3); ii=0;
    for i3=-hh3:hh3
        for i2=-hh2:hh2
            for i1=-hh1:hh1
           ii=ii+1;     
        csh(ii,1)=i1;csh(ii,2)=i2;csh(ii,3)=i3;
            end
        end
    end
    
    for it=1:nit
        display(['Iteration ',num2str(it)]);
        s0=padarray(s0,[hh1,hh2,hh3],'replicate');
        s=zeros(size(s0));
        for i=1:ii
            s=s+circshift(s0,csh(i,:));
        end
        s=s/(h1*h2*h3);
        s0=s(hh1+1:hh1+n1,hh2+1:hh2+n2,hh3+1:hh3+n3);
    end
          
    s=s0;
end
end

function check_vel(vel,n)
nn=ndims(vel);
if nn==2
    if ismatrix(vel)
        
        if ~(n==3 || n==4)
            error('Check the input argument number!');
        end
    else
        error('Check the dimensions of input velocity!');
    end
elseif nn==3
    if ~(n==4 || n==5)
        error('Check the input argument number!');
    end
else
    error('Check the dimensions of input velocity!');
end
end

