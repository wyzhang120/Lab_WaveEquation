function [vel_smooth,refl]=vel_smooth(vel,h1,h2,h3,h4)
check_vel(vel,nargin);

s=1.0./vel;
if ndims(vel)<=2
    if nargin==3
        %     ss=cross_smooth(s,h1,h2,1);
        ss=box_smooth(s,h1,h2,1);
    elseif nargin==4
        %     ss=cross_smooth(s,h1,h2,h3);
        ss=box_smooth(s,h1,h2,h3);
    end
elseif ndims(vel)==3
    if nargin==4
        ss=box_smooth(s,h1,h2,h3);
    elseif nargin==5
        ss=box_smooth(s,h1,h2,h3,h4);
    end
end

vel_smooth=1.0./ss;

% reflectivity-like term for migration is defined as
% 2*s_0*delta_s

%refl=2*ss.*(s-ss);
refl=(s-ss)./s;

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