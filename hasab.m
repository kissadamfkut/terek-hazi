function [x,y] = hasab(bs,s)

r=10;
a=1;
% Create a unit circle centered at (0,0) using four segments.
switch nargin
    case 0
        x = 8; % four edge segments
        return
    case 1
        A = [0,pi/2,pi,3*pi/2,0,0,0,0; % start parameter values
             pi/2,pi,3*pi/2,2*pi,1,1,1,1; % end parameter values
             1,1,1,1,0,0,0,0; % region label to left
             0,0,0,0,1,1,1,1]; % region label to right
        x = A(:,bs); % return requested columns
        return
    case 2
        
        x=zeros(size(s));
        y=zeros(size(s));
        [m,n]=size(bs);
        if m==1 && n==1
          bs=bs*ones(size(s)); % expand bs
        elseif m~=size(s,1) || n~=size(s,2)
          error(message('pde:scatterg:SizeBs'));
        end
        i=find(bs<=4);
    
        if(length(i))
            x(i) = cos(s(i))*r;
            y(i) = sin(s(i))*r;
        end
        i=find(bs==5);
        if (length(i))
            i=find(bs==5);
            x(i) = ones(length(s(i)),1)'.*-a/2;
            y(i) = interp1([0,1],[a, -a],s(i));
        end
        i=find(bs==6);
        if (length(i))
            x(i) = interp1([0,1],[-a/2, a/2],s(i));
            y(i) = ones(length(s(i)),1)'.*-a;
        end
        i=find(bs==7);
        if (length(i))
            x(i) = ones(length(s(i)),1)'.*a/2;
            y(i) = interp1([0,1],[-a, a],s(i));
        end
        i=find(bs==8);
        if (length(i))
            x(i) = interp1([0,1],[a/2, -a/2],s(i));
            y(i) = ones(length(s(i)),1)'.*a;
        end
        x = x ./10;
        y = y ./10;
        x = x +0.5;
        y = y +0.5;
        
end