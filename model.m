function M = model(im_size, x1, x2, y1, y2, sigma, replace_value, threshold, oversample)
% usage: M = model(im_size, x1, x2, y1, y2, sigma, replace_value, threshold, oversample)
% Makes a streak model to compare to data when fitting. 

    if nargin==0, help('radon.model'); return; end
        
    if nargin<7 || isempty(replace_value)
        replace_value = 0;
    end
    
    if nargin<8 || isempty(threshold)
        threshold = 1e-10;
    end
    
    if nargin<9 || isempty(oversample)
        oversample = 4;
    end
    
    if isscalar(im_size)
        im_size = im_size.*[1 1];
    else
        im_size = im_size(1:2);
    end
    
    if oversample
        im_size = im_size.*oversample;
        x1 = (x1-0.5).*oversample + 0.5;
        x2 = (x2-0.5).*oversample + 0.5;
        y1 = (y1-0.5).*oversample + 0.5;
        y2 = (y2-0.5).*oversample + 0.5;
        sigma = sigma.*oversample;
    end
    
    [x,y] = meshgrid(0:im_size(2)-1, 0:im_size(1)-1);
    
    a = (y2-y1)./(x2-x1); % slope parameter
    b = (y1.*x2 - y2.*x1)./(x2-x1); % impact parameter
    
    if x1==x2
        d = abs(x-x1); % distance from vertical line
    else
        d = abs(a.*x - y + b)./sqrt(1+a.^2); % distance from line
    end
            
    M0 = (1./sqrt(2.*pi)./sigma).^1.*exp(-0.5.*d.^2./sigma.^2); % infinit line
    
    % must clip this line! 
    if x1==x2 && y1==y2 % this is extremely unlikely to happen...
        M0 = zeros(size(M0));
    elseif x1==x2 % vertical line (a is infinite)
        if y1>y2
            M0(y>y1) = 0;
            M0(y<y2) = 0;
        else
            M0(y<y1) = 0;
            M0(y>y2) = 0;
        end
    elseif y1==y2 % horizontal line
        if x1>x2
            M0(x>x1) = 0;
            M0(x<x2) = 0;
        else
            M0(x<x1) = 0;
            M0(x>x2) = 0;
        end
    elseif y1<y2
        M0(y<(-1./a.*x+y1+1./a.*x1)) = 0;
        M0(y>(-1./a.*x+y2+1./a.*x2)) = 0;
    else
        M0(y>(-1./a.*x+y1+1./a.*x1)) = 0;
        M0(y<(-1./a.*x+y2+1./a.*x2)) = 0;
    end
    
    M1 = (1./sqrt(2.*pi)./sigma).^1.*exp(-0.5.*((x-x1).^2+(y-y1).^2)./sigma.^2);
    M2 = (1./sqrt(2.*pi)./sigma).^1.*exp(-0.5.*((x-x2).^2+(y-y2).^2)./sigma.^2);
    
    M = max(M0,max(M1,M2));
    
    if oversample
        M = util.img.downsample(M, oversample)./oversample;
    end
    
    M(M<threshold) = replace_value;
    
%     M = M./util.stat.sum2(M);
    
end