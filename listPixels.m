function [xlists_cell, ylists_cell, num_pixels] = listPixels(x1,x2,y1,y2, im_size)
% usage: [xlists_cell, ylists_cell] = listPixels(x1,x2,y1,y2, im_size)
% makes a list of all the pixels of a line going from [x1,y1] to [x2,y2] in 
% an image of size "im_size". If multiple values are given for the x's and
% and y's then multiple lists are returned. 
% 
% INPUTS: position of start and end points of the streak, in units of pixels.
% If the four x/y inputs are vectors (of the same size) then multiple lines 
% are returned. The input x/y may be unrounded. 
% Also input the size of the image the line belongs to. 
%
% OUTPUTS: two cell arrays, one for the list of x values and one for the 
% list of y values where each streak passes. 
% Also returns a vector with the number of pixels in each line. 

    im_size = util.vec.imsize(im_size, 'size');
    
    N = max([length(x1), length(x2), length(y1), length(y2)]);

    if length(x1)<N, x1 = [x1 ones(1, N-length(x1))*x1(end)]; end
    if length(x2)<N, x2 = [x2 ones(1, N-length(x2))*x2(end)]; end
    if length(y1)<N, y1 = [y1 ones(1, N-length(y1))*y1(end)]; end
    if length(y2)<N, y2 = [y2 ones(1, N-length(y2))*y2(end)]; end

    xlists_cell{N} = [];
    ylists_cell{N} = [];
    num_pixels = zeros(1,N);
    
    for ii = 1:N

        % verify line is inside the frame at all!
        if x1(ii)<1 && x2(ii)<1, continue; end
        if x1(ii)>im_size(2) && x2(ii)>im_size(2), continue; end
        if y1(ii)<1 && y2(ii)<1, continue; end
        if y1(ii)>im_size(1) && y2(ii)>im_size(1), continue; end
        
        pix_x1 = x1(ii);
        pix_x2 = x2(ii);
        pix_y1 = y1(ii);
        pix_y2 = y2(ii);

        a = (pix_y2-pix_y1)./(pix_x2-pix_x1);
        
        if isnan(a) || abs(a)>=1 % vertical (or closer to vertical) lines
            
            if pix_y1<pix_y2
                y = pix_y1:pix_y2; % not rounded!
            else
                y = pix_y2:pix_y1; % not rounded!
            end

            if isnan(a) % this happens if pix_x1==pix_x2
                x = ones(1,length(y)).*pix_x1;
            else
                if pix_y1<pix_y2
                    x = pix_x1+(y-pix_y1)./a; % not rounded!
                else
                    x = pix_x2+(y-pix_y2)./a; % not rounded!
                end
            end
            
            ind_low = find(y<0.5, 1, 'last');
            if isempty(ind_low), ind_low=0; end

            ind_high = find(y>im_size(1)+0.5, 1, 'first');
            if isempty(ind_high), ind_high=length(y)+1; end

            y = y(ind_low+1:ind_high-1);
            x = x(ind_low+1:ind_high-1);

            if x(1)<x(end)

                ind_low = find(x<0.5, 1, 'last');
                if isempty(ind_low), ind_low=0; end

                ind_high = find(x>im_size(2)+0.5, 1, 'first');
                if isempty(ind_high), ind_high=length(x)+1; end

            else

                ind_low = find(x<im_size(2)+0.5, 1, 'first');
                if isempty(ind_low), ind_low=0; end

                ind_high = find(x>0.5, 1, 'last');
                if isempty(ind_high), ind_high=length(x)+1; end

            end

            y = y(ind_low+1:ind_high-1);
            x = x(ind_low+1:ind_high-1);

            xlists_cell{ii} = round(x);
            ylists_cell{ii} = round(y);

        elseif abs(a)<1

            if pix_x1<pix_x2
                x = pix_x1:pix_x2; % not rounded!
                y = pix_y1+(x-pix_x1).*a; % not rounded!
            else
                x = pix_x2:pix_x1; % not rounded!
                y = pix_y2+(x-pix_x2).*a; % not rounded!
            end

            ind_low = find(x<0.5, 1, 'last');
            if isempty(ind_low), ind_low=0; end

            ind_high = find(x>im_size(2)+0.5, 1, 'first');
            if isempty(ind_high), ind_high=length(x)+1; end

            x = x(ind_low+1:ind_high-1);
            y = y(ind_low+1:ind_high-1);

            if y(1)<y(end)

                ind_low = find(y<0.5, 1, 'last');
                if isempty(ind_low), ind_low=0; end

                ind_high = find(y>im_size(1)+0.5, 1, 'first');
                if isempty(ind_high), ind_high=length(y)+1; end

            else

                ind_low = find(y<im_size(1)+0.5, 1, 'first');
                if isempty(ind_low), ind_low=0; end

                ind_high = find(y>0.5, 1, 'last');
                if isempty(ind_high), ind_high=length(y)+1; end

            end

            x = x(ind_low+1:ind_high-1);
            y = y(ind_low+1:ind_high-1);

            xlists_cell{ii} = round(x);
            ylists_cell{ii} = round(y);
            
        end

        num_pixels(ii) = numel(xlists_cell{ii});

    end % for ii
    
    
end