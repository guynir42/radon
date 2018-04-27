function [R, finder] = frt(M_in, varargin)
% fast radon transform (in 2D). Usage: [R, finder] = frt(M_in, varargin)
% mandatory argument: M_in is the input image (can be 3 dimensional)
% OPTIONAL ARGUMENTS:
%   -transpose: if you want to transpose before doing the FRT. Default 0.
%   -partial: save the partial layers of the transformation in a cell array
%    output (useful for getting variance). Default 0. 
%   -expand: add margins to the "passive" axis to include streaks with
%    starting points outside the matrix ("corner cutting streaks"). Default 1.
%   -padding: pad the "active axis" to a power of 2. If this is disabled
%    and the matrix is not a power of 2, the function will fail. Default: 1. 
%   -finder: radon.Finder object that can find short streaks, multiple streaks, etc. 
%   -mex: use the mex function defined in "coreFRT.cpp". This is much
%   faster, so the default is 1. 
%   -threads: to be depricated (it doesn't speed things up). 

    import util.text.cs;
    import util.text.parse_bool;
    import util.text.print_size;
    import radon.frt;

    if nargin==0
        help(mfilename);
        return;
    end

    %%%%%%%%%%%%%%%%%%%%% VARARGIN PAIRS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    use_partials = 0;
    use_mex = 1;
    finder = [];
    transpose = 0;
    use_padding = 1;
    use_expand = 0;
    num_threads = 0;
    
    if ~isempty(varargin) && mod(length(varargin),2)==1
        varargin{end+1} = 1; % positive approach
    end
    
    for ii = 1:2:length(varargin)
       
        key = varargin{ii};
        value = varargin{ii+1};
        
        if cs(key, 'transpose')
            transpose = parse_bool(value);       
        elseif cs(key, 'partial')
            use_partials = parse_bool(value);
        elseif cs(key, 'finder')
            finder = value;
        elseif cs(key, {'use_mex', 'mex'})
            use_mex = parse_bool(value);
        elseif cs(key, {'use_padding', 'padding'})
            use_padding = parse_bool(value);
        elseif cs(key, {'expand', 'use_expand'})
            use_expand = parse_bool(value);
        elseif cs(key, {'threads', 'num_threads'})
            num_threads = value;
         end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%% CHECK INPUTS and DEFAULTS %%%%%%%%%%%%%%%%%%%%%
    
    if ndims(M_in)>2
        error(['Cannot handle 3D matrices... ' printsize(M_in)]);
    end
    
    if ~isempty(finder) 
        
        if ~isa(finder, 'radon.Finder')
            error(['must input a Finder. class(finder)= ' class(finder)]);
        end
        
        finder.last_streak = radon.Streak.empty;
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%% PREPARING THE MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%
    M = M_in;

    if transpose
        M = M';
    end

    if ~isempty(finder) % if we have a finder, we MUST get a variance cell array with partial variances

        finder.im_size_tr = size(M); % the size of the transposed or un-transposed image
        finder.im_size = size(M_in);
        
    end

    if use_padding
        M = padMatrix(M); % make sure the active dim is a power of two (if not, zero pad!).
    end    

    if use_expand
        M = expandMatrix(M); % Increase passive dim by 3, to capture streaks coming from the side... 
    end
    
    %%%%%%%%%%%%%%%%%%%%%%% CORE ALGORITHM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Npassive = size(M,1); % number of pixels along the passive axis (remains constant)
    
    R_partial = {};
    dy = 0; % the different shifts for the current M

    for m = 1:log2(size(M,2))

%         disp(['m= ' num2str(m)]);

        Nactive = size(M,2); % number of pixels along the active axis (shrinks by 2 each iteration)

        M_prev = M; % make a copy of the matrix. That copy will be used to shift and add
        dy_prev = dy; % make a copy of the different shifts, relevant for M_prev 

        max_dy = 2^(m)-1; % the maximum sideways shift is equal to the height of each block (at 45 degrees)    

        dy = -max_dy:max_dy; % a vector of the possible dx values for this level of M
        M = zeros(Npassive, Nactive/2, length(dy), 'like', M_prev); % a new, empty matrix that is shorter on the active axis and bigger on the 3rd axis (dy)
                
        % entry point to mex file...
        if use_mex
            
            radon.coreFRT(M, M_prev, dy, num_threads);
            
        else

            Ndy = length(dy); % size of dy axis
            counter = 1;

            for ii = 1:ceil(Nactive/2) % number of foldings of 2 we can do to the data

                for jj = 1:Ndy % number of possible shifts for this size of slice

                    dy_minor = fix((dy(jj))/2) ; % the dy from M_prev required to fit together the current dy
                    idy_minor = dy_minor + ceil((length(dy_prev))/2); % index of dy_minor on the 3rd axis of M_prev
                    
                    gap_y = dy(jj)-dy_minor;% additional pixel gap we need to shift before we add (NOTE: gap is bigger than +-1 to bridge the difference of start to end points)

                    M1 = M_prev(:,counter,idy_minor);
                    M2 = M_prev(:,counter+1,idy_minor);

                    M(:, ii, jj) = shift_add(M1, M2, -gap_y);
                    
                end % for jj (number of shifts for each subsection)

                counter = counter+2;

            end % for ii (number of foldings in one level)

        end
        
        if ~isempty(finder)
            finder.scan(M, transpose);
        end

        if use_partials
            R_partial{m} = M;
        end

    end % for m (logarithmic jumps)

    if use_partials
        R = R_partial;
    else
        R = permute(M, [1,3,2]); % return just the final Radon transform
    end

    % update the finder with any streaks we have found...
    if ~isempty(finder) 
        finder.finalizeFRT(M_in, transpose, permute(M, [1,3,2]));
    end
            
end

function M_out = padMatrix(M_in, pad_value)
% pads the matrix in the row dimension so it is a power of 2. 
% usage: padMatrix(M_in, pad_value). 
% usually pad_value is equal to 0, 
% but you could imagine using nan instead. 

    if nargin<2 || isempty(pad_value)
        pad_value = 0;
    end

    S = size(M_in);
    
    twos_power = ceil(log2(size(M_in, 2)));
    S(2) = 2.^(twos_power);
    
    if pad_value==0
        M_out = zeros(S, 'like', M_in);
    elseif isnan(pad_value)
        M_out = nan(S, 'like', M_in);
    elseif isinf(pad_value)            
        M_out = inf(S, 'like', M_in);
    elseif isnumeric(pad_value)
        M_out = pad_value*ones(S, 'like', M_in);
    end
    
    M_out(:,1:size(M_in,2)) = M_in;

end

function M_out = expandMatrix(M_in, pad_value, pad_passive_factor)
% Expands a matrix: in coloum dimension, expands it by three 
% (to get streaks that start outside the frame). 
% usuage: expandMatrix(M_in, pad_value=0, pad_passive_factor=1) 
% uses zero padding unless you give another value in pad_value
% if pad_passive_factor>0 it will expand the matrix on the passive
% dimension by this factor on either side! default=1

    if nargin<2 || isempty(pad_value)
        pad_value = 0;
    end
    
    if nargin<3 || isempty(pad_passive_factor)
        pad_passive_factor = 1;
    end

    if pad_passive_factor==0
        M_out = M_in;
        return; 
    end
    
    S(1) = size(M_in,1)+floor(pad_passive_factor*2*size(M_in,2)); % add padding to allow lines at 45 degrees    
    S(2) = size(M_in,2);
    
    if pad_value==0
        M_out = zeros(S, 'like', M_in);
    elseif isnan(pad_value)
        M_out = nan(S, 'like', M_in);
    elseif isinf(pad_value)            
        M_out = inf(S, 'like', M_in);
    elseif isnumeric(pad_value)
        M_out = pad_value*ones(S, 'like', M_in);
    end

    M_out(floor(pad_passive_factor*size(M_in,2))+1:end-ceil(pad_passive_factor*size(M_in,2)), :) = M_in;
    
end

function I_out = shift_add(I1, I2, shift)

    if shift==0 % if shift is zero, do nothing! 
        I_out = I1+I2;
    elseif shift>0 && shift<size(I2,1) % positive shifts... 
        I_out = I1 + [zeros(shift,1,size(I2,3), 'like', I2); I2(1:end-shift,1,:)];
    elseif shift<0 && -shift<size(I2,1) % negative shifts...
        I_out = I1 + [I2(1-shift:end,1,:); zeros(-shift,1,size(I2,3), 'like', I2)];
    else
        I_out = I1;
    end
    
end

% This is no longer used, since we are replacing nan's with zeros in the Finder. 
% Also, this is not triggered if using the mex version (which is a lot faster, and the default). 
function I_out = shift_add_nan(I1, I2, shift)

    if shift==0 % if shift is zero, do nothing! 
        I_out = nansum([I1,I2],2);
    elseif shift>0 && shift<size(I2,1) % positive shifts... 
        I_out = nansum([I1, [zeros(shift,1,size(I2,3), 'like', I2); I2(1:end-shift,1,:)]],2);
    elseif shift<0 && -shift<size(I2,1) % negative shifts...
        I_out = nansum([I1,[I2(1-shift:end,1,:); zeros(-shift,1,size(I2,3), 'like', I2)]],2);
    else
        I_out = I1;
        I_out(isnan(I_out)) = 0;
    end
    
end
