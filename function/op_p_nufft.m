function [A, At, G, W, Gw] = op_p_nufft(pppp, N, Nn, No, Ns, ww, param)

% Create the nonuniform gridding matrix and fft operators to be used for
% parallel processing
%
% in:
% p{:}[2] - nonuniformly distributed frequency location points for each
%           cell member which will be treated in parallel
% N[2]    - size of the reconstruction image
% Nn[2]   - size of the kernels (number of neighbors considered on each direction)
% No[2]   - oversampled fft from which to recover the non uniform fft via
%           kernel convolution
% Ns[2]   - fft shift
%
% out:
% A[@]          - function handle for direct operator
% At[@]         - function handle for adjoint operator
% G{:}[:][:]    - convolution kernel matrix (small) associated with each
%               patch in the fourier plane
% W{:}          - mask of the values that contribute to the convolution
% Gw[:][:]      - global convolution kernel matrix

if ~exist('param', 'var')
    param = struct();
end
if ~isfield(param, 'use_nufft_blocks'), param.use_nufft_blocks = 1; end
if ~isfield(param, 'gen_only_fft_op'), param.gen_only_fft_op = 0; end
if ~exist('ww', 'var')
    ww = cell(length(pppp), 1);
    for q=1:length(pppp)
        ww{q} = ones(length(pppp{q}(:, 1)), 1);
    end
end

R = size(pppp, 1);
% R = size222(pppp, 1);%yx
if param.gen_only_fft_op
    [A, At, ~, ~] = op_nufft([0, 0], N, Nn, No, Ns);
    G = [];
    W = [];
    Gw = [];
else
    if ~param.use_nufft_blocks
        %% compute the overall gridding matrix and its associated kernels
        [A, At, Gw, ~] = op_nufft(cell2mat(pppp), N, Nn, No, Ns);

        %% compute small gridding matrices associated with each parallel block
        G = cell(R, 1);
        W = cell(R, 1);

        % block start position
        fprintf('\nComputing block matrices ...\n');
        b_st = 1;
        for q = 1:R
            tstart = tic;
            % current block length
            % the matrix Gw is structured identical to the structure of p thus we 
            % grab it block by block
            b_l = length(pppp{q});

            % get a block out of the large G and trim it
            Gw(b_st:b_st+b_l-1, :) = spdiags(ww{q}, 0, b_l, b_l) * Gw(b_st:b_st+b_l-1, :);
            Gb = Gw(b_st:b_st+b_l-1, :);

            %% now trim the zero rows and store a mask in W

            % preallocate W for speed
            W{q} = false(No(1)*No(2), 1);

            % use the absolute values to speed up the search
            Gb_a = abs(Gb);

            % check if eack line is entirely zero
            W{q} = Gb_a' * ones(b_l, 1) ~= 0;
            
            % store only what we need from G
            G{q} = Gb(:, W{q});

            % iterate among the blocks
            b_st = b_st+b_l;
            tend = toc(tstart);
            fprintf('Block matrix %d: %ds \n', q, ceil(tend));
        end
    else

        %% compute small gridding matrices associated with each parallel block
        
        Gw = spalloc(length(cell2mat(pppp)), No(1)*No(2), 16 * length(cell2mat(pppp)));
        G = cell(R, 1);
        W = cell(R, 1);

        b_st = 1;
        % block start position
        fprintf('\nComputing block matrices ...\n');
        for q = 1:R

            tstart = tic;
            b_l = length(pppp{q});


            %% compute the small gridding matrix and its associated kernels
            [~, ~, Gb, ~] = op_nufft([pppp{q, 1} pppp{q, 2}], N, Nn, No, Ns);

            %% now trim the zero rows and store a mask in W

            % preallocate W for speed
            W{q} = false(No(1)*No(2), 1);

            Gb = spdiags(ww{q}, 0, b_l, b_l) * Gb;
            
            % use the absolute values to speed up the search
            Gb_a = abs(Gb);

            % check if eack line is entirely zero
            W{q} = Gb_a' * ones(size(Gb, 1), 1) ~= 0;

            % store only what we need from G
            G{q} = Gb(:, W{q});

            %% fill the whole Gw
            Gw(b_st:b_st+b_l-1, :) = Gb;

            b_st = b_st+b_l;
            tend = toc(tstart);
            fprintf('Block matrix %d: %ds \n', q, ceil(tend));
        end

        [A, At, ~, ~] = op_nufft([0, 0], N, Nn, No, Ns);
    end
end


end

