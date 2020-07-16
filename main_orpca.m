clear all;

%%
% paramters
perc_miss = 0.00;	% missing entry percentage

sgd = 0;            % don't use the SGD (2nd order)
% sgd = 1;            % use the SGD (1st order)

lambda1 = 1e1;      % constant for the low-rank component X
lambda2 = 1e-3;     % constant for the sparse component E
rho = 1e2;          % step size for SGD for L
epsilon = 1e-3;     % tolerance for convergence

n_epoch = 1;        % number of times to repeat Z
simul_count = 1;    % simulation count to take an average

%%
% access Video with VideoReader
Z_tmp = VideoReader('path to the video file');

% decrease the frame size and frame numbers in order to save the running time
% nframes = Z_tmp.NumberOfFrames;
nframes = 100;
frame_rate = Z_tmp.FrameRate;

% N_tmp = Z_tmp.Height;
% T_tmp = Z_tmp.Width;
N_tmp_start = 520;
T_tmp_start = 300;
T_tmp_end = 500;
N_tmp = Z_tmp.Height - N_tmp_start;
T_tmp = T_tmp_end - T_tmp_start;

% vectorize every frame to form the input matrix Z
Z = zeros(nframes, N_tmp*T_tmp);
N = size(Z,1);
T = size(Z,2);

start_frame = 100;
for k = start_frame+1 : start_frame+nframes
    singleFrame = rgb2gray(read(Z_tmp, k));
    Z(k-start_frame,:) = reshape(singleFrame(N_tmp_start+1:end, T_tmp_start+1:T_tmp_end),[],1);
    figure(1); imshow(singleFrame(N_tmp_start+1:end,T_tmp_start+1:T_tmp_end));
end

%%
% apply the algorithm to the video file
iter = n_epoch*T;
t = linspace(0,iter,iter);
Z = repmat(Z,[1,n_epoch]);
t_transient = (n_epoch-1)*T;    % time not included in performance metric

starting_time = clock;

for i = 1:simul_count
    rng(i);

    L_init = randn(N,nframes);
    
    % generate missing entry indices
    input(i).L_init = L_init;
    input(i).Omega = ones(N,iter);
    input(i).Omega(rand(N*iter,1) < perc_miss) = 0;

    [output(i).L,output(i).R,output(i).E] = ...
        orpca_with_missing(Z,N,iter,L_init,input(i).Omega,lambda1,lambda2,epsilon,nframes,rho,sgd);
    
    X = output(i).L*output(i).R';

    % compute performance (error) metrics excluding the transient portion
    Z_end = Z(:,t_transient+1:end);
    X_end = X(:,t_transient+1:end);
    E_end = output(i).E(:,t_transient+1:end);
    Omega_end = input(i).Omega(:,t_transient+1:end);
    
    rank_Z = rank(Z_end)
    rank_X = rank(X_end)
    
    nmse(i) = norm(Z_end-X_end,'fro')^2 / norm(Z_end,'fro')^2;
    norm_E(i) = norm(E_end, 'fro')
end

avg_nmse = mean(nmse, 1);

%%
% save the low rank matrix X and sparse matrix X as video files
for i = 1:nframes
    Z_video(:,:,i) = reshape(Z_end(i,:),[N_tmp,T_tmp]);
    X_video(:,:,i) = reshape(X_end(i,:),[N_tmp,T_tmp]);
    E_video(:,:,i) = reshape(E_end(i,:),[N_tmp,T_tmp]);
end

data_Z = VideoWriter('data_Z.avi');
low_rank_X = VideoWriter('low_rank_X.avi');
sparse_E = VideoWriter('sparse_E.avi');

open(data_Z);
for i = 1:nframes
    frame = Z_video(:,:,i);
    writeVideo(data_Z,mat2gray(frame));
end
close(data_Z);

open(low_rank_X);
for i = 1:nframes
    frame = X_video(:,:,i);
    writeVideo(low_rank_X,mat2gray(frame));
end
close(low_rank_X);

open(sparse_E);
for i = 1:nframes
    frame = E_video(:,:,i);
    writeVideo(sparse_E,mat2gray(frame));
end
close(sparse_E);

ending_time = clock;
elapsed_time = etime(ending_time, starting_time)/60

% to make a sound after the simulation is finished
load handel
sound(y,Fs)


