clear all;

%%
% % 1. Access Video with VideoReader
% trafficVid = VideoReader('traffic.mj2')
% get(trafficVid)
% 
% % 2. Explore Video with IMPLAY
% implay('/Users/young-hwanlee/Desktop/Young-hwan Lee/Project/archive_YounghwanLee/GitHub/orpca/dataset/VIRAT_S_000200_00_000100_000171.mp4');
% 
% % 3. Develop Your Algorithm
% darkCarValue = 50;
% darkCar = rgb2gray(read(trafficVid, 71));
% noDarkCar = imextendedmax(darkCar, darkCarValue);
% imshow(darkCar)
% figure, imshow(noDarkCar)
% 
% sedisk = strel('disk',2);
% noSmallStructures = imopen(noDarkCar, sedisk);
% imshow(noSmallStructures)
% 
% % 4. Apply the Algorithm to the Video
% nframes = trafficVid.NumberOfFrames;
% I = read(trafficVid, 1);
% taggedCars = zeros([size(I,1) size(I,2) 3 nframes], class(I));
% 
% for k = 1 : nframes
%     singleFrame = read(trafficVid, k);
% 
%     % Convert to grayscale to do morphological processing.
%     I = rgb2gray(singleFrame);
% 
%     % Remove dark cars.
%     noDarkCars = imextendedmax(I, darkCarValue);
% 
%     % Remove lane markings and other non-disk shaped structures.
%     noSmallStructures = imopen(noDarkCars, sedisk);
% 
%     % Remove small structures.
%     noSmallStructures = bwareaopen(noSmallStructures, 150);
% 
%     % Get the area and centroid of each remaining object in the frame. The
%     % object with the largest area is the light-colored car.  Create a copy
%     % of the original frame and tag the car by changing the centroid pixel
%     % value to red.
%     taggedCars(:,:,:,k) = singleFrame;
% 
%     stats = regionprops(noSmallStructures, {'Centroid','Area'});
%     if ~isempty([stats.Area])
%         areaArray = [stats.Area];
%         [junk,idx] = max(areaArray);
%         c = stats(idx).Centroid;
%         c = floor(fliplr(c));
%         width = 2;
%         row = c(1)-width:c(1)+width;
%         col = c(2)-width:c(2)+width;
%         taggedCars(row,col,1,k) = 255;
%         taggedCars(row,col,2,k) = 0;
%         taggedCars(row,col,3,k) = 0;
%     end
% end
% 
% % 5. Visualize Results
% frameRate = trafficVid.FrameRate;
% implay(taggedCars,frameRate);

%%
% v = VideoReader('VIRAT_S_000200_00_000100_000171.mp4');
% v.CurrentTime = 15;     % 15 second
% currAxes = axes;
% 
% while hasFrame(v)
%     vidFrame = readFrame(v);
%     image(vidFrame, 'Parent', currAxes);
%     currAxes.Visible = 'off';
%     pause(1/v.FrameRate);
% end
% 
% while v.CurrentTime <= 20
%     vidFrame = readFrame(v);
%     image(vidFrame, 'Parent', currAxes);
%     currAxes.Visible = 'off';
% %     pause(1/v.FrameRate);
% end
% 
% implay('/Users/young-hwanlee/Desktop/Young-hwan Lee/Project/archive_YounghwanLee/GitHub/orpca/dataset/VIRAT_S_000200_00_000100_000171.mp4');

%%
% Z = VideoReader('/Users/young-hwanlee/Desktop/Young-hwan Lee/Project/archive_YounghwanLee/GitHub/orpca/dataset/VIRAT_S_000200_00_000100_000171.mp4');
% 
% % Z.CurrentTime = 15;
% % currAxes = axes;
% % while Z.CurrentTime <= 20
% %     vidFrame = readFrame(Z);
% %     image(vidFrame, 'Parent', currAxes);
% % %     currAxes.Visible = 'off';
% % end
% 
% N = Z.Height;
% T = Z.Width;
% 
% iter = n_epoch*T;

%%
% paramters
perc_miss = 0.00;	% missing entry percentage

sgd = 0;            % don't use the SGD (2nd order)
% sgd = 1;            % use the SGD (1st order)

lambda1 = 1e1;     % constant for the low-rank component X
lambda2 = 1e-3;     % constant for the sparse component E
rho = 1e2;          % step size for SGD for L
epsilon = 1e-3;     % tolerance for convergence

n_epoch = 1;        % number of times to repeat Z
simul_count = 1;    % simulation count to take an average

%%
% access Video with VideoReader
Z_tmp = VideoReader('/Users/young-hwanlee/Desktop/Young-hwan Lee/Project/archive_YounghwanLee/GitHub/orpca/dataset/VIRAT_S_000200_00_000100_000171.mp4');

% % vectorize every frame to form the input matrix Z
% videof = read(Z_tmp);
% for i = 1:649
%     frg = rgb2gray(videof(:,:,:,i));
%     videog(:,i) = frg(:);
% end

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

Z = zeros(nframes, N_tmp*T_tmp);
% Z = zeros(N_tmp*T_tmp, nframes);
N = size(Z,1);
T = size(Z,2);

start_frame = 100;
for k = start_frame+1 : start_frame+nframes
    singleFrame = rgb2gray(read(Z_tmp, k));
%     Z(k-start_frame,:) = reshape(singleFrame,[],1);
%     Z(:,k-start_frame) = reshape(singleFrame,[],1);
    Z(k-start_frame,:) = reshape(singleFrame(N_tmp_start+1:end, T_tmp_start+1:T_tmp_end),[],1);
    figure(1); imshow(singleFrame(N_tmp_start+1:end,T_tmp_start+1:T_tmp_end));
end

%%
% apply the algorithm to the Video
iter = n_epoch*T;
t = linspace(0,iter,iter);
Z = repmat(Z,[1,n_epoch]);
t_transient = (n_epoch-1)*T;    % time not included in performance metric

starting_time = clock;

for i = 1:simul_count
    rng(i);

%     L_init = randn(N,r);
    L_init = randn(N,nframes);
    
    % generate missing entry indices
    input(i).L_init = L_init;
    input(i).Omega = ones(N,iter);
    input(i).Omega(rand(N*iter,1) < perc_miss) = 0;
    
%     [output(i).L,output(i).R,output(i).E,output(i).Delta_Lt_traj,output(i).fx_traj,output(i).fx_accumulate] = ...
%         orpca_with_missing(Z,N,iter,L_init,input(i).Omega,lambda1,lambda2,epsilon,r,rho,sgd);
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
for i = 1:nframes
    Z_video(:,:,i) = reshape(Z_end(i,:),[N_tmp,T_tmp]);
    X_video(:,:,i) = reshape(X_end(i,:),[N_tmp,T_tmp]);
    E_video(:,:,i) = reshape(E_end(i,:),[N_tmp,T_tmp]);
end

% Z = reshape(Z_end,[N_tmp,T_tmp,nframes]);
% X = reshape(X_end,[N_tmp,T_tmp,nframes]);
% E = reshape(E_end,[N_tmp,T_tmp,nframes]);

original_data_Z = VideoWriter('original_data_Z.avi');
low_rank_X = VideoWriter('low_rank_X.avi');
sparse_E = VideoWriter('sparse_E.avi');

open(original_data_Z);
for i = 1:nframes
    frame = Z_video(:,:,i);
    writeVideo(original_data_Z,mat2gray(frame));
end
close(original_data_Z);

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

%%
ending_time = clock;
elapsed_time = etime(ending_time, starting_time)/60

load handel
sound(y,Fs)


