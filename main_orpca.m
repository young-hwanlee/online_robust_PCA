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
% paramters
perc_miss = 0.00;	    % missing entry percentage

sgd = 0;                % use the SGD or not
if sgd == 0             % 2nd order
    r = 20;             % maximum rank of C matrix
    lambda1 = 5e-2;     % constant for the low-rank component X
    lambda2 = 1e-3;     % constant for the sparse component E
elseif sgd == 1         % 1st order
    r = 12;             % maximum rank of C matrix
    lambda1 = 1e-2;     % constant for the low-rank component X
    lambda2 = 1e-6;     % constant for the sparse component E
end

rho = 1e3;              % step size for SGD for A
epsilon = 1e-4;         % tolerance for convergence
n_epoch = 1;            % number of times to repeat Z
iteration = 1;          % when taking an average
outlier_threshold = 1e-3;

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
% 1. Access Video with VideoReader
Z = VideoReader('/Users/young-hwanlee/Desktop/Young-hwan Lee/Project/archive_YounghwanLee/GitHub/orpca/dataset/VIRAT_S_000200_00_000100_000171.mp4');
get(Z)

% 2. Explore Video with IMPLAY
implay('/Users/young-hwanlee/Desktop/Young-hwan Lee/Project/archive_YounghwanLee/GitHub/orpca/dataset/VIRAT_S_000200_00_000100_000171.mp4');

% 3. Develop Your Algorithm
darkCarValue = 50;
darkCar = rgb2gray(read(Z, 500));    % 500 -> frame index
noDarkCar = imextendedmax(darkCar, darkCarValue);
figure(1), imshow(darkCar)
figure(2), imshow(noDarkCar)

% darkCar = rgb2gray(read(Z, 71));
% for i = 1:5
%     darkCarValue = i*10;
%     noDarkCar = imextendedmax(darkCar, darkCarValue);
%     figure(2), subplot(1,5,i), imshow(noDarkCar)
% end

sedisk = strel('disk',2);
noSmallStructures = imopen(noDarkCar, sedisk);
figure(3), imshow(noSmallStructures)

% 4. Apply the Algorithm to the Video
nframes = Z.NumberOfFrames;
I = read(Z, 1);
taggedCars = zeros([size(I,1) size(I,2) 3 nframes], class(I));

% for k = 1 : nframes
for k = 501 : 550
    singleFrame = read(Z, k);

    % Convert to grayscale to do morphological processing.
    I = rgb2gray(singleFrame);

    % Remove dark cars.
    noDarkCars = imextendedmax(I, darkCarValue);

    % Remove lane markings and other non-disk shaped structures.
    noSmallStructures = imopen(noDarkCars, sedisk);

    % Remove small structures.
    noSmallStructures = bwareaopen(noSmallStructures, 150);

    % Get the area and centroid of each remaining object in the frame. The
    % object with the largest area is the light-colored car.  Create a copy
    % of the original frame and tag the car by changing the centroid pixel
    % value to red.
    taggedCars(:,:,:,k) = singleFrame;

    stats = regionprops(noSmallStructures, {'Centroid','Area'});
    if ~isempty([stats.Area])
        areaArray = [stats.Area];
        [junk,idx] = max(areaArray);
        c = stats(idx).Centroid;
        c = floor(fliplr(c));
        width = 2;
        row = c(1)-width:c(1)+width;
        col = c(2)-width:c(2)+width;
        taggedCars(row,col,1,k) = 255;
        taggedCars(row,col,2,k) = 0;
        taggedCars(row,col,3,k) = 0;
    end
end

% 5. Visualize Results
frameRate = Z.FrameRate;
implay(taggedCars,frameRate);

%%
starting_time = clock;

for i = 1:iteration
    rng(i);

    L_init = randn(N,r);
    
    % generate missing entry indexes
    input(i).L_init = L_init;
    input(i).Omega = ones(N,iter);
    input(i).Omega(rand(N*iter,1) < perc_miss) = 0;
    
%     [output(i).L,output(i).R,output(i).E,output(i).Delta_Lt_traj,output(i).fx_traj,output(i).fx_accumulate] = ...
%         orpca_with_missing_v2(Z.*input(i).Omega,N,iter,L_init,input(i).Omega,lambda1,lambda2,epsilon,r,rho,sgd);
    [output(i).L,output(i).R,output(i).E,output(i).Delta_Lt_traj,output(i).fx_traj,output(i).fx_accumulate] = ...
        orpca_with_missing(Z,N,iter,L_init,input(i).Omega,lambda1,lambda2,epsilon,r,rho,sgd);
    
    X = output(i).L*output(i).R';

    % compute performance (error) metrics excluding the transient portion
    Z_end = Z(:,t_transient+1:end);
    X_end = X(:,t_transient+1:end);
    E_end = output(i).E(:,t_transient+1:end);
    Omega_end = input(i).Omega(:,t_transient+1:end);
    outlier_end = outlier(:,t_transient+1:end);
    no_outlier_end = no_outlier(:,t_transient+1:end);
    
    rank_Z = rank(Z_end)
    rank_X = rank(X_end)
    
    nmse(i) = norm(Z_end-X_end,'fro')^2 / norm(Z_end,'fro')^2
    nmse_miss(i) = norm((Z_end-X_end).*(1-Omega_end),'fro')^2 / norm(Z_end.*(1-Omega_end),'fro')^2
    nmse_miss_no_outlier(i) = norm((Z_end-X_end).*(1-Omega_end).*no_outlier_end,'fro')^2 / norm(Z_end.*(1-Omega_end).*no_outlier_end,'fro')^2
    norm_E(i) = norm(E_end.*outlier_end, 'fro')

    true_nmse(i) = norm(Z_end.*outlier_end,'fro')^2 / norm(Z_end,'fro')^2;
    true_norm_E(i) = norm(Z_end.*outlier_end, 'fro');
    true_DC_nmse(i) = norm(Z_end.*no_outlier_end - X_end,'fro')^2 / norm(Z_end.*no_outlier_end,'fro')^2

    detected_inlier = zeros(size(E_end));
    detected_inlier(E_end <= outlier_threshold) = 1;
    detected_outlier = 1 - detected_inlier;

    false_alarm(i) = length(nonzeros(detected_outlier.*no_outlier_end)) / length(nonzeros(no_outlier_end))
    missed_detection(i) = length(nonzeros(detected_inlier.*outlier_end)) / length(nonzeros(outlier_end))

    % plots
    figure(2);
    subplot(2,1,2);
    plot(t, X);
    xlim([(n_epoch-1) n_epoch]*T);
    grid on;
    xlabel('Number of samples');
    ylabel('X');

    subplot(2,1,1);
    plot(1:T, Z_end.*no_outlier_end);
    xlim([0 T]);
    ylim([-0.05 0.05]);
    xlabel('Number of samples');
    ylabel('True X');
    grid on;
    
    figure(3);
    subplot(2,1,2);
    plot(1:iter, output(i).E);
    xlim([(n_epoch-1) n_epoch]*T);
    ylim([-0.3 0.1]);
    grid on;
    xlabel('Number of samples');
    ylabel('E');
    
    subplot(2,1,1);
    plot(1:T, Z_end.*outlier_end);
    xlim([0 T]);
    ylim([-0.3 0.1]);
    xlabel('Number of samples');
    ylabel('True E');
    grid on;

    % check convergence
    figure(98);
    subplot(3,1,1);
    plot(1:iter,output(i).Delta_Lt_traj);
    xlim([0 iter]);
    title('||L_t-L_{t-1}||_F');
    xlabel('Iteration');
    grid on;

    subplot(3,1,2);
    plot(1:iter,output(i).fx_traj);
    xlim([0 iter]);
    title('(1/2)*||P_{\Omega_t}(z_t)-L_tr_t-e_t||_2^2 + (\lambda_1/2)*(||L_t||_F^2+||r_t||_2^2) + \lambda_2*||P_{\Omega_t}(e_t)||_1');
    xlabel('Iteration');
    grid on;

    for jj = 1:length(output(i).fx_traj)
%         fx_conv(jj) = sum(output(i).fx_traj(1:jj))/jj;
        fx_conv(jj) = output(i).fx_accumulate(jj)/jj;          % faster
    end
    subplot(3,1,3);
    plot(1:iter,fx_conv);
    xlim([0 iter]);
    title('1/T\Sigma(1/2)*||P_{\Omega_t}(z_t)-L_tr_t-e_t||_2^2 + (\lambda_1/2)*(||L_t||_F^2+||r_t||_2^2) + \lambda_2*||P_{\Omega_t}(e_t)||_1');
    xlabel('Iteration');
    grid on;
    
%     % % save files
%     save('input_orpca_lambda_5e_1_mu_5e_3_miss_5perc_20200125.mat', 'input', '-v7.3');
%     save('output_orpca_lambda_5e_1_mu_5e_3_miss_5perc_20200125.mat', 'output', '-v7.3');
end

avg_nmse = mean(nmse, 1);
avg_nmse_miss = mean(nmse_miss, 1);
avg_nmse_miss_no_outlier = mean(nmse_miss_no_outlier, 1);

%%
% % plots
% % true measurement
figure(1);
plot(1:iter, Z);
xlim([(n_epoch-1) n_epoch]*T);
ylim([-0.3 0.1]);
grid on;
xlabel('Number of samples');
ylabel('Z');

ending_time = clock;
elapsed_time = etime(ending_time, starting_time)/60

% load handel
% sound(y,Fs)


