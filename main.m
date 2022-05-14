close all
clear
clc

#import 2d geometry utils
source "tools/utilities/geometry_helpers_2d.m"
source "functions/odometry.m"
source "functions/get_initial_guess.m"
source "functions/initial_guess_matching.m"
source "functions/get_obs_ass.m"
source "functions/build_obs.m"
addpath("functions")
addpath("tools/g2o_wrapper")
addpath("Datasets")
addpath("LS")
source "LS/Optimization_2d.m"
source "LS/Optimization_2d_odometry.m"
source "LS/Optimization_2d_edges.m"
source "tools/visualization/plotter.m"
source "tools/utilities/preprocessTransitions.m"

# h = figure(1); #attenzione

%%------------------------------------------%%
%%--------------LOAD DATASETS---------------%%
%%------------------------------------------%%

% The dataset is composed by a g2o file which contain poses and range observations. 
% The file contain also odometry edges that are used to construct the initial guess
% for the problem.

% LOADING BOTH DATASETS
[~, poses, transitions, observations] = loadG2o('slam2d_range_only_initial_guess.g2o');
[landmarks_ground_truth, poses_ground_truth, transitions_ground_truth, observations_ground_truth] = loadG2o('slam2d_range_only_ground_truth.g2o');

%% ------------------------------------ %%
%% ----------- INITIAL GUESS ---------- %%
%% ------------------------------------ %%
printf("\n\n%% ------------------------------------ %%\n%% ----- Generating an initial guess -- %%\n%% ------------------------------------ %%\n\n");
fflush(stdout);

% here i'm getting a LIST of HOW MANY landarmsks and their ID
available_landmarks = zeros(length(landmarks_ground_truth),1);
for l=1:length(landmarks_ground_truth)
	available_landmarks(l) = landmarks_ground_truth(l).id;
endfor

%% ------------------------------------ %%
%% ----------CALIBRATION INIT-----------%%
%% ------------------------------------ %%

dataset='Datasets';
guess_algo='latemax'; 
plots_path='plots';

variance_threshold = 0.6; 
max_range_obs = 6; 
l_noise = 9e-3; % Noise parameter for lateration algorithm 

global num_poses = length(poses);
global num_lms=length(available_landmarks)
global num_trs=length(transitions)
global obs_num=length(observations_ground_truth)

data="poses"; #(write "poses" or "edges")
method="edges";
tf_vec = [-0.2 -1.8 0.22]; % Transformation vector for a better visualization of the result (simulated calibration)

% do this to calibrate the odometry. With these data, you have coincident actual and ground truth transitions, so it cannot be done
% we cannot use the same measurements derived from pose-pose displacements of this dataset, we have to retrieve them separately
%{
edges=transitions;
trajectory_edges=build_obs(edges, transitions_ground_truth, num_trs);
Sensor_bias=ls_calibrate_odometry(trajectory_edges);
tf_vec2=t2v(Sensor_bias)
pause
%}

%
%% ----------------STEP 1------------------- %%
%%--------BUILD XR_GUESS and XR_true-------- %%

% the graph edges are the pose-pose relative displacements of the robot X^i_{i-1}, with X=(x,y,θ)
% encapsulated in the variable transitions

XR_guess = zeros(3, 3, num_poses);
XR_true = zeros(3, 3, num_poses);

%% Case 1: odometry without graph edges

if data=="poses"
  XR_true=compute_odometry_trajectory(poses_ground_truth);
  XR_guess=compute_odometry_trajectory(poses);
end

%% case 2: odometry using graph edges

if data=="edges"
   XR_true=compute_odometry_trajectory(poses_ground_truth);
   [XR_guess,poses_xy]=poses_odom_SE2(poses, transitions, num_poses, num_trs);
end


%% ----------------STEP 2------------------- %%
%%--------BUILD XL_GUESS and XL_true-------- %%

% get an INITIAL GUESS of the landmarks given their id, poses and obs from init_guess_dataset
[landmarks, obs_new] = get_initial_guess(available_landmarks, poses, observations, 
                                         guess_algo, l_noise, max_range_obs, variance_threshold);

global num_landmarks = length(landmarks);
% xy position of landmarks
XL_guess = zeros(2,num_landmarks);
XL_true = zeros(2,num_landmarks);
for l=1:num_landmarks
	XL_guess(1:2,l) = landmarks(l).landmark_position;
	XL_true(1:2,l) = [landmarks_ground_truth(l).x_pose; landmarks_ground_truth(l).y_pose];
endfor

% Evaluating how good is your XL guess
eval_guess = initial_guess_matching(XL_true,XL_guess);
fflush(stdout);

%% ----------------STEP 3------------------- %%
%%-----------GET THE OBSERVATIONS----------- %%

[Z,associations]=get_obs_ass(observations); %careful

%

%%-----------------------------------------%%
%%-------------OPTIMIZATION----------------%%
%%-----------------------------------------%%

printf("\n\n%% ------------------------------------ %%\n%% ------ Starting Least Squares ------ %%\n%% ------------------------------------ %%\n\n");
fflush(stdout);

% test for the calibration part
%{
for i=1:num_poses
        XR_guess2(:,:,i) = v2t(tf_vec2) * XR_guess(:,:,i);
        XR_guess(:,:,i) = v2t(tf_vec) * XR_guess(:,:,i);		
endfor
%}

%
num_iterations = 10;
damping = 0.01;
kernel_threshold = 1.0;

## save and load variables ##
save("Datasets/variables", "XR_guess", "XL_guess", "XL_true", "XR_true", "Z", "associations", "eval_guess");
%{
S=load("Datasets/variables");
XR_guess=S.XR_guess;
XR_true=S.XR_true;
XL_guess=S.XL_guess;
XL_true=S.XL_true;
global num_landmarks=size(XL_guess, 2)
Z=S.Z;
associations=S.associations;
eval_guess=S.eval_guess;
%}

method
if method=="poses" 
disp("poses");
[XR, XL, chi_stats, num_inliers]=doLeastSquares(XR_guess, XL_guess, Z, 
							associations, 
							num_poses, 
							num_landmarks, 
							num_iterations=10, 
							damping=0.15, 
							kernel_threshold=0.2);
end

if method=="edges"
T=preprocessTransitions(transitions, poses);
[XR, XL, chi_stats, num_inliers, e_record]=doLeastSquaresEdges(XR_guess, XL_guess, T, poses, Z, 
							associations, 
							num_trs,
							num_poses, 
							num_landmarks, 
							num_iterations=10, 
							damping=0.15, 
							kernel_threshold=0.2);
end

printf("\n\n%% ------------------------------------ %%\n%% ------- Optimization DONE --------- %%\n%% ------------------------------------ %%\n\n");
fflush(stdout);
pause(1);

% Now we evaluate the correction in the same way as for the initial guess
eval_correction = initial_guess_matching(XL_true,XL);
save("Datasets/solution", "eval_correction", "XR", "XL", "chi_stats", "num_inliers", "e_record");
fflush(stdout);
printf("\nAt guess stage %i were verified landmarks\nAt correction stage %i are verified landmarks\n",[eval_guess,eval_correction]);
fflush(stdout);

%

%%-----------------------------------------%%
%%---------------TRAJECTORY----------------%%
%%-----------------------------------------%%

% test for calibrated odometry
%{
figure();
%
hold on;
title("Robot trajectory");

r = 1;
c = num_poses;
    
plot(reshape(XR_guess(1,3,:), r, c), reshape(XR_guess(2,3,:), r, c), 'g-', 'linewidth', 3);
plot(reshape(XR_guess2(1,3,:), r, c), reshape(XR_guess2(2,3,:), r, c), 'r-', 'linewidth', 3);
plot(reshape(XR_true(1,3,:), r, c), reshape(XR_true(2,3,:), r, c), 'k-', 'linewidth', 3);
legend("initial guess", "initial guess 2", "ground truth");
legend("initial guess", "ground truth");
pause

%}

% check error history over iterations and transitions
%{
figure()
plot(e_record(1,:), 'r-', 'linewidth', 3);
pause
plot(e_record(end,:), 'r-', 'linewidth', 3);
pause
%}

% XR=0;
% Map plot: Landmarks IG, OPT, GT
plot_map(XL_guess, XL, XL_true, plots_path, guess_algo, dataset);
pause

% Trajectory plot: Poses IG, OPT, GT
plot_traj(XR_guess, XR, XR_true, num_poses, plots_path, guess_algo, dataset);
pause

% Trajectory plot w/ simulated odometry calibration: Poses IG, OPT, GT
plot_cal_traj(XR_guess, XR, XR_true, num_poses, tf_vec, plots_path, guess_algo, dataset);
pause
%

%{
figure()
title("FINAL RESULTS");
hold on;
grid minor;
xlim([-15,15]);
ylim([-15,10]);
plot(XR_true(1,3,:),XR_true(2,3,:), 'g-.', 'linewidth', 3);
#hold on;
plot(XR_guess(1,3,:),XR_guess(2,3,:), 'k-.', 'linewidth', 3);
#legend("ground truth", "guess");
#hold on;
plot(XR(1,3,:),XR(2,3,:), 'b-', 'linewidth', 3);
legend("ground truth","guess","correction");
pause
%}