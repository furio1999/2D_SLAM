close all
clear
clc

#import 2d geometry utils
source "tools/utilities/geometry_helpers_2d.m"
source "functions/odometry.m"
source "functions/get_initial_guess.m"
source "functions/initial_guess_matching.m"
source "functions/get_obs_ass.m"
addpath("functions")
addpath("tools/g2o_wrapper")
addpath("Datasets")
addpath("ICP")
#source "ICP/MultiICP.m"
source "ICP/multi_ICP_2d.m"
source "tools/visualization/plotter.m"

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

%% get an INITIAL GUESS of the landmarks given their id, poses and obs from init_guess_dataset
% FIELDS: 
% 	landmarks(i).id
% 	landmarks(i).landmark_position(1) and landmarks(i).landmark_position(2)
landmarks = get_initial_guess(available_landmarks, poses, observations);

%% ------------------------------------ %%
%% ----------- ICP INIT --------------- %%
%% ------------------------------------ %%

global num_poses = length(poses);
global num_landmarks = length(landmarks);

%% ----------------STEP 1------------------- %%
%%--------BUILD XR_GUESS and XR_true-------- %%

XR_guess = zeros(3, 3, num_poses);
XR_true = zeros(3, 3, num_poses);

XR_true=compute_odometry_trajectory(poses_ground_truth);
XR_guess=compute_odometry_trajectory(poses);

%% ----------------STEP 2------------------- %%
%%--------BUILD XL_GUESS and XL_true-------- %%

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

[Z,associations]=get_obs_ass(observations);

%

%%-----------------------------------------%%
%%-------------ICP OPTIMIZATION------------%%
%%-----------------------------------------%%

printf("\n\n%% ------------------------------------ %%\n%% ------ Starting ICP optimization --- %%\n%% ------------------------------------ %%\n\n");
fflush(stdout);

num_iterations = 10;
damping = 0.01;
kernel_threshold = 1.0;
[XR, XL, chi_stats, num_inliers]=doMultiICP(XR_guess, XL_guess, Z, 
							associations, 
							num_poses, 
							num_landmarks, 
							num_iterations=10, 
							damping=0.15, 
							kernel_threshold=0.2);

printf("\n\n%% ------------------------------------ %%\n%% ------- ICP optimization DONE ------ %%\n%% ------------------------------------ %%\n\n");
fflush(stdout);
pause(1);

% Now we evaluate the correction in the same way as for the initial guess
eval_correction = initial_guess_matching(XL_true,XL);
fflush(stdout);
printf("\nAt guess stage %i were verified landmarks\nAt correction stage %i are verified landmarks\n",[eval_guess,eval_correction]);
fflush(stdout);

%

%%-----------------------------------------%%
%%---------------TRAJECTORY----------------%%
%%-----------------------------------------%%

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

poses_num=length(poses);
dataset='Datasets';
guess_algo='naive'; % to be defined
plots_path='plots';
% XR=0;
% Map plot: Landmarks IG, OPT, GT
plot_map(XL_guess, XL, XL_true, plots_path, guess_algo, dataset);
pause

% Trajectory plot: Poses IG, OPT, GT
plot_traj(XR_guess, XR, XR_true, poses_num, plots_path, guess_algo, dataset);
pause

tf_vec = [-0.2 -1.8 0.22]; % Transformation vector for a better visualization of the result (simulated calibration)
% Trajectory plot w/ simulated odometry calibration: Poses IG, OPT, GT
plot_cal_traj(XR_guess, XR, XR_true, poses_num, tf_vec, plots_path, guess_algo, dataset);
pause
