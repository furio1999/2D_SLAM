source 'tools/utilities/geometry_helpers_2d.m'
warning('off', 'all');


% Plots the map with landmarks IG/OPT/GT
function plot_map(ig, opt, gt, path, algo, name)
    figure();
    title("Landmarks World Map");
    hold on;

    plot(ig(1,:), ig(2,:),'g*', "linewidth", 2);
    plot(opt(1,:), opt(2,:), 'bx', "linewidth", 2);
    plot(gt(1,:), gt(2,:), 'k*', "linewidth", 2);
    legend("Initial guess","LS optimization","Ground truth");

    print(fullfile(path, cstrcat(name, "_", algo, "_map_plot.png")));

end


% Plots the trajectories IG/OPT/GT
function plot_traj(ig, opt, gt, poses_num, path, algo, name)
    figure();
    hold on;
    title("Robot trajectory");

    r = 1;
    c = poses_num;
    
    plot(reshape(ig(1,3,:), r, c), reshape(ig(2,3,:), r, c), 'g-', 'linewidth', 3);
    plot(reshape(opt(1,3,:), r, c), reshape(opt(2,3,:), r, c), 'b-', 'linewidth', 3);
    plot(reshape(gt(1,3,:), r, c), reshape(gt(2,3,:), r, c), 'k-', 'linewidth', 3);
    legend("Initial guess","LS optimization","Ground truth");

    print(fullfile(path, cstrcat(name, "_", algo, "_traj_plot.png")));

end


% Plots the calibrated trajectories IG/OPT/GT
function plot_cal_traj(ig, opt, gt, poses_num, tf, path, algo, name);
    % Apply transformation vector
    for i=1:poses_num
        ig(:,:,i) = v2t(tf) * ig(:,:,i);
        opt(:,:,i) = v2t(tf) * opt(:,:,i);
    endfor

    figure();
    hold on;
    title("Robot trajectory + Simulated odometry calibration");

    r = 1;
    c = poses_num;

    plot(reshape(ig(1,3,:), r, c), reshape(ig(2,3,:), r, c), 'g-', 'linewidth', 3);
    plot(reshape(opt(1,3,:), r, c), reshape(opt(2,3,:), r, c), 'b-', 'linewidth', 3);
    plot(reshape(gt(1,3,:), r, c), reshape(gt(2,3,:), r, c), 'k-', 'linewidth', 3);
    legend("Initial guess","LS optimization","Ground truth");

    print(fullfile(path, cstrcat(name, "_", algo, "_cal_traj_plot.png")));

end