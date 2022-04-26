source 'tools/utilities/geometry_helpers_2d.m'

function T=compute_odometry_trajectory(U)
	T=zeros(3,3,size(U));

	for p=1:size(U, 2)   
	% robot coordinates for that pose p
	[x,y,theta] = deal(U(p).x, U(p).y, U(p).theta);
	u=[x,y,theta]';

	% R
	T(1,1,p) = cos(theta);
	T(1,2,p) = -sin(theta);
	T(2,1,p) = sin(theta);
	T(2,2,p) = cos(theta);
	
	% t
	T(1,3,p) = x;
	T(2,3,p) = y;
	
	% 1
	T(3,3,p) = 1;
	end
end

% Computes ground truth and initial guess SE2 poses with odometry
function [Xr_ig_odom, poses_ig_odom_xy] = poses_odom_SE2(poses_ig, transitions_ig, poses_num, trs_num)
    Xr_ig_odom = zeros(3, 3, poses_num);
    poses_ig_odom_xy = zeros(2, poses_num);
    
    % First pose is fixed (same for xy) -> Prior
    Xr_ig_odom(:,:,1)      = v2t([poses_ig(1).x poses_ig(1).y poses_ig(1).theta]);
    poses_ig_odom_xy(:, 1) = [poses_ig(1).x; poses_ig(1).y];
    
    for i=1:trs_num
        tr = transitions_ig(i).v;   
        u_x = tr(1);
        u_theta = tr(3); 
        
        % Actual pose of the robot
        pose_act = Xr_ig_odom(:,:,i);

        % Compute pose displacement and get next pose
        pose_d = v2t([u_x 0 u_theta]);
        pose_next = pose_act * pose_d; 
        Xr_ig_odom(:,:,i+1) = pose_next; 

        % Fill xy poses vector
        poses_ig_odom_xy(:, i+1) = pose_next(1:2, 3);
    endfor

end

function T=compute_odometry_trajectory_edges(U, trs_num)
	T=zeros(size(U,1),3);
	current_T=v2t(zeros(1,3));
	for i=1:size(U,1),
		u=U(i,1:3)';
		current_T*=v2t(u);
		T(i,1:3)=t2v(current_T)';
	end
end

function C=apply_odometry_correction(X, U)
	C=zeros(size(U,1),3);
	for i=1:size(U,1),
		u=U(i,1:3)';
		uc=X*u;
		C(i,:)=uc;
	end
end

