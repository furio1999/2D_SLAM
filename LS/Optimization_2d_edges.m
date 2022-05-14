source "tools/utilities/geometry_helpers_2d.m"
source "LS/LS_Utils.m"

# implementation of the optimization loop with robust kernel
# applies a perturbation to a set of landmarks and robot poses
# input:
#   XR: the initial robot poses (4x4xnum_poses: array of homogeneous matrices)
#   XL: the initial landmark estimates (3xnum_landmarks matrix of landmarks)
#   Z:  the measurements (3xnum_measurements)
#   associations: 2xnum_measurements. 
#                 associations(:,k)=[p_idx,l_idx]' means the kth measurement
#                 refers to an observation made from pose p_idx, that
#                 observed landmark l_idx
#   num_poses: number of poses in XR (added for consistency)
#   num_landmarks: number of landmarks in XL (added for consistency)
#   num_iterations: the number of iterations of least squares
#   damping:      damping factor (in case system not spd)
#   kernel_threshod: robust kernel threshold

# output:
#   XR: the robot poses after optimization
#   XL: the landmarks after optimization
#   chi_stats: array 1:num_iterations, containing evolution of chi2
#   num_inliers: array 1:num_iterations, containing evolution of inliers

function [XR, XL, chi_stats, num_inliers, e_record]=doLeastSquaresEdges(XR, XL, edges, poses, Z, 
							associations, 
              num_edges,
							num_poses, 
							num_landmarks, 
							num_iterations, 
							damping, 
							kernel_threshold)
  global pose_dim;
  global landmark_dim;

  chi_stats=zeros(1,num_iterations);
  chi_stats_2=zeros(1,num_iterations);
  num_inliers=zeros(1,num_iterations);
  # size of the linear system
  system_size=pose_dim*num_poses+landmark_dim*num_landmarks;  % 3*301 + 2*61
  first_pose_id=poses(1).id

  for (iteration=1:num_iterations)
	
		printf('\nIterations processed : %i%%',iteration*10);
		fflush(stdout);

	
		% we iterate num_iterations times
	
    H=zeros(system_size, system_size);
    b=zeros(system_size,1);
		
    chi_stats(iteration)=0; % vettore degli errori
    chi_stats_2(iteration)=0;

    %pose-pose
    e_record=zeros(num_iterations,length(edges));
    for trans_num=1:length(edges)

      % in main, v2t(edge) to give in form X_i^{i+1} (R|t)
      transition=edges(trans_num); 
      v=transition.v;
      idx_i=(transition.id_from - first_pose_id)+1; #trans_num
      idx_j=(transition.id_to - first_pose_id)+1;   #trans_num+1

    
      Xri=XR(:,:,idx_i);
      #to refine, you can store pose-pose indices and retrive there to find Xj
      # in this case, there is a 1-1 correspondance with the successive pose index (in edges, every pose in relative to the previous one)
      Xrj=XR(:,:, idx_j);
      [e2,Ji,Jj]=errorAndJacobianOdometry(Xri, Xrj, v);
      e_record(iteration, trans_num)=norm(e2);
      
      chi2=e2'*e2;
      if (chi2>kernel_threshold)
      	e2*=sqrt(kernel_threshold/chi2);
      	chi2=kernel_threshold;
      else
      	num_inliers(iteration)++;
      endif;
      chi_stats_2(iteration)+=chi2;

      Hii=Ji'*Ji;
      Hij=Ji'*Jj;
      Hjj=Jj'*Jj;
      bi=Ji'*e2;
      bj=Jj'*e2;

      idx_mat_i=poseMatrixIndex(idx_i, num_poses, num_landmarks);
      idx_mat_j=poseMatrixIndex(idx_j, num_poses, num_landmarks);
      
      H(idx_mat_i:idx_mat_i+pose_dim-1, idx_mat_i:idx_mat_i+pose_dim-1)+=Hii;
      H(idx_mat_i:idx_mat_i+pose_dim-1, idx_mat_j:idx_mat_j+pose_dim-1)+=Hij;
      H(idx_mat_j:idx_mat_j+pose_dim-1, idx_mat_i:idx_mat_i+pose_dim-1)+=Hij';
      H(idx_mat_j:idx_mat_j+pose_dim-1, idx_mat_j:idx_mat_j+pose_dim-1)+=Hjj;

      b(idx_i:idx_i+pose_dim-1)+=bi;
      b(idx_j:idx_j+pose_dim-1)+=bj;

    endfor
		
		
    % pose - landmark
    for (measurement_num=1:size(Z,2))
			
			% for all the measure we have 1780 
			
			% get pose index and landmark measured
      pose_index=associations(1,measurement_num);
      landmark_index=associations(2,measurement_num);
			
			% take the measure, robot and landamrk coordinates 
      z=Z(:,measurement_num);
      Xr=XR(:,:,pose_index);
      Xl=XL(:,landmark_index);
			
      [e,Jr,Jl] = errorAndJacobian(Xr, Xl, z);
      chi=e'*e;
			
      if (chi>kernel_threshold)
      	e*=sqrt(kernel_threshold/chi);
      	chi=kernel_threshold;
      else
      	num_inliers(iteration)++;
      endif;
      chi_stats(iteration)+=chi;

      Hrr=Jr'*Jr;
      Hrl=Jr'*Jl; 
      Hll=Jl'*Jl;
			
      br=Jr'*e;
      bl=Jl'*e;
		
      pose_matrix_index=poseMatrixIndex(pose_index, num_poses, num_landmarks);
      landmark_matrix_index=landmarkMatrixIndex(landmark_index, num_poses, num_landmarks);

			H(pose_matrix_index:pose_matrix_index+pose_dim-1,
			pose_matrix_index:pose_matrix_index+pose_dim-1)+=Hrr;

			H(pose_matrix_index:pose_matrix_index+pose_dim-1,
			landmark_matrix_index:landmark_matrix_index+landmark_dim-1)+=Hrl;

			H(landmark_matrix_index:landmark_matrix_index+landmark_dim-1,
			landmark_matrix_index:landmark_matrix_index+landmark_dim-1)+=Hll;

			H(landmark_matrix_index:landmark_matrix_index+landmark_dim-1,
			pose_matrix_index:pose_matrix_index+pose_dim-1)+=Hrl';

			b(pose_matrix_index:pose_matrix_index+pose_dim-1)+=br; %'
			b(landmark_matrix_index:landmark_matrix_index+landmark_dim-1)+=bl;

    endfor
		
    H+=eye(system_size)*damping;
    dx=zeros(system_size,1);

    % we solve the linear system, blocking the first pose
    % this corresponds to "remove" from H and b the locks
    % of the 1st pose, while solving the system

    dx(pose_dim+1:end)=-(H(pose_dim+1:end,pose_dim+1:end)\b(pose_dim+1:end,1));
    [XR, XL]=boxPlus(XR,XL,num_poses, num_landmarks, dx);
  endfor
endfunction


# plot landmarks and poses
#
#
#
function i = plotState(XL, XL_guess, XL_gt)
#plot landmarks
hold on;
plot3(XL(1,:),XL(2,:),XL(3,:),'b*',"linewidth",2);
hold on;
plot3(XL_guess(1,:),XL_guess(2,:),XL_guess(3,:),'ro',"linewidth",2);
hold on;
plot3(XL_gt(1,:),XL_gt(2,:),XL_gt(3,:),'g*',"linewidth",2);
hold on;
legend("estimate","initial guess","ground truth")
i = 1;
endfunction
