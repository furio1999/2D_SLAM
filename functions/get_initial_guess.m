%% ------------------------------------ %%
%% ------ MAIN FUNCTION --------------- %%
%% ------------------------------------ %%

% inputs
% 	available_landmarks is a lsit of all landmaks id
% 	poses are all the poses
% 	observations are all the observations

function [my_initial_condition, obs_aux]= get_initial_guess(available_landmarks, poses, observations, algorithm, 
                                                  l_noise, max_range_obs, variance_threshold)

	my_initial_condition = [];
	obs_aux=[];

	printf('\nDataset processed : %d%%',0);
	for l=1:length(available_landmarks)

		if mod( round(l/length(available_landmarks)*100), 10)<1.1
			printf('\nDataset processed : %d%%',[round(l/length(available_landmarks)*100)]);
			fflush(stdout);
		endif
		
		% now i'm optimizing this landmark id
		searched_landmark = available_landmarks(l);
		
		% retrieve landmark info from the dataset 
		% you get landmark_position, range, id, how many times you got it
		landmark_info = find_landmark_references(searched_landmark,poses,observations, max_range_obs, variance_threshold); 
		obs_aux(l).observation=landmark_info;
		% start optimization process for this landmark based on the info you have
		if landmark_info(1).id ~= -1
		landmark_pos = triangulate_landmark_pos(landmark_info);
		end
		
		% update initial_condition
		my_initial_condition(end+1).id = searched_landmark; %non sono entrambi end o end+1??
		my_initial_condition(end).landmark_position = landmark_pos;

	endfor
	# printf('\nDataset processed : %d%% \n',100);
	# fflush(stdout);
	pause(0.1);


end

%% ----------------------------------------- %%
%%  GETTING INFO FROM THE SEARCHED LANDMRK - %%
%% ----------------------------------------- %%
function id_pose_range = find_landmark_references(searched_landmark,poses,observations, max_range_obs, variance_threshold)
	
	% THESE ARE THE INFOS I WILL SEND BACK FOR EACH SEARCHED_LMK
	id_pose_range = [];
	land_poses_id = [];
	land_poses_xy = [];
	land_ranges = [];
	count = 0;

	for s=1:length(observations)
	% for each sample in observations
		sample = observations(s);
		
		for o=1:length(sample.observation)
		
			landmark_and_range = sample.observation(o);
			landmark_seen = landmark_and_range.id;
			
			if landmark_seen == searched_landmark
				count = count+1;
				% add pose and range to the lists
				land_poses_id(end+1) = sample.pose_id;
				land_ranges(end+1) = landmark_and_range.range;
			
			endif
		endfor
	endfor

	% printf('\nlandmark %i has been seen from %i poses\n', [searched_landmark, count]);
	
	% retrieve each xy couple from any pose_id in land_poses_id in poses from G2o

	for i=1:length(land_poses_id)
		for p=1:length(poses)

			if poses(p).id == land_poses_id(i)
			% tell me xy
				land_poses_xy(end+1).x = poses(p).x;
				land_poses_xy(end).y = poses(p).y;
			endif

		endfor
	endfor
	
	% printf('\nland_poses_id has length %i\n', length(land_poses_id));
	% printf('\nland_poses_xy has length %i\n', length(land_poses_xy));
	% printf('\nland_ranges has length %i\n', length(land_ranges));
		
	% save a prettier data structure
	r=land_ranges;
	variance=var(r);
    step_s=1;
	idx_s=1;
	%
	if count<3
       id_pose_range(1).id=-1;
	elseif length(count) > max_range_obs  && variance < variance_threshold
	        printf("variance regulation")
	        count=max_range_obs;
	        step_s = floor(length(r) / max_range_obs);
            idx_s = 1;
	        for e=1:max_range_obs
		       id_pose_range(end).id = land_poses_id(idx_s);
		       id_pose_range(end).pose = land_poses_xy(idx_s);
		       id_pose_range(end).range=land_ranges(idx_s);
			   idx_s+=step_s;
			endfor
	else	
	%	
	  for e=1:count
		 id_pose_range(e).id = land_poses_id(e);
		 id_pose_range(e).pose = land_poses_xy(e);
		 id_pose_range(e).range = land_ranges(e);
		
	    endfor
	end
	% add info how many times you've seen it (just first position)
	id_pose_range(1).counts = count;
end

% lateration 

function landmark_position = triangulate_landmark_pos(info)
	
	num_points = length(info);
	last_x = info(num_points).pose.x;
	last_y = info(num_points).pose.y;
	r_last = info(num_points).range;

	A = zeros(num_points-1,2);
	b = zeros(num_points-1,1);
	
	for(i=1:num_points-1)
		xi = info(i).pose.x;
		yi = info(i).pose.y;
		ri = info(i).range;

		A(i,1) = xi - last_x;
		A(i,2) = yi - last_y;
		b(i) = xi^2-last_x^2 + yi^2-last_y^2 + r_last^2-ri^2;
	endfor
	b = 0.5*b;
	damp_pinv = inv(A'*A+eye(2)*0.001)*A';

	landmark_position = damp_pinv*b;
end