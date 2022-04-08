function guess = myguess(landmarks, poses, observations)
        matrix=data_ass(landmarks, poses, observations)
        disp(matrix)
        #
        # landmark infos
    endfor

end

# data association
function matrix = data_ass(landmarks, poses, observations)
    lms_num = length(landmarks);
    poses_num=length(poses);

    % Vectors to map landmarks indices with
    % their ids.
    lm_indices_all = ones(10000, 1) * -1;
    lm_ids_all = ones(lms_num, 1) * -1;

    % Pair landmarks and their observers in a matrix. Access it in o(1) 
    observed_lms = ones(lms_num, poses_num) * -1;
    idx = 1;

    % Matrix construction
    first_pose_id = landmarks(1);
    for i=1:obs_num
        obs = observations(i);
        pose_id = observations(i).pose_id;
        p_idx = pose_id - first_pose_id + 1;
        
        for j=1:length(obs.observation)
            lm = obs.observation(j).id;
            rn = obs.observation(j).range;
    
            % Check lm index
            if lm_indices_all(lm) == -1
                lm_indices_all(lm) = idx;
                lm_ids_all(idx) = lm;
                idx += 1;
            endif
            
            observed_lms(lm_indices_all(lm), p_idx) = rn;

        endfor
    endfor

    % gating techniques
    

end