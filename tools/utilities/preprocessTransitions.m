function pose=poseFromId(id,poses)
    %cost O(1)
    %this work because each pose has id = 1 + previous id 
    first_id=poses(1).id;
    pose=poses(id-first_id+1);

    %{
    %real loop used to check for real if they correspond, now not needed
    %it will be needed if the pose has not id=previous_id+1
    for i=1:length(poses)
        if poses(i).id==id
            pose=poses(i);
        endif
    endfor
    %}
endfunction


function T=preprocessTransitions(transitions,poses)

    for i=1:length(transitions)

        pose_prev=poseFromId(transitions(i).id_from,poses);
        pose_next=poseFromId(transitions(i).id_to,poses);
        T(end+1).pose_prev=v2t([pose_prev.x pose_prev.y pose_prev.theta]');
        T(end).pose_new=v2t([pose_next.x pose_next.y pose_next.theta]');
        T(end).v=transitions(i).v;

        T(end).id_from=transitions(i).id_from;
        T(end).id_to=transitions(i).id_to;
    endfor

endfunction