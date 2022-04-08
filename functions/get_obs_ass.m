function [Z,associations]=get_obs_ass(observations)

num_measurements = 0;
for o=1:length(observations)
	% num_measurements is the total amount of measuremtns you got
	num_measurements += length(observations(o).observation);
endfor

% building associations and Z
% associations is a 2x num_measurements vmatrix in which we have [ID_pose, ID_landmark]
% so for each pose you have all the landmark_id you see from there

% Z is just all the ranges you measure
associations = zeros(2,num_measurements);
Z = zeros(1,num_measurements);

% init idices
poseID = 1;
measure_num = 1;
% get a ordered list of landmark_id and their order position
allIDs = getIDs_and_index(observations);

for o=1:length(observations)
		
		for i=1:length(observations(o).observation)
		
			associations(1,measure_num) = poseID;
			% real id must be converted in [1,..61] to match ICP indices
			% 90 is not acceptable for example (see getIDs_and_index.m)
			real_id = observations(o).observation(i).id;
			associations(2,measure_num) = get_formatted_index(allIDs,real_id);

			Z(1,measure_num) = observations(o).observation(i).range;
			
			measure_num+=1;
		endfor
		poseID+=1;
	
endfor