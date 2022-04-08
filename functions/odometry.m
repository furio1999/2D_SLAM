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