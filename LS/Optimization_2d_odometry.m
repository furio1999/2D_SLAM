source "tools/utilities/geometry_helpers_2d.m"

function X=ls_calibrate_odometry(Z)
	#accumulator variables for the linear system
	H=zeros(9,9);
	b=zeros(9,1);
	#initial solution (the identity transformation)
	X=eye(3); 
	
	#loop through the measurements and update the
	#accumulators
	size(Z,1)
	for i=1:size(Z,1)
		e=error_function(i,X,Z);
		A=jacobian(i,Z);
		H=H+A'*A;
		b=b+A'*e;
	end
	#solve the linear system
	deltaX=-H\b;
	#this reshapes the 9x1 increment vector in a 3x3 atrix
	dX=reshape(deltaX,3,3)';
	#computes the cumulative solution
	X=X+dX;
end

function e=error_function(i,X,Z)
	uprime=Z(i,1:3)';
	u=Z(i,4:6)';
	e=uprime-X*u;
end

#derivative of the error function for the ith measurement in Z
#does not depend on the state
#i:	the measuement number
#Z:	the measurement matrix
#A:	the jacobian of the ith measurement
function A=jacobian(i,Z)
	u=Z(i,4:6);
	A=zeros(3,9);
	A(1,1:3)=-u;
	A(2,4:6)=-u;
	A(3,7:9)=-u;
end
