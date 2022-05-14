%(minimal) size of pose and landmarks
global pose_dim=3;
global landmark_dim=2;

function v_idx=poseMatrixIndex(pose_index, num_poses, num_landmarks)
  global pose_dim;
  global landmark_dim;

  if (pose_index>num_poses)
    v_idx=-1;
    return;
  endif;
  v_idx=1+(pose_index-1)*pose_dim;
endfunction;



function v_idx=landmarkMatrixIndex(landmark_index, num_poses, num_landmarks)
  global pose_dim;
  global landmark_dim;
  if (landmark_index>num_landmarks)
    v_idx=-1;
    return;
  endif;
  v_idx=1 + (num_poses)*pose_dim + (landmark_index-1) * landmark_dim;
endfunction;

% Retrieves indices of poses/landmarks in the Hessian matrix
function [p_matrix_idx, l_matrix_idx] = matrix_index(p_index, lm_index, num_poses, lms_num)
    global pose_dim;
    global lm_dim;
    
    if p_index > num_poses || p_index <= 0
        p_matrix_idx = -1;
    else
        p_matrix_idx = 1 + (p_index-1)*pose_dim;
    endif


    if lm_index > lms_num || lm_index <= 0
        l_matrix_idx = -1;
    else
        l_matrix_idx = 1 + (num_poses)*pose_dim + (lm_index-1)*lm_dim;
    endif
    
end

# error and jacobian of a measured landmark
# input:
#   Xr: the robot pose (4x4 homogeneous matrix)
#   Xl: the landmark pose (3x1 vector, 3d pose in world frame)
#   z:  measured position of landmark
# output:
#   e: 3x1 is the difference between prediction and measurement
#   Jr: 3x6 derivative w.r.t the error and a perturbation on the
#       pose
#   Jl: 3x3 derivative w.r.t the error and a perturbation on the
#       landmark
function [e,Jr,Jl]=errorAndJacobian_2(Xr,Xl,z)

	[xr,yr] = deal(Xr(1,3), Xr(2,3));
	[xl,yl] = deal(Xl(1), Xl(2));
	
	z_hat = sqrt( (xr-xl)^2 + (yr-yl)^2 ) ;
	 
	Jr = [ (xr-xl)/z_hat    (yr-yl)/z_hat    0];
	Jl = [-(xr-xl)/z_hat   -(yr-yl)/z_hat    ];

	if(z_hat==0)
		%disp("z_hat è 0")
		z_hat = 0.01;
	endif

	e=(z_hat-z);
endfunction;

# another error and jacobian function

# error and jacobian of a measured landmark
# input:
#   Xr: the robot pose (4x4 homogeneous matrix)
#   Xl: the landmark pose (3x1 vector, 3d pose in world frame)
#   z:  measured position of landmark
# output:
#   e: 3x1 is the difference between prediction and measurement
#   Jr: 3x6 derivative w.r.t the error and a perturbation on the
#       pose
#   Jl: 3x3 derivative w.r.t the error and a perturbation on the
#       landmark

function [e,Jr,Jl]=errorAndJacobian(Xr,Xl,z)
   R=Xr(1:2,1:2);
   t=Xr(1:2,3);
   p_hat = R'*(Xl -t);  % Xr^(-1)*Xl
   z_hat=norm(p_hat); % prediction
   e=z_hat-z;
	 
   Jr = zeros(1,3);
   J_icp = zeros(2,3);
   J_icp(1:2,1:2) = -R';
	 
   J_icp(1:2,3) = R'*[0 1;-1 0]*Xl;
   Jr = (1/norm(p_hat))*p_hat'* J_icp;
   Jl= (1/norm(p_hat))*p_hat'*R';
	 
endfunction;

function [e, Ji, Jj]=errorAndJacobianEdge(Xi,Xj,z)
   Ri=Xi(1:2,1:2);
   ti=Xi(1:2,3);
   Rj=Xj(1:2,1:2);
   tj=Xj(1:2,3);
   z_hat=inv(Xi)*Xj;
   z_hat=reshape(z_hat', 1, []);
   v=z.v;
   z=v2t(v);
   z=reshape(z',1,[]);

   e=z_hat-z;
   % prebuild derivative here?
   c=Rj(1,1);
   s=Rj(1,2);
   dRj=[-s, -c;
          c,  -s];

   dgtheta=[Ri'*dRj*Rj , Ri'*dRj*tj];
   dgtheta=reshape(dgtheta',1,[]);
   dgx=[zeros(2,2), Ri'*[1 0]'];
   dgx=reshape(dgx', 1, []);
   dgy=[zeros(2,2), Ri'*[0 1]'];
   dgy=reshape(dgy', 1, []);

   # 6x3 matrix
   Jj=[dgx', dgy', dgtheta']; 
   Ji=-Jj;

endfunction


# implementation of the boxplus
# applies a perturbation to a set of landmarks and robot poses
# input:
#   XR: the robot poses (4x4xnum_poses: array of homogeneous matrices)
#   XL: the landmark pose (3xnum_landmarks matrix of landmarks)
#   num_poses: number of poses in XR (added for consistency)
#   num_landmarks: number of landmarks in XL (added for consistency)
#   dx: the perturbation vector of appropriate dimensions
#       the poses come first, then the landmarks
# output:
#   XR: the robot poses obtained by applying the perturbation
#   XL: the landmarks obtained by applying the perturbation

function [XR, XL]=boxPlus(XR, XL, num_poses, num_landmarks, dx)
  global pose_dim;
  global landmark_dim;

  for(pose_index=1:num_poses)
    pose_matrix_index=poseMatrixIndex(pose_index, num_poses, num_landmarks);
    dxr=dx(pose_matrix_index:pose_matrix_index+pose_dim-1);
    XR(:,:,pose_index)=v2t(dxr)*XR(:,:,pose_index);
  endfor;

  for(landmark_index=1:num_landmarks)
    landmark_matrix_index=landmarkMatrixIndex(landmark_index, num_poses, num_landmarks);
    dxl=dx(landmark_matrix_index:landmark_matrix_index+landmark_dim-1,:);
    XL(:,landmark_index)+=dxl(1:2); % edit (1:2)
  endfor;
  
endfunction;

function [ep,J1,J2]=errorAndJacobianOdometry(P1,P2,trans)

  p1=t2v(P1);
  p2=t2v(P2);

  %some components to avoid giant jacobian
  x1=p1(1); y1=p1(2); theta1=p1(3);
  x2=p2(1); y2=p2(2); theta2=p2(3);
  c1=cos(theta1); s1=sin(theta1); c2=cos(theta2); s2=sin(theta2); 
  s12=sin(theta1 - theta2); c12=cos(theta1 - theta2);

  %odometry measurement
  #v=trans(1); t=trans(3);
  v1=trans(1); v2=trans(2); t=trans(3);

  %estimate next pose using current pose and odometry measurements
  #p2_est=[x1+v*c1 y1+v*s1 theta1+t]';
  p2_est=[x1+v1 y1+v2 theta1+t]';
  P2_est=v2t(p2_est);

  %error= displacement I have - displacement I measure 
  #actual_disp=inv(P1)*P2;
  actual_disp=v2t(trans); #improvement from 23 correct to 39 correct out of 61
  measure_disp=inv(P1)*P2_est;
  ep=actual_disp - measure_disp;

  
  %reshaping error to have a vector so it cna be used to compute J'*e and so B
  %useless values are removed and it's possible to compute the jacobian using symbolic operations
  ep=reshape(ep(1:2,:),6,1);

  
  %derivative of error (6x1) w.r.t previoud and next state (6 variables)
  %so jacobian is 6x6
  Jp=[
    0,     0,              -s12,               0,  0,  s12;
    0,     0,              -c12,               0,  0,  c12;
    0,     0,               c12,               0,  0, -c12;
    0,     0,              -s12,               0,  0,  s12;
    -c1, -s1, y2*c1 - y1*c1 + x1*s1 - x2*s1,  c1, s1,   0;
     s1, -c1, x1*c1 - x2*c1 + y1*s1 - y2*s1, -s1, c1,   0];

%6x3
J1=Jp(:,1:3);

%6x3
J2=Jp(:,4:6);

endfunction

%boxplus operator to add perturbation in manifold space
function [XR, XL]=boxPlus2(XR, XL, num_poses, num_landmarks, dx)
  global pose_dim;
  global landmark_dim;
  for(pose_index=1:num_poses)
    pose_matrix_index=poseMatrixIndex(pose_index, num_poses, num_landmarks);
    dxr=dx(pose_matrix_index:pose_matrix_index+pose_dim-1);
    %for poses we can't do simple addition, here is why boxplus is needed
    XR(:,:,pose_index)=v2t(dxr)*XR(:,:,pose_index);
    
  endfor;
  for(landmark_index=1:num_landmarks)
    %this is needed because some landmarks were not considered
    if landmark_index<=length(XL) && XL(landmark_index)
      landmark_matrix_index=landmarkMatrixIndex(landmark_index, num_poses, num_landmarks);
      dxl=dx(landmark_matrix_index:landmark_matrix_index+landmark_dim-1,:);
      %nothing special needed for landmark
      XL(:,landmark_index)+=dxl(1:2);
    endif
  endfor;
endfunction;