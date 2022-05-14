source "tools/utilities/geometry_helpers_2d.m"


function [e, Ji, Jj]=errorAndJacobianEdge(Xi,Xj,z)
   Ri=Xi(1:2,1:2);
   ti=Xi(1:2,3);
   Rj=Xj(1:2,1:2);
   tj=Xj(1:2,3);
   z_hat=inv(Xi)*Xj;
   z_hat=reshape(z_hat', 1, []);
   z=reshape(z',1,[]);

   e=z_hat-z;
   % prebuild derivative here?
   c=Rj(1,1);
   s=Rj(1,2);
   dRj=[-s, -c;
          c,  -s];
   %{
   dgtheta=[Ri'*dRj*Rj , Ri'*dRj*tj];
   dgtheta=reshape(dgtheta',1,[]);
   dgx=[zeros(2,2), Ri'*[1 0]'];
   dgx=reshape(dgx', 1, []);
   dgy=[zeros(2,2), Ri'*[0 1]'];
   dgy=reshape(dgy', 1, []);
   %}
   r=reshape(Ri'*dRj*Rj, [], 1)

   # 6x3 matrix
   #Jj=[dgx', dgy', dgtheta']
   Jj=[zeros(4,1), Ri'; r, Ri'*tj]
   Ji=-Jj;

endfunction


function plot_H(H, num_measurements, pose_dim, size_grid)
  rows=size(H,1);
  cols=size(H,2);
  result=zeros(num_measurements*size_grid, num_measurements*size_grid);

  for i=1:num_measurements
     for j=1:num_measurements
       upleft=[(i+1)*size_grid, j*size_grid];
       downleft=[(i)*size_grid, j*size_grid];
       upright=[(i+1)*size_grid, (j+1)*size_grid];
       downright=[i*size_grid, (j+1)*size_grid];
       rectangle("Position", [(i)*size_grid, j*size_grid, (i+1)*size_grid, (j+1)*size_grid], "FaceColor", "w")
       if H(i:i+pose_dim-1, j:j+pose_dim-1) != zeros(pose_dim, pose_dim)
       result(i,j)=1;
       rectangle("Position", [(i)*size_grid, j*size_grid, (i+1)*size_grid, (j+1)*size_grid], "FaceColor", "b")
       endif
      endfor
  endfor

endfunction

%{
x=[0, 0, 2 ;
    1, 5, 6];
x=reshape(x', 1, [])
pose_dim=2;
num=5;
size_grid=2;
H=zeros(num*pose_dim, num*pose_dim);
H(1:1+pose_dim-1, 1:1+pose_dim-1)=1;
H(1:1+pose_dim-1, 3:3+pose_dim-1)=1;
H(3:3+pose_dim-1, 1:1+pose_dim-1)=1
plot_H(H, num, pose_dim, size_grid);
pause
%}
z=[0 -1  5;
   1   0   6;
   0    0   1]
X(:,:, 1)=z;
X(:,:, 2)=z;
save("Datasets/prova", "X")
S=load("Datasets/prova");
X=S.X

X1=v2t([0 0 3.14]);
X2=v2t([2 6 0.6]);
#[e,Ji,Jj]=errorAndJacobianEdge(X1,X2,z);
X21=reshape(X1', 1, [])
