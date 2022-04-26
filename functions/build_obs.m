function Z=build_obs(edges, edges_gt, nitems)
   Z=zeros(nitems, 6);
   count=0
   for i=1:nitems
      v=edges(i).v;
      v_gt=edges(i).v;
      Z(i,1:3)=v_gt;
      Z(i, 4:6)=v;
      if isequal(v,v_gt)
        count=count+1;
      end
   counter=count 
   endfor
end

%{
function trans=build_edge(poses, num_poses)
   for i=1:num_poses-1
        temp=poses(i+1)-poses(i);
   endfor

end
%}