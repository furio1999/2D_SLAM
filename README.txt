Project 03 - Range only SLAM

The dataset is composed by a g2o file which contain poses and range observations. The file contain also odometry edges that are used to construct the initial guess
for the problem. The use of those edges is your own choice if you think it's necessary.

Hint :
    
    - Parse the whole dataset and initialize the landmarks by using at least 3 range observations with the proper parallax
    - Setup a LS optimization that involves all the poses and landmarks that you have initialized
    - In the range edges the ID of the landmarkd is reported, use it to identify them for both the initialization and the global optimization
         
        -EDGE_RANGE_SE2_XY id_pose id_landmark range

Expected output :
  - Robot trajectory and Map


For any question concerning the initialization part feel free to contact me
