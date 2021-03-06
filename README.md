# 2D_SLAM
Project for the course in [Probabilistic Robotics](https://sites.google.com/diag.uniroma1.it/probabilistic-robotics-2021-22/home) at Sapienza University of Rome\
\
<a href="https://www.dis.uniroma1.it/"><img src="http://www.dis.uniroma1.it/sites/default/files/marchio%20logo%20eng%20jpg.jpg" width="500"></a>

## Background
How do we navigate an autonomous robot in an indoor scenario? 
One of the main techniques is to fill the environment with known landmarks, estimating the robot position by simple heuristics, such as latemax.
But this naive approach leads to the common phenomenon of "dead recking", which is the generation of a trajectory which diverges from the actual one.
In this work we will show how to handle this scenario using a Least Square Approach to optimize the trajectory, closing the trajectory loop.

## Material
- Dataset (in .g2o format) composed of : 2D range measurements, landmarks positions (only ground-truth), edge transitions, Robot Poses
- Wheeled Mobile Robot to provide the data

## Approach
choose whether to use absolute positions or pose-pose relative displacements
```
# in main.m
data="poses" # compute initial trajectory using absolute poses
data="edges" # compute initial trajectory using transition measurements
method="poses" # Do Least Squares on range measurements only
method="edges" # Do least Squares using both range and odometry measurements
```

Launch the main file
```
octave main.m
```
## Range-Only
### Landmark Map
Firstly, we generate a map of the real and estimated landmark positions
![Datasets_naive_map_plot](https://user-images.githubusercontent.com/63920397/163470460-0323ed6d-b63f-46ba-b879-2e283c8630ff.png)

### Trajectories
Then, we plot the initial guess and the optimzed trajectory.
![Datasets_naive_traj_plot](https://user-images.githubusercontent.com/63920397/163470490-d0d9fa7d-7b14-4d61-9414-6983f1a34ffc.png)

As we see, it's better to apply a lest squares approach to generate a more accurate trajectory. In particular, given the state (x,y,θ) ∈ SE(2), we used a Least Squares approach equipped with a BoxPlus operator to update the optimal solution respecting the manifold gemometry relations.
In the next paragraph, we show more clearly the difference between the optimized and naive trajectories.
### Trajectories with Calibrated Odometry
From the previous graph, it seems that the computed optimized trajectory follows the ground-truth, except for a rototranslation.
Given the hypothesis of calibrated odometry, I provide a transformation vector to map the computed trajecories into their original pose.
![Datasets_naive_cal_traj_plot](https://user-images.githubusercontent.com/63920397/163470520-119055b9-f461-4d72-b3a4-c034b566e139.png)

## Range and Odometry
Here I show the results using odometry measurements too
### Landmark Map
![Datasets_latemax_map_plot](https://user-images.githubusercontent.com/63920397/168425780-f8ed3bba-eadc-4753-8ae1-62e36c438c1d.png)

### Trajectories
![Datasets_latemax_traj_plot](https://user-images.githubusercontent.com/63920397/168425798-3ea61b47-65f0-42b0-85d8-fc21c7290054.png)
