# Spatiotemporal-Colorization-of-Video
Matlab  implementation of the paper [Spatiotemporal Colorization of Video Using 3D Steerable Pyramids](https://ieeexplore.ieee.org/document/7428858/)

Scripts `level2marked_3d.m` and `level1marked_3d.m` generates partially colorized video volumes at levels 2 and 1 as referenced in the paper. Script `level0marked_3d.m` generates a fully colored video volume. 

Example implementation colors 10 frames from the *Foreman* sequence where the first and last frames are taken as scribbled keyframes. 

Run scripts in following order:  
`filters_3d.m`  
`edges.m`  
`level2marked_3d.m`  
`level1marked_3d.m`  
`level0marked_3d.m`  
