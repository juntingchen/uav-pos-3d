# 3D Urban UAV Relay Placement: Linear Complexity Algorithm and Analysis

This repository hosts the data and code for studying the UAV positioning problem under obstructive environment. This study focuses on the scenario of the low-altitude UAV to ground communication in a dense urban environment. There could be a lot of local structure, such as buildings and trees, that blocks the communication signal. As a result, the UAV should be optimized to smartly explore a good propagation condition to communicate with the user. At the same time, the UAV also needs to balance the communication link with the BS. 

## Geographical data
We use the geographical data captured in central Washington DC, USA in 2013 to study the communication environment. The original geographical dataset is available at USGS database: http://ngmdb.usgs.gov.

## Data pre-processing
Based on the satellite image, street map, and height map in the original dataset, we pre-process the data and extract a 800 meter by 800 meter area of interest for the study. The data is stored in a MATLAB MAT-file 'urbanMap_DC'mat' with the following structure:

          BldArea: {1265×4 cell}
             Maps: {2×1 cell}
             Xvec: [1×834 double]
             Yvec: [1×834 double]
    meterPerPixel: 0.9594
         stepsize: 0.9594
               x0: [0 0]

- 'Maps' contains two maps: Maps{1} is a 834x834 matrix for the vegetation map, each entry recording the height of the vegetation at the corresponding pixel. Each pixel spans meterPerPixel = 0.9594 meter. Maps{2} is a 834x834 matrix for the building map, each entry recording the height of the building. 
- 'Xvec' and 'Yvec' respectively record the x and y coordinates of the pixels. 
- 'BldArea' contains the coordinates of the contours of 1265 buildings in the area of interest.


## Ray-tracing statistics and user deployment
We use a direct ray-tracing method to analyze the statistics of the air-to-ground propagation under the 'urbanMap_DC' dataset described above. Specifically, if a propagation path is blocked by vegetation, it is classified as obstructed-LOS (OLOS); if it is blocked by building structure, it is classified as NLOS; if it is not blocked by any structure or vegetation, it is classified as LOS.

In addition, we deploy and record 10,000 user locations in the street area in 'urbanMap_DC' for the reproduction of our experiments. The information is stored in 'topology_DC.mat' with the following structure:

      Angles: [1×89 double]
     LosFreq: [3×89 double]
    PosArray: [10000×2 double]

- Angles: a vector of elevation angle in radius
- LosFreq: the relative frequency (probability) of {LOS, OLOS, NLOS} under the corresponding elevation angle. The three rows of LosFreq represent LOS, OLOS, and NLOS, respectively. For example, LosFreq(1, i) represents the empirical probability density (or relative frequency) of LOS under elevation angle Angles(i).
- PosArray: the xy coordinates of the 10,000 user locations.

## Algorithm 
(in progress)

Contact information:
- Author: Junting CHEN
- Email: juntingc@cuhk.edu.cn
- Date: Nov 1, 2020