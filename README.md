# compmodel-

1. Introduction

This project implements to learn the compositional model from edges. 
The compositional model assumes that all objects are composed by parts, and parts are composed by subparts, and so on. The lower level parts can be shared across all objects, while the higher level larger parts are specific to certain object.

The basic idea is to use Gabor filters to find edges (or points) from all images, and each edge can be thought as the primary part. For the current layer, for each part, we look around its neighbors, and combine them together into larger parts. Then, we select the most frequent parts into the next layer. We repeat the two steps (combination and selection) until we reach the defined number of layers. 


2. Functionality 

complearn.cpp, main file for part learning

compinference.cpp, main file for part inference (bottom up)

comptopdown.cpp, main file for top down inference

mexc_complearn2.cpp, matlab mex function, interface of complearn, main file for learning

mexc_compinference.cpp,  matlab mex function, interface of compinference, main file for inference (bottom up)

mexc_compinference_batch.cpp, matlab mex function, with batch processing

mexc_comptopdown.cpp,  matlab mex function, interface of comptopdown, main file for top down inference

filter.cpp, main file for convolution

chechparts.cpp, main file for part selection

util2.cpp, main file for utility functions, such as sorting, filtering, etc.


3. References:

[1] Joachim Utans. Learning in Compositional Hierarchies: Inducing the Structure of Objects from Data, In NIPS, 1993.  
[2] Sanja Fidler, Ales Leonardis. Towards scalable representations of object categories: Learning a hierarchy of parts, In CVPR,2007.
[3] Gang Chen, Sanja Fidler and Jason J. Corso. Grammar Structure Leearning for Object Detection, 2014.

