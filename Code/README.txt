* Please use a version of GCC that supports OpenMP 4.5, such as GCC-7.4.0, which can be loaded on Hamilton using: module load gcc/7.4.0. This is to ensure that the array-reduction used in Step 6 compiles.

* To compile, please use the following command:
	g++ [fileName] -fopenmp -o [executableName] --std=c++11 -O3

* During my experimentation phase I implemented two different versions of the parallelisation in Step 6, hence why there are two versions of the file. The original version (solution-step6.c) is the version that uses atomic updates, whereas the second version (solution-step6-1.c) is the version that uses an array reduction. If you only accept one version, please mark the original version. Both are discussed in Step 7.