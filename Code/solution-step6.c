// Translate this file with
//
// g++ -O3 --std=c++11 assignment-2018.c -o assignment
//
// Run it with
//
// ./assignment
//
// There should be a result.pvd file that you can open with Paraview.
// Sometimes, Paraview requires to select the representation "Point Gaussian"
// to see something meaningful.
//
// (C) 2018 Tobias Weinzierl

#include <cmath>
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <math.h>
#include <limits>
#include <omp.h>


double t          = 0;
double tFinal     = 0;
double tPlot      = 0;
double tPlotDelta = 0;

int NumberOfBodies = 0;

/**
 * Pointer to pointers. Each pointer in turn points to three coordinates, i.e.
 * each pointer represents one molecule/particle/body.
 */
double** x;

/**
 * Equivalent to x storing the velocities.
 */
double** v;

/**
 * One mass entry per molecule/particle.
 */
double*  mass;

/**
 * Global time step size used.
 */
double   timeStepSize = 0.1;

/**
 * Maximum velocity of all particles.
 */
double   maxV;

/**
 * Minimum distance between two elements.
 */
double   minDx;


/**
 * Set up scenario from the command line.
 *
 * This operation is not to be changed in the assignment.
 */
void setUp(int argc, char** argv) {
  NumberOfBodies = (argc-2) / 7;

  x    = new double*[NumberOfBodies];
  v    = new double*[NumberOfBodies];
  mass = new double [NumberOfBodies];

  int readArgument = 1;

  tPlotDelta  = std::stof(argv[readArgument]); readArgument++;
  tFinal      = std::stof(argv[readArgument]); readArgument++;

  for (int i=0; i<NumberOfBodies; i++) {
    x[i] = new double[3];
    v[i] = new double[3];

    x[i][0] = std::stof(argv[readArgument]); readArgument++;
    x[i][1] = std::stof(argv[readArgument]); readArgument++;
    x[i][2] = std::stof(argv[readArgument]); readArgument++;

    v[i][0] = std::stof(argv[readArgument]); readArgument++;
    v[i][1] = std::stof(argv[readArgument]); readArgument++;
    v[i][2] = std::stof(argv[readArgument]); readArgument++;

    mass[i] = std::stof(argv[readArgument]); readArgument++;

    if (mass[i]<=0.0 ) {
      std::cerr << "invalid mass for body " << i << std::endl;
      exit(-2);
    }
  }

  std::cout << "created setup with " << NumberOfBodies << " bodies" << std::endl;
  
  if (tPlotDelta<=0.0) {
    std::cout << "plotting switched off" << std::endl;
    tPlot = tFinal + 1.0;
  }
  else {
    std::cout << "plot initial setup plus every " << tPlotDelta << " time units" << std::endl;
    tPlot = 0.0;
  }
}


std::ofstream videoFile;


/**
 * This operation is not to be changed in the assignment.
 */
void openParaviewVideoFile() {
	videoFile.open( "result.pvd" );
	videoFile << "<?xml version=\"1.0\"?>" << std::endl
		<< "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">" << std::endl
		<< "<Collection>";
}


/**
 * This operation is not to be changed in the assignment.
 */
void closeParaviewVideoFile() {
	videoFile << "</Collection>"
		<< "</VTKFile>" << std::endl;
}


/**
 * The file format is documented at http://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf
 *
 * This operation is not to be changed in the assignment.
 */
void printParaviewSnapshot() {
	static int counter = -1;
	counter++;
	std::stringstream filename;
	filename << "result-" << counter <<  ".vtp";
	std::ofstream out( filename.str().c_str() );
	out << "<VTKFile type=\"PolyData\" >" << std::endl
		<< "<PolyData>" << std::endl
		<< " <Piece NumberOfPoints=\"" << NumberOfBodies << "\">" << std::endl
		<< "  <Points>" << std::endl
		<< "   <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">";

	for (int i=0; i<NumberOfBodies; i++) {
		out << x[i][0]
			<< " "
			<< x[i][1]
			<< " "
			<< x[i][2]
			<< " ";
	}

	out << "   </DataArray>" << std::endl
		<< "  </Points>" << std::endl
		<< " </Piece>" << std::endl
		<< "</PolyData>" << std::endl
		<< "</VTKFile>"  << std::endl;

	videoFile << "<DataSet timestep=\"" << counter << "\" group=\"\" part=\"0\" file=\"" << filename.str() << "\"/>" << std::endl;
}


/**
 * This is the only operation you are allowed to change in the assignment.
 */
void updateBody() {
	const double epsilon4 = 1.65 * pow(10, -21) * -4.0;
	const double sigma = 3.4 * pow(10, -10);

	maxV   = 0.0;
	minDx  = std::numeric_limits<double>::max();
	int i = 0, j = 0, k = 0;
	double force[NumberOfBodies][3];

#pragma omp parallel
	{
		//Create our private variables inside the omp parallel region - this ensures they are private to each thread
		double sr, sr3, sr6, sr12, dx, dy, dz, r, fx, fy, fz, lj;

#pragma omp for collapse(2) private(i, k)
		//Collapse the nested loop into one loop space, so that each i,k iteration is distributed evenly among the threads
		//Throughout, I explicitly privatise the counter variables even though it is done so internally by OMP, for semantic 
		//clarity
		for(i = 0; i < NumberOfBodies; i++){
			for(k = 0; k < 3; k++){
				force[i][k] = 0.0;
			}
		}

#pragma omp for collapse(2) private(i, j) reduction(min:minDx) 
		//To ease the parallelisation process, I 'de-optimised' this loop - i.e, instead of iterating through every combination
		//of the particles and using Newton's 3rd Law to decrease the number of computation requred, I iterate through every
		//permutation of the particles. 
		for(i = 0; i < NumberOfBodies; i++){ 
			for(j = 0; j < NumberOfBodies; j++){

				//All of the following computations are done privately per thread, as all of these variable were declared inside
				//the parallel region
				dx = x[i][0] - x[j][0];
				dy = x[i][1] - x[j][1];
				dz = x[i][2] - x[j][2];

				r = sqrt((dx * dx) + (dy * dy) + (dz * dz));

				sr = sigma / r;
				sr3 = sr * sr * sr;
				sr6 = sr3 * sr3;
				sr12 = sr6 * sr6;

				lj = epsilon4 * ((-12.0 * sr12) + (6.0 * sr6)) / r;

				fx = (abs(dx) < 1e-17) ? 0 : dx * lj / r;
				fy = (abs(dy) < 1e-17) ? 0 : dy * lj / r;
				fz = (abs(dz) < 1e-17) ? 0 : dz * lj / r;

				//The only shared variable is the force array, which without the atomic data accesses would lead to race conditions
				//By default, #pragma omp atomic is equivalent to #pragma omp atomic update, which is the functionality we need
				#pragma omp atomic
				force[i][0] += fx;
				#pragma omp atomic
				force[i][1] += fy;
				#pragma omp atomic
				force[i][2] += fz;

				//Each thread maintains its own version of minDx, which it uses across all of the iterations it computes. At the
				//end, these private minDx values are reduced using the 'min' function due to the reduction clause in the #pragma omp for
				minDx = std::min(minDx, r);
			}		
		}

#pragma omp for collapse(2) private(i, k) reduction(max:maxV)
		for(i = 0; i < NumberOfBodies; i++){ 
			for(k = 0; k < 3; k++){
				x[i][k] = timeStepSize * v[i][k] + x[i][k];
				v[i][k] = v[i][k] + timeStepSize * force[i][k] / mass[i];

				maxV = std::max(maxV, std::sqrt( v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]));
			}

		}
	}

	t += timeStepSize;
}


/**
 * Main routine.
 *
 * Not to be changed in assignment.
 */
int main(int argc, char** argv) {
  if (argc==1) {
    std::cerr << "usage: " + std::string(argv[0]) + " snapshot final-time objects" << std::endl
              << "  snapshot        interval after how many time units to plot. Use 0 to switch off plotting" << std::endl
              << "  final-time      simulated time (greater 0)" << std::endl
              << std::endl
              << "Examples:" << std::endl
              << "0.01  100.0   0 0 0 1.0   0   0 1.0 \t One body moving form the coordinate system's centre along x axis with speed 1" << std::endl
              << "0.01  100.0   0 0 0 1.0   0   0 1.0 0 1.0 0 1.0 0   0 1.0 \t One spiralling around the other one" << std::endl
              << "0.01  100.0 3.0 0 0   0 1.0   0 0.4 0   0 0   0 0   0 0.2 2.0 0 0 0 0 0 1.0 \t Three body setup from first lecture" << std::endl
              << std::endl
              << "In this naive code, only the first body moves" << std::endl;

    return -1;
  }
  else if ( (argc-3)%7!=0 ) {
    std::cerr << "error in arguments: each planet is given by seven entries (position, velocity, mass)" << std::endl;
    return -2;
  }

  setUp(argc,argv);

  openParaviewVideoFile();

  int snapshotCounter = 0;
  if (t > tPlot) {
    printParaviewSnapshot();
    std::cout << "plotted initial setup" << std::endl;
    tPlot = tPlotDelta;
  }

  int timeStepCounter = 0;
  while (t<=tFinal) {
    updateBody();
    timeStepCounter++;
    if (t >= tPlot) {
      printParaviewSnapshot();
      std::cout << "plot next snapshot"
    		    << ",\t time step=" << timeStepCounter
    		    << ",\t t="         << t
				<< ",\t dt="        << timeStepSize
				<< ",\t v_max="     << maxV
				<< ",\t dx_min="    << minDx
				<< std::endl;

      tPlot += tPlotDelta;
    }
  }

  closeParaviewVideoFile();

  return 0;
}
