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
	double sr, sr3, sr6, sr12, dx, dy, dz, r, fx, fy, fz, lj;
	int i = 0, j = 0, k = 0;
	double force[NumberOfBodies][3];

	maxV   = 0.0;
	minDx  = std::numeric_limits<double>::max();

#pragma omp parallel
	{

#pragma omp for collapse(2) private(i, k)
		for(i = 0; i < NumberOfBodies; i++){
			for(k = 0; k < 3; k++){
				force[i][k] = 0.0;
			}
		}


		//In this version of the parallelisation, I moved the variable declarations out of the parallel region and thus had to explicitly
		//privatise them in this clause. I also changed the way each thread access the force array - each thread will now receive an individual
		//copy of the array (this behaviour is implicit in the reduction), and upon exiting the loop the array will be reduced component-wise
		//using the addition operator. This removes the need for atomic updates, but may significantly increase parallelisation overhead
#pragma omp for collapse(2) private(i, j, sr, sr3, sr6, sr12, dx, dy, dz, r, fx, fy, fz, lj) reduction(min:minDx) reduction(+:force[0:NumberOfBodies * 3])
		for(i = 0; i < NumberOfBodies; i++){ 
			for(j = 0; j < NumberOfBodies; j++){
				
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

				//Atomic updates are now no longer required
				force[i][0] = fx;
				force[i][1] = fy;
				force[i][2] = fz;

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