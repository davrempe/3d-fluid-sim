//----------------------------------------------------------------------
// Title: 2D Fluid Simulation
// Author: Davis Rempe
//
// Implementation based on "Fluid Simulation for Computer Graphics"
// by Robert Bridson.
//----------------------------------------------------------------------

#include <string>
#include <fstream>
#include <iostream>
#include <cmath>

#include "FluidSolver3d.h"
#include "FluidRenderer3d.h"
#include "SimUtil.h"

//----------------------------------------------------------------------
// Execution Options
//----------------------------------------------------------------------

// whether to run the simulation
const bool RUN_SIM = true;
// whether to run rendering
const bool RUN_RENDERING = true;

//----------------------------------------------------------------------
// Simulation Parameters
//----------------------------------------------------------------------

// resolution of the grid to use (width, height, depth) = (x, y, z)
const int GRID_WIDTH = 100;
const int GRID_HEIGHT = 50;
const int GRID_DEPTH = 50;
// grid cell width (in meters)
const float GRID_CELL_WIDTH = 0.005f;
// simulation time step (in seconds)
const float TIME_STEP = 0.01f;

//----------------------------------------------------------------------
// I/O Parameters
//----------------------------------------------------------------------

// input file for initial system state - an .obj file created in Maya with "solid" and "fluid" group mesh defined
const std::string INITIAL_GEOMETRY_FILE_IN = "obj_test.obj";
// output file for particle data
const std::string PARTICLE_DATA_FILE_OUT = "particle_out.csv";
// the number of frames to simulate
const int NUM_SIM_FRAMES = 100;
// frame rate for render (fps)
const float FRAME_RATE = 25.0f;
// time step between outputted frames
const float FRAME_TIME_STEP = 1.0f / FRAME_RATE;

//----------------------------------------------------------------------
// Global Variables
//----------------------------------------------------------------------



int main(int argc, char** argv) {
	if (RUN_SIM) {
		// open and clear output file
		std::ofstream *particleOut = new std::ofstream(PARTICLE_DATA_FILE_OUT, std::ofstream::trunc);

		FluidSolver3D solver(GRID_WIDTH, GRID_HEIGHT, GRID_CELL_WIDTH, TIME_STEP);
		std::cout << "Simulating Frame 1" << std::endl;
		solver.init(INITIAL_GEOMETRY_FILE_IN);
		
		// run simulation
		// save initial frame
		solver.saveParticleData(particleOut);
		float t = 0.0f;
		for (int framesOut = 1; framesOut < NUM_SIM_FRAMES; framesOut++) {
			std::cout << "Simulating Frame " << framesOut + 1 << std::endl;
			
			float subTime = 0.0f;
			while (subTime < FRAME_TIME_STEP) {
				// perform sim time step
				solver.step();

				subTime += TIME_STEP;
				t += TIME_STEP;
			}
			// render at current time
			solver.saveParticleData(particleOut);
		}

		// cleanup
		particleOut->close();
		delete particleOut;
	}

	if (RUN_RENDERING) {
		FluidRenderer3D renderer(INITIAL_GEOMETRY_FILE_IN, PARTICLE_DATA_FILE_OUT, FRAME_RATE, GRID_WIDTH, GRID_HEIGHT, GRID_CELL_WIDTH);
		renderer.init(argc, argv);
		renderer.render();
	}


	return 0;
}