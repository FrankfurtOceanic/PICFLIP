#ifndef GRID_SOLVER_H
#define GRID_SOLVER_H

#define BIG_GRID

#define GLEW_STATIC

#include <vector>
#include <deque>
#include <memory>
#include <GL/glew.h>


#include "Vertex.h"

#ifndef BIG_GRID
const int cols = 3, rows = 3;
const int cellSize = 60; // cell size in display
#else
const int cols = 40, rows = 30;
const int cellSize = 30; // cell size in display
#endif

//Particle variables
const float gravity = -9.81f; //0.0f;
const int pSystem = 3;
const int particlesX = 100;
const int particlesY = 100;
const int numParticles = particlesX * particlesY; //number of particles

const bool usingLedge = false;

const int frameWidth = cols * cellSize, frameHeight = rows * cellSize;

const float penetrationCoef = 0.5f; //penetration of boundary cells

extern float nv, kappa, buoyancy, dt, velScale, PICratio;
extern float simTime;
extern float dx;
extern float refinementThreshold;
extern int iterations;

extern double tpSrc[rows + 2][cols + 2];
extern double vxSrc[rows + 2][cols + 2];
extern double vySrc[rows + 2][cols + 2];


extern double (*cur_tp)[cols + 2];

extern std::deque<std::shared_ptr<std::vector<Vertex>>> filaments;
extern std::deque<double> ages;
extern float maxAge;


extern Vertex gridVertices[rows + 1][cols + 1];
extern GLuint gridIndices[6 * rows * cols];
extern Vertex velVertices[rows][cols][2];

extern Vertex particleVertices[numParticles];

void initGrid();
double getReferenceTemperature(double tp[rows + 2][cols + 2]);
void step();
void updateGrid();
void addFilament(double x, double y);
void updateFilament();

void displayParticles(int winW, int winH);

#endif
