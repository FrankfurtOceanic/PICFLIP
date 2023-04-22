#include <vector>

#define GLEW_STATIC
#include <GL/glew.h>
#include <GLFW/glfw3.h>

#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "GLSL.h"

#include "Particle.hpp"
#include "GridSolver.h"

#include <Eigen/Dense>
using Eigen::MatrixXf;
using Eigen::VectorXf;

using namespace std;
/**
 * Implementation of a simple particle system
 * @author kry
 */
class ParticleSystem {

public:
    std::vector<Particle*> particles;

    /**
     * Creates an empty particle system
     */
    ParticleSystem() {
        // do nothing
    }

    void simulateParticles(glm::vec2 accel, float dt) {
        //symplectic euler updating the velocity first then position
        for (Particle* p : particles) {
            
            p->v += accel * dt;
            
            //p->v.x += force.x * p->mass * dt;
            //p->v.y += force.y * p->mass * dt;

            p->p += p->v * dt;

            //fixing positions in case of boundary penetration
            //move it outside by a fraction of the cell size
            if (p->p.x < 0 + penetrationCoef * dx) {
                p->p.x = penetrationCoef * dx;
            }
            if (p->p.x > cols * dx - penetrationCoef * dx) {
                p->p.x = dx * (cols - penetrationCoef);
            }

            if (p->p.y < 0 + penetrationCoef * dx) {
                p->p.y = penetrationCoef * dx;
            }
            if (p->p.y > rows * dx - penetrationCoef * dx) {
                p->p.y = dx * (rows - penetrationCoef);
            }
        }
    }


    /**
     * Creates one of a number of simple test systems.
     * @param which
     */
    void createSystem(int which) {
        
        //grid of particles system
        if (which == 1) {

            //creates a grid of particles (particlesX * particleY) with spacing (space)
            //int particlesX = 10;
            //int particlesY = 10; 
            float spacing = dx / 6;

            float offsetX = 0;
            float offsetY = rows * dx - dx/2;

            for (int i = 0; i < particlesX; i++) {
                for (int j = 0; j < particlesY; j++) {
                    Particle* p1 = new Particle(offsetX + i * spacing, offsetY - j * spacing, 0, 0);
                    p1->index = particles.size();
                    particles.push_back(p1);
                }
            }
        }

        //single particle system
        if (which == 2) {
            float offsetX = dx * cols / 2;
            float offsetY = dx * rows / 2;
            Particle* p1 = new Particle(offsetX, offsetY, 0, 0);
            p1->index = particles.size();
            particles.push_back(p1);

        }

        if (which == 3) {

            //upper fluid

            int count = 0.01 * particlesX * particlesY;
            int rootCount = sqrt(count);
            float spacing = dx / 4;

            float offsetX = cols * dx/2 - (rootCount-1) * spacing / 2;
            float offsetY = (rows - 1) * dx ;

            for (int i = 0; i < rootCount; i++) {
                for (int j = 0; j < rootCount; j++) {
                    Particle* p1 = new Particle(offsetX + i * spacing, offsetY - j * spacing, 0, 0);
                    p1->index = particles.size();
                    particles.push_back(p1);
                }
            }


           
            offsetX = 0;
            offsetY = rows * dx * 0.25;

            count = particlesX * particlesY - count;

            int countX = sqrt(count * (4 * cols) / rows);
            int countY = float(rows) / float(cols * 4) * countX;

            float spacingX = (cols * dx)/ (countX-1);
            float spacingY =  offsetY/ (countY - 1);

            for (int i = 0; i < countX; i++) {
                for (int j = 0; j < countY; j++) {
                    Particle* p1 = new Particle(offsetX + i * spacingX, offsetY - j * spacingY, 0, 0);
                    p1->index = particles.size();
                    particles.push_back(p1);
                }
            }
            //resting fluid


        }
        
        /*
        if (which == 1) {
            glm::vec2 p(100, 100);
            glm::vec2 d(20, 0);
            Particle* p1 = new Particle(p.x - d.y, p.y + d.x, 0, 0);
            p1->index = particles.size();
            particles.push_back(p1);
            Particle* p2 = new Particle(p.x + d.y, p.y - d.x, 0, 0);
            p2->index = particles.size();
            particles.push_back(p2);
            springs.push_back(new Spring(p1, p2));
            p1->pinned = true;
            p2->pinned = true;
            p += d;
            p += d;
            int N = 10;
            for (int i = 1; i < N; i++) {
                //d.set( 20*Math.cos(i*Math.PI/N), 20*Math.sin(i*Math.PI/N) );                
                Particle* p3 = new Particle(p.x - d.y, p.y + d.x, 0, 0);
                p3->index = particles.size();
                particles.push_back(p3);
                Particle* p4 = new Particle(p.x + d.y, p.y - d.x, 0, 0);
                p4->index = particles.size();
                particles.push_back(p4);
                springs.push_back(new Spring(p3, p1));
                springs.push_back(new Spring(p3, p2));
                springs.push_back(new Spring(p3, p2));
                springs.push_back(new Spring(p4, p1));
                springs.push_back(new Spring(p4, p2));
                springs.push_back(new Spring(p4, p3));
                p1 = p3;
                p2 = p4;
                p += d;
                p += d;
            }
        }
        else if (which == 2) {
            Particle* p1 = new Particle(320, 100, 0, 0);
            p1->index = particles.size();
            particles.push_back(p1);
            Particle* p2 = new Particle(320, 200, 0, 0);
            p2->index = particles.size();
            particles.push_back(p2);
            p1->pinned = true;
            springs.push_back(new Spring(p1, p2));
        }
        else if (which == 3) {
            float ypos = 100;
            Particle* p0 = NULL;
            Particle* p1 = new Particle(320, ypos, 0, 0);
            Particle* p2;
            p1->index = particles.size();
            p1->pinned = true;
            particles.push_back(p1);
            int N = 10;
            for (int i = 0; i < N; i++) {
                ypos += 20;
                p2 = new Particle(320, ypos, 0, 0);
                p2->index = particles.size();
                particles.push_back(p2);
                springs.push_back(new Spring(p1, p2));
                // Hum.. this is not great in comparison to a proper bending energy...
                // use Maple to generate some code though, as it is painful to write by hand! :(
                if (p0 != NULL) springs.push_back(new Spring(p2, p0));
                p0 = p1;
                p1 = p2;
            }

        }
        else if (which == 4) {

            float radius = 80;
            float ypos = 100;
            float xpos = 100;
            float yvel = 500;
            float xvel = 500;
            int n = 5;

            float pi = 3.14159265358979323846;
            Particle* p0 = new Particle(xpos, ypos, xvel, yvel); //center
            p0->index = particles.size();
            particles.push_back(p0);
            Particle* prev = NULL;
            Particle* p = NULL;

            Particle* p1 = new Particle(xpos + radius, ypos, xvel, yvel);
            p1->index = particles.size();
            particles.push_back(p1);

            springs.push_back(new Spring(p1, p0));

            prev = p1;
            for (int i = 1; i < n; i++) {
                p = new Particle(xpos + radius * cos(2 * pi * i / n), ypos + radius * sin(2 * pi * i / n), xvel, yvel);
                p->index = particles.size();
                particles.push_back(p);

                springs.push_back(new Spring(p, p0));

                springs.push_back(new Spring(p, prev));
                prev = p;
            }
            springs.push_back(new Spring(p1, prev)); //add last spring between last particle and second particle


        }
        */
    }


    /**
     * Resets the positions of all particles
     */
    void resetParticles() {
        for (Particle* p : particles) {
            p->reset();
        }
        time = 0;
    }

    /**
     * Deletes all particles
     */
    void clearParticles() {
        for (Particle* p : particles) { delete p; }
        particles.clear();
    }

    /**
     * Gets the phase space state of the particle system
     * @param phaseSpaceState
     */
    
    /**
     * Fixes positions and velocities after a step to deal with collisions
     */
    void postStepFix() {
        // do wall collisions
        float r = restitution;
        for (Particle* p : particles) {
            if (p->p.x <= 0) {
                p->p.x = 0;
                if (p->v.x < 0) p->v.x = -p->v.x * r;
                if (p->f.x < 0) p->f.x = 0;
            }
            if (p->p.x >= width) {
                p->p.x = width;
                if (p->v.x > 0) p->v.x = -p->v.x * r;
                if (p->f.x > 0) p->f.x = 0;
            }

            if (p->p.y >= height) {
                p->p.y = height;
                if (p->v.y > 0) p->v.y = -p->v.y * r;
                if (p->f.y > 0) p->f.y = 0;
            }
            if (p->p.y <= 0) {
                p->p.y = 0;
                if (p->v.y < 0) p->v.y = -p->v.y * r;
                if (p->f.y < 0) p->f.y = 0;
            }
        }
    }

    /** Elapsed simulation time */
    double time = 0;


    VectorXf state;
    VectorXf stateOut;

    // these get created in init() and are probably useful for Backward Euler computations
    //ConjugateGradientMTJ CG;
    MatrixXf A;
    MatrixXf dfdx;
    MatrixXf dfdv;
    VectorXf deltaxdot;
    VectorXf b;
    VectorXf f;
    VectorXf xdot;

    /**
     * Fills in the provided vector with the particle velocities.
     * @param xd
     */


    /**
     * Sets the velocities of the particles given a vector
     * @param xd
     */


    /**
     *  Evaluates derivatives for ODE integration.
     * @param t time
     * @param p phase space state (don't modify)
     * @param dydt to be filled with the derivative
     */

    /*void derivs(float t, VectorXf& p, VectorXf& dpdt) {
        // set particle positions to given values
        setPhaseSpace(p);

        // TODO: Objective 2, for explicit integrators, compute forces, and accelerations, and set dpdt

        //Compute force: 

        //for each particle: first zero forces then apply gravity and viscous air damping
        for (Particle* p : particles) {
            p->clearForce();

            if (useGravity) {
                p->addForce(glm::vec2(0, p->mass * gravity)); //apply gravity if using
            }
            //apply viscous damping 

            p->addForce(-viscousDamping * p->v);

        }


        //filling out dpdt
        int count = 0;
        for (Particle* p : particles) {

            //calculate acceleration and velocity
            glm::vec2 acceleration = p->f / p->mass;      // acel = f/m
            glm::vec2 vel = p->v + float(t - time) * acceleration;   // v(1) = v0 + t * acel 

            dpdt[count++] = vel.x;
            dpdt[count++] = vel.y;
            dpdt[count++] = acceleration.x;
            dpdt[count++] = acceleration.y;
        }

    }
    */
    /** Time in seconds that was necessary to advance the system */
    float computeTime;

    /**
     * Advances the state of the system
     * @param elapsed
     */
    /*
    void advanceTime(float elapsed) {

        int n = getPhaseSpaceDim();

        double now = glfwGetTime();

        if (useExplicitIntegration) {
            if (n != state.size()) {
                state.resize(n);;
                stateOut.resize(n);
            }
            // TODO: See explicit stepping here
            getPhaseSpace(state);
            integrator->step(state, n, time, elapsed, stateOut, this);
            setPhaseSpace(stateOut);
        }
        else {
            if (f.size() != n) {
                init();
            }

            // TODO: Objective 8, your backward Euler implementation will go here!
            // Note that the init() method called above creates a bunch of very 
            // useful MTJ working variables for you, and the ConjugateGradientMTJ object.
            // Go look at that code now!




        }
        time = time + elapsed;
        postStepFix();
        computeTime = (glfwGetTime() - now);
    }
    */


    /**
     * Creates a new particle and adds it to the system
     * @param x
     * @param y
     * @param vx
     * @param vy
     * @return the new particle
     */
    Particle* createParticle(float x, float y, float vx, float vy) {
        Particle* p = new Particle(x, y, vx, vy);
        p->index = particles.size();
        particles.push_back(p);
        return p;
    }

    void remove(Particle* p) {

        particles.erase(std::remove(particles.begin(), particles.end(), p));
        // reset indices of each particle :(
        for (int i = 0; i < particles.size(); i++) {
            particles[i]->index = i;
        }
    }




    void init() {
        // do nothing
    }

    int height;
    int width;

    void display(int winW, int winH) {

        glPointSize(10);
        glBegin(GL_POINTS);
        for (Particle* p : particles) {
            double alpha = 0.5;
                //glColor4d( p->color.x, p->color.y, p->color.z, alpha );

                //custom coloring
                glm::vec2 uv = abs(glm::normalize(p->v));
                if (p->v.x == 0 && p->v.y == 0) {
                    glColor4d(1, 1, 1, 1); //color white if stationary
                }
                else {
                    glColor4d(uv.x, uv.y, 0, 1); //color uv if moving enough
                }
            glVertex2d(p->p.x + (winW - frameWidth)/2, p->p.y + (winH - frameHeight) / 2);
        }
        glEnd();


    }

    bool useGravity = true;
    float gravity = 9.8;
    float viscousDamping = 0.01; //known as c in the player
    /** should only go between 0 and 1 for bouncing off walls */
    float restitution = 0;
    int solverIterations = 100;
    bool useExplicitIntegration = true;
};
