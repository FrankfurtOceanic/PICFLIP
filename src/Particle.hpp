#pragma once
#define _USE_MATH_DEFINES
#include <cmath>

#include <vector>

#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/fwd.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

class Particle {
public:
    /** Identifies this particles position in the particle list */
    int index;


    glm::vec3 color{ 0.1f, 0.1f, 1.0f };

    float size = 10;

    float mass = 1;

    glm::vec2 p;
    glm::vec2 v;
    glm::vec2 p0; //initial position
    glm::vec2 v0; //initial velocity
    glm::vec2 f;



    /** Default constructor */
    Particle() {}

    /**
     * Creates a particle with the given position and velocity
     * @param x
     * @param y
     * @param vx
     * @param vy
     */
    Particle(float x, float y, float vx, float vy) {
        p0 = glm::vec2(x, y);
        v0 = glm::vec2(vx, vy);
        reset();
    }

    /**
     * Resets the position of this particle
     */
    void reset() {
        p = p0;
        v = v0;
        f = glm::vec2(0, 0);
    }



    /**
     * Computes the distance of a point to this particle
     * @param x
     * @param y
     * @return the distance
     */
    float distance(float x, float y) {
        glm::vec2 diff = p - glm::vec2(x, y);

        return (float)sqrt(diff.x * diff.x + diff.y * diff.y);
    }
};
