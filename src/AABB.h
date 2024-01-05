#include "Vec3.h"
#include "Ray.h"
#include "Mesh.h"
#include <GL/gl.h>
#include <cfloat>
#include <vector>
#include <algorithm>

class AABB {
    public:

        AABB() = default;

        AABB(const Vec3 vmin, const Vec3 vmax) {
            min = vmin;
            max = vmax;
        }

        ~AABB() = default;

        // Version de base
        // // Intersection avec une AABB
        // bool intersect(const Ray &ray) {
        //     Vec3 origin = ray.origin();
        //     Vec3 direction = ray.direction();
        //     direction.normalize();

        //     float txmin = (min[0] - origin[0]) / direction[0];
        //     float txmax = (max[0] - origin[0]) / direction[0];

        //     if(txmin > txmax) std::swap(txmin,txmax);

        //     float tymin = (min[1] - origin[1]) / direction[1];
        //     float tymax = (max[1] - origin[1]) / direction[1];

        //     if(tymin > tymax) std::swap(tymin,tymax);

        //     float tzmin = (min[2] - origin[2]) / direction[2];
        //     float tzmax = (max[2] - origin[2]) / direction[2];

        //     if(tzmin > tzmax) std::swap(tzmin,tzmax);

        //     if(txmin > tzmax || tzmin > txmax) return false;

        //     if(tzmin > txmin) txmin = tzmin;

        //     if(tzmax < txmax) txmax = tzmax;

        //     return true;
        // }

        // Version optimisée
        // Intersection avec une AABB
        bool intersect(const Ray &ray) {
            Vec3 origin = ray.origin();
            Vec3 direction = ray.direction();
            Vec3 invDir;

            for (int i = 0; i < 3; ++i) {
                if (direction[i] != 0.0f) {
                    invDir[i] = 1.0f / direction[i];
                } else {
                    invDir[i] = direction[i] > 0 ? FLT_MAX : -FLT_MAX;
                }
            }

            float t[6]; // txmin, txmax, tymin, tymax, tzmin, tzmax

            for (int i = 0; i < 3; ++i) {
                if (std::abs(invDir[i]) < FLT_EPSILON) {
                    // Le rayon est parallèle à la face de l'AABB sur cet axe
                    t[i * 2] = (min[i] - origin[i]) > 0 ? -FLT_MAX : FLT_MAX;
                    t[i * 2 + 1] = (max[i] - origin[i]) > 0 ? FLT_MAX : -FLT_MAX;
                } else {
                    t[i * 2] = (min[i] - origin[i]) * invDir[i];
                    t[i * 2 + 1] = (max[i] - origin[i]) * invDir[i];
                }

                if (t[i * 2] > t[i * 2 + 1]) std::swap(t[i * 2], t[i * 2 + 1]);
            }

            float tmin = std::max({t[0], t[2], t[4]});
            float tmax = std::min({t[1], t[3], t[5]});

            return tmax >= tmin && tmax >= 0.0f;
        }

        static AABB computeGlobalBox(const std::vector<Triangle>& triangles) {
            Vec3 min(FLT_MAX, FLT_MAX, FLT_MAX);
            Vec3 max(-FLT_MAX, -FLT_MAX, -FLT_MAX);

            for (Triangle tri : triangles) {
                for (Vec3& vertex : tri.getVertices()) {
                    min = Vec3::min(min, vertex);
                    max = Vec3::max(max, vertex);
                }
            }

            return AABB(min, max);
        }

        static AABB combine(const AABB& a, const AABB& b) {
            Vec3 newMin(std::min(a.min[0], b.min[0]), std::min(a.min[1], b.min[1]), std::min(a.min[2], b.min[2]));
            Vec3 newMax(std::max(a.max[0], b.max[0]), std::max(a.max[1], b.max[1]), std::max(a.max[2], b.max[2]));
            return AABB(newMin, newMax);
        }


        void drawBox() {
            Vec3 corners[8];
            corners[0] = Vec3(min[0], min[1], min[2]);
            corners[1] = Vec3(max[0], min[1], min[2]);
            corners[2] = Vec3(min[0], max[1], min[2]);
            corners[3] = Vec3(max[0], max[1], min[2]);
            corners[4] = Vec3(min[0], min[1], max[2]);
            corners[5] = Vec3(max[0], min[1], max[2]);
            corners[6] = Vec3(min[0], max[1], max[2]);
            corners[7] = Vec3(max[0], max[1], max[2]);

            
            glLineWidth(1.0f);
            glBegin(GL_LINES);

            // Bas de la boîte
            glVertex3f(corners[0][0], corners[0][1], corners[0][2]);
            glVertex3f(corners[1][0], corners[1][1], corners[1][2]);

            glVertex3f(corners[1][0], corners[1][1], corners[1][2]);
            glVertex3f(corners[3][0], corners[3][1], corners[3][2]);

            glVertex3f(corners[3][0], corners[3][1], corners[3][2]);
            glVertex3f(corners[2][0], corners[2][1], corners[2][2]);

            glVertex3f(corners[2][0], corners[2][1], corners[2][2]);
            glVertex3f(corners[0][0], corners[0][1], corners[0][2]);

            // Haut de la boîte
            glVertex3f(corners[4][0], corners[4][1], corners[4][2]);
            glVertex3f(corners[5][0], corners[5][1], corners[5][2]);

            glVertex3f(corners[5][0], corners[5][1], corners[5][2]);
            glVertex3f(corners[7][0], corners[7][1], corners[7][2]);

            glVertex3f(corners[7][0], corners[7][1], corners[7][2]);
            glVertex3f(corners[6][0], corners[6][1], corners[6][2]);

            glVertex3f(corners[6][0], corners[6][1], corners[6][2]);
            glVertex3f(corners[4][0], corners[4][1], corners[4][2]);

            // Côtés de la boîte
            glVertex3f(corners[0][0], corners[0][1], corners[0][2]);
            glVertex3f(corners[4][0], corners[4][1], corners[4][2]);

            glVertex3f(corners[1][0], corners[1][1], corners[1][2]);
            glVertex3f(corners[5][0], corners[5][1], corners[5][2]);

            glVertex3f(corners[2][0], corners[2][1], corners[2][2]);
            glVertex3f(corners[6][0], corners[6][1], corners[6][2]);

            glVertex3f(corners[3][0], corners[3][1], corners[3][2]);
            glVertex3f(corners[7][0], corners[7][1], corners[7][2]);

            glEnd();
    }

    private:
        Vec3 min, max;
};