#ifndef TRIANGLE_H
#define TRIANGLE_H
#include "Vec3.h"
#include "Ray.h"
#include "Plane.h"
#include <optional>

struct RayTriangleIntersection{
    bool intersectionExists;
    float t;
    float w0,w1,w2;
    unsigned int tIndex;
    Vec3 intersection;
    Vec3 normal;
};

class Triangle {
private:
    Vec3 m_c[3] , m_normal;
    float area;
public:
    Triangle() {}
    Triangle( Vec3 const & c0 , Vec3 const & c1 , Vec3 const & c2 ) {
        m_c[0] = c0;
        m_c[1] = c1;
        m_c[2] = c2;
        updateAreaAndNormal();
    }
    void updateAreaAndNormal() {
        Vec3 nNotNormalized = Vec3::cross( m_c[1] - m_c[0] , m_c[2] - m_c[0] );
        float norm = nNotNormalized.length();
        m_normal = nNotNormalized / norm;
        area = norm / 2.f;
    }
    void setC0(Vec3 const & c0) { m_c[0] = c0; } // remember to update the area and normal afterwards!
    void setC1(Vec3 const & c1) { m_c[1] = c1; } // remember to update the area and normal afterwards!
    void setC2(Vec3 const & c2) { m_c[2] = c2; } // remember to update the area and normal afterwards!
    Vec3 const & normal() const { return m_normal; }
    Vec3 projectOnSupportPlane( Vec3 const & p ) const {
        Vec3 result;
        //TODO completer
        return result;
    }
    float squareDistanceToSupportPlane( Vec3 const & p ) const {        
        float result;
        //TODO completer
        return result;
    }
    float distanceToSupportPlane(Vec3 const & p) const { return sqrt( squareDistanceToSupportPlane(p) ); }

    bool isParallelTo( Line const & L ) const {
        Plane supportPlane(m_c[0],m_normal); // On prend n'importe quel point du triangle et la normale au triangle pour définir le plan support
        return supportPlane.isParallelTo(L);
    }

    Vec3 getIntersectionPointWithSupportPlane(Line const & L) const {
        // you should check first that the line is not parallel to the plane!

        Vec3 result = Vec3(-100.0f,-100.0f,-100.0f);
        
        if (isParallelTo(L)) return result;

        Vec3 P0 = L.origin();
        Vec3 P1 = m_c[0];
        Vec3 D = L.direction();
        float t;

        t = Vec3::dot(P1 - P0,m_normal) / Vec3::dot(D,m_normal);

        result = P0 + t * D;

        return result;
        
    }

    void computeBarycentricCoordinates(Vec3 const &p, float &u0, float &u1, float &u2) const {

        Vec3 A = m_c[0];
        Vec3 B = m_c[1];
        Vec3 C = m_c[2];

        Vec3 v0 = B - A;
        Vec3 v1 = C - A;
        Vec3 v2 = p - A;
        float d00 = Vec3::dot(v0, v0);
        float d01 = Vec3::dot(v0, v1);
        float d11 = Vec3::dot(v1, v1);
        float d20 = Vec3::dot(v2, v0);
        float d21 = Vec3::dot(v2, v1);
        float denom = d00 * d11 - d01 * d01;

        u1 = (d11 * d20 - d01 * d21) / denom;
        u2 = (d00 * d21 - d01 * d20) / denom;
        u0 = 1.0f - u1 - u2;

    }

    // void computeBarycentricCoordinates(Vec3 const &p , float &u0 , float &u1 , float &u2) const {
    //     //TODO Complete
    //     Vec3 A = m_c[0];
    //     Vec3 B = m_c[1];
    //     Vec3 C = m_c[2];

    //     float areaABC = Triangle(A,B,C).area;
    //     float areaPCB = Triangle(p,C,B).area;
    //     float areaPCA = Triangle(p,C,A).area;
    //     float areaPBA = Triangle(p,B,A).area;

    //     u0 = areaPCB/areaABC;
    //     u1 = areaPCA/areaABC;
    //     u2 = areaPBA/areaABC;
    
    // }

    RayTriangleIntersection getIntersection( Ray const & ray ) const {

        RayTriangleIntersection result;
        result.intersectionExists = false;
        // 1) check that the ray is not parallel to the triangle: // Déjà fait dans getIntersectionPointWithSupportPlane

        // if(isParallelTo(ray)) return result;

        // 2) check that the triangle is "in front of" the ray:

        Vec3 intersectionPoint = getIntersectionPointWithSupportPlane(ray);
        Vec3 originToIntersection = intersectionPoint - ray.origin();
        if(Vec3::dot(ray.direction(),originToIntersection) < 0) return result;

        // 3) check that the intersection point is inside the triangle:
        // CONVENTION: compute u,v such that p = w0*c0 + w1*c1 + w2*c2, check that 0 <= w0,w1,w2 <= 1

        float u0,u1,u2;
        computeBarycentricCoordinates(intersectionPoint, u0, u1, u2);
        if(u0 < 0 || u0 > 1 || u1 < 0 || u1 > 1 || u2 < 0 || u2 > 1) return result;

        // 4) Finally, if all conditions were met, then there is an intersection! :

        //std::cout<<intersectionPoint<<std::endl;

        result.intersectionExists = true;
        result.intersection = intersectionPoint;
        result.w0 = u0;
        result.w1 = u1;
        result.w2 = u2;
        result.t = originToIntersection.norm();
        Vec3 normal = Vec3::cross(m_c[1] - m_c[0], m_c[2] - m_c[0]);
        normal.normalize();
        result.normal = normal;

        return result;
    }
};
#endif
