#ifndef PLANE_H
#define PLANE_H
#include "Vec3.h"
#include "Line.h"
#include <cmath>
#define UNVALID_INTERSECTION Vec3(-M_PI,-M_PI,-M_PI)
class Plane {
private:
    Vec3 m_center , m_normal;
public:
    Plane() {}
    Plane( Vec3 const & c , Vec3 const & n ) {
        m_center = c;
        m_normal = n; m_normal.normalize();
    }

    void setCenter( Vec3 const & c ) { m_center = c; }

    void setNormal( Vec3 const & n ) { m_normal = n; m_normal.normalize(); }

    Vec3 const & center() const { return m_center; }

    Vec3 const & normal() const { return m_normal; }

    Vec3 project( Vec3 const & p ) const {
        Vec3 result;
        float distance = Vec3::dot(p - m_center,m_normal);
        result = p - distance * m_normal;
        return result;
    }

    float squareDistance( Vec3 const & p ) const { return (project(p) - p).squareLength(); }

    float distance( Vec3 const & p ) const { return sqrt( squareDistance(p) ); }

    bool isParallelTo( Line const & L ) const {
        
        return fabs(Vec3::dot(L.direction(),m_normal)) < 1e-6; // Pas 0 car imprécision sur les flottants en machine
        
    }

    Vec3 getIntersectionPoint( Line const & L ) const {
        
        Vec3 result = Vec3(-100.0f,-100.0f,-100.0f);

        if(!isParallelTo(L)) {

            float D = Vec3::dot(m_center,m_normal);
            float t = (D - Vec3::dot(L.origin(),m_normal)) / Vec3::dot(L.direction(),m_normal);

            if(t > 0) {

                result = L.origin() + t * L.direction();

            }

        }

        return result;

    }
};
#endif
