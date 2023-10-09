#ifndef PLANE_H
#define PLANE_H
#include "Vec3.h"
#include "Line.h"
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
        //TODO completer
        return result;
    }

    float squareDistance( Vec3 const & p ) const { return (project(p) - p).squareLength(); }

    float distance( Vec3 const & p ) const { return sqrt( squareDistance(p) ); }

    bool isParallelTo( Line const & L ) const {
        
        bool result = true;

        if(fabs(Vec3::dot(L.direction(),m_normal)) < 1e-6) { // Pas 0 car imprécision sur les flottants en machine

            result = false;

        }

        return result;

    }

    Vec3 getIntersectionPoint( Line const & L ) const {
        // you should check first that the line is not parallel to the plane!
        if(!isParallelTo(L)) {
            Vec3 result;
            float D = Vec3::dot(m_center,m_normal);
            float t = (D - Vec3::dot(L.origin(),m_normal)) / Vec3::dot(L.direction(),m_normal);

            if(t > 0) {

                result = L.origin() + t * L.direction();

            }

            return result;

        }
    }
};
#endif
