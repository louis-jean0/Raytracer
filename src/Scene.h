#ifndef SCENE_H
#define SCENE_H

#include <cfloat>
#include <cstdlib>
#include <limits>
#include <locale>
#include <vector>
#include <string>
#include "Camera.h"
#include "Material.h"
#include "Mesh.h"
#include "Sphere.h"
#include "Square.h"
#include "Triangle.h"
#include "Vec3.h"
#include <GL/glut.h>

enum LightType {
    LightType_Spherical,
    LightType_Quad
};

enum MeshType {
    MeshType_Mesh,
    MeshType_Sphere,
    MeshType_Square
};

struct Light {
    Vec3 material;
    bool isInCamSpace;
    LightType type;

    Vec3 pos;
    float radius;

    Mesh quad;

    float powerCorrection;
    float ambientPower;
    float diffusePower;
    float specularPower;

    Light() : powerCorrection(1.0) {}

};

struct RaySceneIntersection{
    bool intersectionExists;
    unsigned int typeOfIntersectedObject;
    unsigned int objectIndex;
    float t;
    RayTriangleIntersection rayMeshIntersection;
    RaySphereIntersection raySphereIntersection;
    RaySquareIntersection raySquareIntersection;
    RaySceneIntersection() : intersectionExists(false) , t(FLT_MAX) {}
};



class Scene {
    std::vector< Mesh > meshes;
    std::vector< Sphere > spheres;
    std::vector< Square > squares;
    std::vector< Light > lights;

public:

    Scene() {
    }

    void draw() {
        // iterer sur l'ensemble des objets, et faire leur rendu :
        for( unsigned int It = 0 ; It < meshes.size() ; ++It ) {
            Mesh const & mesh = meshes[It];
            mesh.draw();
        }
        for( unsigned int It = 0 ; It < spheres.size() ; ++It ) {
            Sphere const & sphere = spheres[It];
            sphere.draw();
        }
        for( unsigned int It = 0 ; It < squares.size() ; ++It ) {
            Square const & square = squares[It];
            square.draw();
        }
    }

    Vec3 sampleAreaLight(const Light &light, float squareSide) {

        Vec3 center = light.pos;

        Vec3 normal(0.0f,0.0f,-1.0f);
        Vec3 upVector(0.0f,1.0f,0.0f);
        Vec3 rightVector(1.0f,0.0f,0.0f);

        upVector = upVector * (squareSide / 2.0f);
        rightVector = rightVector * (squareSide / 2.0f);

        float u = (float)rand()/(float)(RAND_MAX) - 0.5f;
        float v = (float)rand()/(float)(RAND_MAX) - 0.5f;

        Vec3 sampledPoint = center + (u * rightVector) + (v * upVector);

        return sampledPoint;

    }

    RaySceneIntersection computeIntersection(Ray const & ray) {

        RaySceneIntersection result;
        result.intersectionExists = false;
        float closestIntersection = FLT_MAX;
        int nbSpheres = spheres.size();
        int nbSquares = squares.size();
        int nbMeshes = meshes.size();

        // Traitement des sphères

        for(int i = 0; i < nbSpheres; i++) {

            const Sphere sphere = spheres[i];

            RaySphereIntersection sphereIntersection = sphere.intersect(ray);

            if(sphereIntersection.intersectionExists && sphereIntersection.t < closestIntersection) {

                closestIntersection = sphereIntersection.t;
                result.intersectionExists = true;
                result.t = sphereIntersection.t;
                result.objectIndex = i;
                result.typeOfIntersectedObject = MeshType_Sphere;
                result.raySphereIntersection = sphereIntersection;

            }

        }

        // Traitement des carrés

        for(int j = 0; j < nbSquares; j++) {

            const Square square = squares[j];

            RaySquareIntersection squareIntersection = square.intersect(ray);

            if(squareIntersection.intersectionExists && squareIntersection.t < closestIntersection) {

                closestIntersection = squareIntersection.t;
                result.intersectionExists = true;
                result.t = squareIntersection.t;
                result.objectIndex = j;
                result.typeOfIntersectedObject = MeshType_Square;
                result.raySquareIntersection = squareIntersection;

            }

        }

        // Traitement des meshes

        for(int k = 0; k < nbMeshes; k++) {

            const Mesh mesh = meshes[k];
            RayTriangleIntersection triangleIntersection = mesh.intersect(ray);

            if(triangleIntersection.intersectionExists && triangleIntersection.t < closestIntersection) {

                closestIntersection = triangleIntersection.t;
                result.intersectionExists = true;
                result.t = triangleIntersection.t;
                result.objectIndex = k;
                result.typeOfIntersectedObject = MeshType_Mesh;
                result.rayMeshIntersection = triangleIntersection;

            }


        }

        return result;

    }

    Material getMaterial(RaySceneIntersection r) {

        int type = r.typeOfIntersectedObject;
        
        if(type == MeshType_Mesh){
            return meshes[0].material;
        }

        if(type == MeshType_Sphere) {
            return spheres[r.objectIndex].material;
        }
        
        if(type == MeshType_Square) {
            return squares[r.objectIndex].material;
        }

        return Material();

    }

    Vec3 rayTraceRecursive(Ray ray, int NRemainingBounces) {
        
        RaySceneIntersection raySceneIntersection = computeIntersection(ray);
        Vec3 intersectionPoint(0.0f,0.0f,0.0f);
        Vec3 normalToIntersectionPoint(0.0f,0.0f,0.0f);
        Vec3 reflectedDirection(0.0f,0.0f,0.0f);
        Vec3 color(1.0f,1.0f,1.0f); // Fond blanc
        //Vec3 color(0.0f,0.0f,0.0f); // Fond noir
        Vec3 ambient(0.0f,0.0f,0.0f);
        Vec3 diffuse(0.0f,0.0f,0.0f);
        Vec3 specular(0.0f,0.0f,0.0f);
        Vec3 N(0.0f,0.0f,0.0f);
        Vec3 L(0.0f,0.0f,0.0f);
        Vec3 R(0.0f,0.0f,0.0f);
        Vec3 V(0.0f,0.0f,0.0f);
        Sphere currentSphere;
        Square currentSquare;
        Mesh currentMesh;
        Ray shadowRay;
        RaySceneIntersection shadowIntersection;
        float epsilon = 0.001f;
        float distanceToLight;
        int nbSamplesSmoothShadows = 5;
        float shadowCounter = 0.0f;
        float squareSide = 1.0f;
        Vec3 sampledPointPosition(0.0f,0.0f,0.0f);
        Vec3 lightDirectionForSampledPoint(0.0f,0.0f,0.0f);

        if(!raySceneIntersection.intersectionExists) {
            return color;
        }

        switch(raySceneIntersection.typeOfIntersectedObject) {

            case MeshType_Sphere:
                currentSphere = spheres[raySceneIntersection.objectIndex];
                intersectionPoint = raySceneIntersection.raySphereIntersection.intersection;
                normalToIntersectionPoint = raySceneIntersection.raySphereIntersection.normal;
            break;

            case MeshType_Square:
                currentSquare = squares[raySceneIntersection.objectIndex];
                intersectionPoint = raySceneIntersection.raySquareIntersection.intersection;
                normalToIntersectionPoint = raySceneIntersection.raySquareIntersection.normal;
            break;

            case MeshType_Mesh:
            {
                currentMesh = meshes[raySceneIntersection.objectIndex];
                intersectionPoint = raySceneIntersection.rayMeshIntersection.intersection;
                normalToIntersectionPoint = raySceneIntersection.rayMeshIntersection.normal;
            }
            break;

            default:
            break;

        }

        N = normalToIntersectionPoint;
        L = lights[0].pos - intersectionPoint;
        R = (2 * Vec3::dot(L,N) * N) - L;
        V = ray.origin() - intersectionPoint;
        N.normalize();
        L.normalize();
        R.normalize();
        V.normalize();

        // Ombres douces
        for(int i = 0; i < nbSamplesSmoothShadows; i++) {
            sampledPointPosition = sampleAreaLight(lights[0], squareSide);
            lightDirectionForSampledPoint = sampledPointPosition - intersectionPoint;
            distanceToLight = lightDirectionForSampledPoint.length();
            lightDirectionForSampledPoint.normalize();
            shadowRay = Ray(intersectionPoint + epsilon * N, lightDirectionForSampledPoint);
            shadowIntersection = computeIntersection(shadowRay);
            if(!(shadowIntersection.intersectionExists && shadowIntersection.t < distanceToLight)) {
                    shadowCounter += 1.0f;
            }
        }

        // Calcul de la direction du rayon réfléchi
        Vec3 direction = ray.direction();
        direction.normalize();
        reflectedDirection = direction - (2 * Vec3::dot(direction, normalToIntersectionPoint) * normalToIntersectionPoint);
        reflectedDirection.normalize();

        Material m = getMaterial(raySceneIntersection);

        if (m.type == Material_Mirror) {

            if (NRemainingBounces > 0) {
                color = rayTraceRecursive(Ray(intersectionPoint + epsilon * N, reflectedDirection), NRemainingBounces - 1);
            } else {
                color = Vec3(0.0f, 0.0f, 0.0f);
            }

        } 

        else {
    
            ambient = m.ambient_material * lights[0].ambientPower;
            diffuse = m.diffuse_material * lights[0].diffusePower * Vec3::dot(L,N);
            specular = m.specular_material * lights[0].specularPower * pow(Vec3::dot(R,V),m.shininess);
            shadowCounter /= nbSamplesSmoothShadows;

            diffuse *= shadowCounter;
            specular *= shadowCounter;

            color = ambient + diffuse + specular;
            
            if(NRemainingBounces > 0) {
                rayTraceRecursive(Ray(intersectionPoint,reflectedDirection),NRemainingBounces - 1);
            }

        }

        return color;
        
    }

    Vec3 rayTrace(Ray const & rayStart) {

        return rayTraceRecursive(rayStart,2);

    }

    void setup_single_sphere() {
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();

        {
            lights.resize( lights.size() + 1 );
            Light & light = lights[lights.size() - 1];
            light.pos = Vec3(-5,5,5);
            light.radius = 2.5f;
            light.powerCorrection = 2.f;
            light.ambientPower = 1.0f;
            light.diffusePower = 1.0f;
            light.specularPower = 0.4f;
            light.type = LightType_Spherical;
            light.material = Vec3(1,1,1);
            light.isInCamSpace = false;
        }
        {
            spheres.resize( spheres.size() + 1 );
            Sphere & s = spheres[spheres.size() - 1];
            s.m_center = Vec3(1.0 , 0. , 0.);
            s.m_radius = 1.f;
            s.build_arrays();
            s.material.type = Material_Mirror;
            s.material.diffuse_material = Vec3( 1.0,0.0,0.0 );
            s.material.specular_material = Vec3( 0.2,0.2,0.2 );
            s.material.shininess = 20;
        }
        {
            spheres.resize(spheres.size() + 1);
            Sphere &s2 = spheres[spheres.size() - 1];
            s2.m_center = Vec3(-1.0f,0.0f,0.0f);
            s2.m_radius = 1.0f;
            s2.build_arrays();
            s2.material.type = Material_Diffuse_Blinn_Phong;
            s2.material.diffuse_material = Vec3(1.0f,1.0f,0.0f);
            s2.material.specular_material = Vec3(0.2f,0.2f,0.2f);
            s2.material.shininess = 20;
        }
    }

    void setup_single_square() {
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();

        {
            lights.resize( lights.size() + 1 );
            Light & light = lights[lights.size() - 1];
            light.pos = Vec3(-5, 5, 5);
            light.radius = 2.5f;
            light.powerCorrection = 2.f;
            light.ambientPower = 0.5f;
            light.diffusePower = 0.5f;
            light.specularPower = 1.0f;
            light.type = LightType_Spherical;
            light.material = Vec3(1,1,1);
            light.isInCamSpace = false;
        }
        {
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.scale(Vec3(2.0f,2.0f,1.0f));
            s.build_arrays();
            s.material.diffuse_material = Vec3( 1,0.4,0.4 );
            s.material.specular_material = Vec3( 0.8,0.8,0.8 );
            s.material.shininess = 20;
        }

        {
            // squares.resize( squares.size() + 1 );
            // Square & s2 = squares[squares.size() - 1];
            // s2.setQuad(Vec3(0., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            // s2.rotate_x(-45);
            // s2.build_arrays();
            // s2.material.diffuse_material = Vec3( 0.2,1,0.8 );
            // s2.material.specular_material = Vec3( 0.8,0.8,0.8 );
            // s2.material.shininess = 20;
        }

        }

    void setup_cornell_box(){
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();

        {
            lights.resize( lights.size() + 1 );
            Light & light = lights[lights.size() - 1];
            light.pos = Vec3(0,1.5f,0);
            light.radius = 2.5f;
            light.powerCorrection = 2.f;
            light.ambientPower = 1.0f;
            light.diffusePower = 1.0f;
            light.specularPower = 0.6f;
            light.type = LightType_Spherical;
            light.material = Vec3(1,1,1);
            light.isInCamSpace = false;
        }

        { //Back Wall
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.scale(Vec3(2., 2., 1.));
            s.translate(Vec3(0., 0., -2.));
            s.build_arrays();
            s.material.diffuse_material = Vec3( 0.3,0.2,0.5);
            s.material.specular_material = Vec3( 1.,1.,1. );
            s.material.shininess = 16;
            s.material.type = Material_Mirror;
        }

        { //Left Wall

            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.scale(Vec3(2., 2., 1.));
            s.translate(Vec3(0., 0., -2.));
            s.rotate_y(90);
            s.build_arrays();
            s.material.diffuse_material = Vec3( 1.,0.,0. );
            s.material.specular_material = Vec3( 1.,0.,0. );
            s.material.shininess = 16;
        }

        { //Right Wall
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_y(-90);
            s.build_arrays();
            s.material.diffuse_material = Vec3( 0.0,1.0,0.0 );
            s.material.specular_material = Vec3( 0.0,1.0,0.0 );
            s.material.shininess = 16;
        }

        { //Floor
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_x(-90);
            s.build_arrays();
            s.material.diffuse_material = Vec3( 1.0,0.0,1.0 );
            s.material.specular_material = Vec3( 1.0,1.0,1.0 );
            s.material.shininess = 16;
        }

        { //Ceiling
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_x(90);
            s.build_arrays();
            s.material.diffuse_material = Vec3( 0.0,0.0,1);
            s.material.specular_material = Vec3( 1.0,1.0,1.0 );
            s.material.shininess = 16;
        }

        // { //Front Wall
        //     squares.resize( squares.size() + 1 );
        //     Square & s = squares[squares.size() - 1];
        //     s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
        //     s.translate(Vec3(0., 0., -3));
        //     s.scale(Vec3(2., 2., 1.));
        //     s.rotate_y(180);
        //     s.build_arrays();
        //     s.material.diffuse_material = Vec3( 1.0,0.0,0.0 );
        //     s.material.specular_material = Vec3( 1.0,1.0,1.0 );
        //     s.material.shininess = 16;
        // }


        { //GLASS Sphere

            spheres.resize( spheres.size() + 1 );
            Sphere & s = spheres[spheres.size() - 1];
            s.m_center = Vec3(1.0, -1.25, 0.5);
            s.m_radius = 0.75f;
            s.build_arrays();
            s.material.type = Material_Mirror;
            s.material.diffuse_material = Vec3( 1.,0.,0. );
            s.material.specular_material = Vec3( 1.,0.,0. );
            s.material.shininess = 5;
            s.material.transparency = 1.0;
            s.material.index_medium = 1.4;
        }


        { //MIRRORED Sphere
            spheres.resize( spheres.size() + 1 );
            Sphere & s = spheres[spheres.size() - 1];
            s.m_center = Vec3(-1.0, -1.25, -0.5);
            s.m_radius = 0.75f;
            s.build_arrays();
            s.material.type = Material_Mirror;
            s.material.diffuse_material = Vec3( 1.,1.,1. );
            s.material.specular_material = Vec3(  1.,1.,1. );
            s.material.shininess = 10;
            s.material.transparency = 0.;
            s.material.index_medium = 0.;
        }

    }

    void setup_single_triangle() {
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();
        {
            lights.resize( lights.size() + 1 );
            Light & light = lights[lights.size() - 1];
            light.pos = Vec3(-5, 5, 5);
            light.radius = 2.5f;
            light.powerCorrection = 2.f;
            light.ambientPower = 0.5f;
            light.diffusePower = 0.5f;
            light.specularPower = 1.0f;
            light.type = LightType_Spherical;
            light.material = Vec3(1,1,1);
            light.isInCamSpace = false;
        } 
        {
            MeshVertex v0,v1,v2;
            v0.position = Vec3(0.0f,0.0f,0.0f);
            v1.position = Vec3(1.0f, 0.0f, 0.0f);
            v2.position = Vec3(0.0f, 1.0f, 0.0f);
            v0.normal = Vec3(0.0f, 0.0f, 1.0f);
            v1.normal = Vec3(0.0f, 0.0f, 1.0f);
            v2.normal = Vec3(0.0f, 0.0f, 1.0f);
            MeshTriangle triangle;
            triangle.v[0] = 0;
            triangle.v[1] = 1;
            triangle.v[2] = 2;
            meshes.resize(meshes.size() + 1);
            Mesh &mesh = meshes[meshes.size() - 1];
            mesh.vertices.push_back(v0);
            mesh.vertices.push_back(v1);
            mesh.vertices.push_back(v2);
            mesh.triangles.push_back(triangle);
            mesh.build_arrays();
        }
    }

    void setup_simple_mesh() {
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();
        {
            lights.resize( lights.size() + 1 );
            Light & light = lights[lights.size() - 1];
            light.pos = Vec3(-5, 5, 5);
            light.radius = 2.5f;
            light.powerCorrection = 2.f;
            light.ambientPower = 0.5f;
            light.diffusePower = 0.5f;
            light.specularPower = 1.0f;
            light.type = LightType_Spherical;
            light.material = Vec3(1,1,1);
            light.isInCamSpace = false;
        }
        // { //Back Wall
        //     squares.resize( squares.size() + 1 );
        //     Square & s = squares[squares.size() - 1];
        //     s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
        //     s.scale(Vec3(2., 2., 1.));
        //     s.translate(Vec3(0., 0., -2.));
        //     s.build_arrays();
        //     s.material.diffuse_material = Vec3( 0.3,0.2,0.5);
        //     s.material.specular_material = Vec3( 1.,1.,1. );
        //     s.material.shininess = 16;
        // }
        {
            meshes.resize(meshes.size() + 1);
            Mesh &mesh = meshes[meshes.size() - 1];
            mesh.loadOFF("models/pipe.off");
            mesh.build_arrays();
            // mesh.recomputeNormals();
            // mesh.centerAndScaleToUnit();
            mesh.material.ambient_material = Vec3(0.4f,0.2f,0.1f);
            mesh.material.diffuse_material = Vec3(0.5f,0.5f,0.5f);
            mesh.material.specular_material = Vec3(0.1f,0.2f,0.2f);
        }
    }

            

};



#endif
