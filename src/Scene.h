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
#include <algorithm>
#include "KdTree.h"

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
    KdTree kdTree;

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

        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        glColor3f(1.0f, 0.0f, 0.0f); // Couleur rouge pour les AABB
        glEnable(GL_LINE_SMOOTH);

        // Dessinez les AABBs du KdTree
        kdTree.drawAABBs();

        // Restaurer les paramètres de OpenGL
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        glDisable(GL_LINE_SMOOTH);
        
    }

    void buildKdTreeForScene() {
        kdTree.build(meshes);
    }

    // Crée un carré autour de la lumière et choisi un point au hasard dans ce carré (utile pour les ombres douces)
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

            RaySphereIntersection sphereIntersection = spheres[i].intersect(ray);

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

            RaySquareIntersection squareIntersection = squares[j].intersect(ray);

            if(squareIntersection.intersectionExists && squareIntersection.t < closestIntersection) {

                closestIntersection = squareIntersection.t;
                result.intersectionExists = true;
                result.t = squareIntersection.t;
                result.objectIndex = j;
                result.typeOfIntersectedObject = MeshType_Square;
                result.raySquareIntersection = squareIntersection;

            }

        }

        //Traitement des meshes avec kd-tree

        RayTriangleIntersection triangleIntersection;
        if(!kdTree.isEmpty()) {
            if (kdTree.intersect(ray, triangleIntersection)) {
                if (triangleIntersection.intersectionExists && triangleIntersection.t < closestIntersection) {
                    closestIntersection = triangleIntersection.t;
                    result.intersectionExists = true;
                    result.t = triangleIntersection.t;
                    result.objectIndex = 0;
                    result.typeOfIntersectedObject = MeshType_Mesh;
                    result.rayMeshIntersection = triangleIntersection;
                }
            }
        }
        //Méthode avant kd-tree
        // for(int k = 0; k < nbMeshes; k++) {

        //     RayTriangleIntersection triangleIntersection = meshes[k].intersect(ray);

        //     if(triangleIntersection.intersectionExists && triangleIntersection.t < closestIntersection) {

        //         closestIntersection = triangleIntersection.t;
        //         result.intersectionExists = true;
        //         result.t = triangleIntersection.t;
        //         result.objectIndex = k;
        //         result.typeOfIntersectedObject = MeshType_Mesh;
        //         result.rayMeshIntersection = triangleIntersection;

        //     }

        // }

        return result;

    }

    // Pour récupérer le material d'un objet (glass ou mirror)
    Material getMaterial(RaySceneIntersection r) {

        int type = r.typeOfIntersectedObject;
        
        if(type == MeshType_Mesh){
            return meshes[r.objectIndex].material;
        }

        if(type == MeshType_Sphere) {
            return spheres[r.objectIndex].material;
        }
        
        if(type == MeshType_Square) {
            return squares[r.objectIndex].material;
        }

        return Material();

    }

    // Quand il n'y a aucune intersection (couleur de fond donc)
    Vec3 handleNoIntersection() {
        return Vec3(1.0f,1.0f,1.0f); // Fond blanc
        //return Vec3(0.0f,0.0f,0.0f); // Fond noir
    }

    // Switch sur les intersections pour déterminer le type d'objet intersecté
    std::pair<Vec3, Vec3> computeIntersectionDetails(const RaySceneIntersection& intersection) {

        Vec3 intersectionPoint;
        Vec3 normalToIntersectionPoint;

        switch (intersection.typeOfIntersectedObject) {
            case MeshType_Sphere:
                intersectionPoint = intersection.raySphereIntersection.intersection;
                normalToIntersectionPoint = intersection.raySphereIntersection.normal;
                break;
            case MeshType_Square:
                intersectionPoint = intersection.raySquareIntersection.intersection;
                normalToIntersectionPoint = intersection.raySquareIntersection.normal;
                break;
            case MeshType_Mesh:
                intersectionPoint = intersection.rayMeshIntersection.intersection;
                normalToIntersectionPoint = intersection.rayMeshIntersection.normal;
                break;
            default:
                // Gérer les cas non pris en charge si nécessaire
                break;
        }

        return {intersectionPoint, normalToIntersectionPoint};

    }

    // Calcul des composantes de Phong
    void computePhongComponents(const Vec3& N, const Vec3& L, const Vec3& R, const Vec3& V, const RaySceneIntersection& intersection, Vec3& ambient, Vec3& diffuse, Vec3& specular) {
        
        Material m = getMaterial(intersection);
        ambient = m.ambient_material * lights[0].ambientPower;
        diffuse = m.diffuse_material * lights[0].diffusePower * Vec3::dot(L, N);
        specular = m.specular_material * lights[0].specularPower * pow(Vec3::dot(R, V), m.shininess);
        // std::cout<<"Ambient : "<<ambient<<std::endl;
        // std::cout<<"Diffuse : "<<diffuse<<std::endl;
        // std::cout<<"Specular : "<<specular<<std::endl;

    }

    // Calcul des ombres douces
    float computeSoftShadows(const Vec3& intersectionPoint, const Vec3& N, const RaySceneIntersection& intersection) {
        
        float shadowCounter = 0.0f;
        int nbSamplesSmoothShadows = 5;
        float squareSide = 1.0f;
        Vec3 sampledPointPosition;
        Vec3 lightDirectionForSampledPoint;
        float epsilon = 0.001f;

        for (int i = 0; i < nbSamplesSmoothShadows; i++) {
            sampledPointPosition = sampleAreaLight(lights[0], squareSide);
            lightDirectionForSampledPoint = sampledPointPosition - intersectionPoint;
            float distanceToLight = lightDirectionForSampledPoint.length();
            lightDirectionForSampledPoint.normalize();
            
            Ray shadowRay = Ray(intersectionPoint + epsilon * N, lightDirectionForSampledPoint);
            RaySceneIntersection shadowIntersection = computeIntersection(shadowRay);
            
            if (!(shadowIntersection.intersectionExists && shadowIntersection.t < distanceToLight)) {
                shadowCounter += 1.0f;
            }
        }

        shadowCounter /= nbSamplesSmoothShadows;
        return shadowCounter;
    }

    // Calcule la direction réfléchie du rayon considéré
    Vec3 computeReflectedDirection(Vec3& direction, const Vec3& normal) {
        direction.normalize();
        Vec3 normalizedDirection = direction;
        Vec3 reflectedDirection = normalizedDirection - (2 * Vec3::dot(normalizedDirection, normal) * normal);
        reflectedDirection.normalize();
        return reflectedDirection;
    }

    Vec3 computeRefractedDirection(const Vec3 &incident, const Vec3 &normal, float n1, float n2) {
        
        float cosI = std::clamp(Vec3::dot(incident, normal), -1.0f, 1.0f);
        Vec3 n = normal;

        if (cosI < 0) {
            cosI = -cosI;
        } else {
            std::swap(n1, n2);
            n = -1*n;
        }

        float eta = n1 / n2;
        float k = 1 - eta * eta * (1 - cosI * cosI);

        if (k < 0) {
            return Vec3(0.0f, 0.0f, 0.0f); // Réflexion totale interne
        } else {
            return eta * incident + (eta * cosI - sqrtf(k)) * n;
        }
    }

    float computeFresnelEffect(float cosi, float eta) {
        float r0 = (1 - eta) / (1 + eta);
        r0 = r0 * r0;
        return r0 + (1 - r0) * pow(1 - cosi, 5);
    }

    Vec3 rayTraceRecursive(Ray ray, int NRemainingBounces) {

        RaySceneIntersection raySceneIntersection = computeIntersection(ray);

        if (!raySceneIntersection.intersectionExists) {
            return handleNoIntersection();
        }

        auto [intersectionPoint, normalToIntersectionPoint] = computeIntersectionDetails(raySceneIntersection);

        normalToIntersectionPoint.normalize();
        Vec3 N = normalToIntersectionPoint;
        Vec3 L = lights[0].pos - intersectionPoint;
        Vec3 R = (2 * Vec3::dot(L, N) * N) - L;
        Vec3 V = ray.origin() - intersectionPoint;
        L.normalize(); R.normalize(); V.normalize();

        Vec3 ambient, diffuse, specular;
        computePhongComponents(N, L, R, V, raySceneIntersection, ambient, diffuse, specular);

        float shadowCounter = computeSoftShadows(intersectionPoint, N, raySceneIntersection);
        diffuse *= shadowCounter;
        specular *= shadowCounter;

        Vec3 color = handleNoIntersection();

        Material m = getMaterial(raySceneIntersection);

        Vec3 rayDirection = ray.direction();
        rayDirection.normalize();
        Vec3 reflectedDirection = computeReflectedDirection(rayDirection, normalToIntersectionPoint);
        float airRefractionIndice = 1.00029f;
        Vec3 bias = 0.001f * normalToIntersectionPoint;

        if (m.type == Material_Mirror && NRemainingBounces > 0) {
            Vec3 reflectedRayOrigin = intersectionPoint + bias;
            Ray reflectedRay(reflectedRayOrigin,reflectedDirection);
            color = rayTraceRecursive(Ray(reflectedRayOrigin, reflectedDirection), NRemainingBounces - 1);
        } else if (m.type == Material_Glass && NRemainingBounces > 0) {

            Vec3 unitIncident = ray.direction();
            unitIncident.normalize();
            Vec3 unitNormal = normalToIntersectionPoint;
            unitNormal.normalize();
            float cosI = -std::clamp(Vec3::dot(unitIncident, unitNormal), -1.0f, 1.0f);
            float sinI2 = std::max(0.0f, 1.0f - cosI * cosI);
            float eta = airRefractionIndice / m.index_medium;
            float sinT2 = eta * eta * sinI2;

            Vec3 refraction;
            bool test = false;

            if (sinT2 > 1.0) {
                refraction = Vec3(0.0f, 0.0f, 0.0f); // Réflexion totale interne
                test = true;
            }

            float cosT = std::sqrt(std::max(0.0f, 1.0f - sinT2));
            if (cosI > 0) {
                cosT = -cosT;
            }

            if(!test) {
                refraction = eta * unitIncident + (eta * cosI - cosT) * unitNormal;
            }

            reflectedDirection = computeReflectedDirection(refraction, unitNormal);
            
            if(refraction != Vec3(0.0f,0.0f,0.0f)) {
                Ray refractionRay = Ray(intersectionPoint,refraction);
                //std::cout<<"OK2"<<std::endl;
                color = m.transparency * rayTraceRecursive(refractionRay, NRemainingBounces - 1);
            }
            // else {
            //     // Réflexion totale interne
            //     Ray internalReflectionRay(intersectionPoint, reflectedDirection);
            //     color = rayTraceRecursive(internalReflectionRay, NRemainingBounces - 1);
            // }

        } else if (NRemainingBounces > 0) {
            color = ambient + diffuse + specular;
        }

        return color;

    }

    Vec3 rayTrace(Ray const & rayStart) {

        return rayTraceRecursive(rayStart,3);

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
            //s.material.type = Material_Mirror;
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
            s.material.type = Material_Glass;
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
            mesh.loadOFF("models/SeaMonster.off");
            mesh.build_arrays();
            // mesh.recomputeNormals();
            // mesh.centerAndScaleToUnit();
            mesh.material.ambient_material = Vec3(0.4f,0.2f,0.1f);
            mesh.material.diffuse_material = Vec3(0.1f,0.2f,0.3f);
            mesh.material.specular_material = Vec3(0.1f,0.2f,0.2f);
        }
    }

    void setup_large_scene() {

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
        // Back Wall
        {
            squares.resize(squares.size() + 1);
            Square &s = squares[squares.size() - 1];
            s.setQuad(Vec3(-10.0, -1.0, -10.0), Vec3(1.0, 0.0, 0.0), Vec3(0.0, 1.0, 0.0), 40.0, 40.0);
            s.build_arrays();
            s.material.diffuse_material = Vec3(0.75, 0.75, 0.75);
            s.material.specular_material = Vec3(0.0, 0.0, 0.0);
            s.material.shininess = 0;
        }

        // Floor
        {
            squares.resize(squares.size() + 1);
            Square &s = squares[squares.size() - 1];
            s.setQuad(Vec3(-10.0, -1.0, -10.0), Vec3(1.0, 0.0, 0.0), Vec3(0.0, 0.0, 1.0), 40.0, 40.0);
            s.build_arrays();
            s.material.diffuse_material = Vec3(0.75, 0.75, 0.75);
            s.material.specular_material = Vec3(0.0, 0.0, 0.0);
            s.material.shininess = 0;
        }
        {
            spheres.resize(spheres.size() + 1);
            Sphere &s = spheres[spheres.size() - 1];
            s.m_center = Vec3(0.0, 0.5, 0.0);
            s.m_radius = 0.75f;
            s.build_arrays();
            s.material.type = Material_Glass;
            s.material.diffuse_material = Vec3(1.0, 1.0, 1.0); 
            s.material.specular_material = Vec3(1.0, 1.0, 1.0);
            s.material.shininess = 32; 
            s.material.transparency = 1.0;
            s.material.index_medium = 1.5; 
        }
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
         {
            meshes.resize(meshes.size() + 1);
            Mesh &mesh = meshes[meshes.size() - 1];
            mesh.loadOFF("models/pipe.off");
            mesh.translate(Vec3(4.0f,0.0f,0.0f));
            mesh.build_arrays();
            // mesh.recomputeNormals();
            // mesh.centerAndScaleToUnit();
            mesh.material.ambient_material = Vec3(0.4f,0.2f,0.1f);
            mesh.material.diffuse_material = Vec3(0.5f,0.5f,0.5f);
            mesh.material.specular_material = Vec3(0.1f,0.2f,0.2f);
        }

   
        // {
        //     squares.resize(squares.size() + 1);
        //     Square &sq = squares[squares.size() - 1];
        //     sq.setQuad(Vec3(-0.5, -0.5, 0.5), Vec3(1.0, 0.0, 0.0), Vec3(0.0, 1.0, 0.0), 1.0, 1.0); // Adjust size as needed
        //     sq.build_arrays();
        //     sq.material.type = Material_Glass;
        //     sq.material.diffuse_material = Vec3(0.5, 0.5, 0.5); 
        //     sq.material.specular_material = Vec3(1.0, 1.0, 1.0);
        //     sq.material.shininess = 32; 
        //     sq.material.transparency = 1.0;
        //     sq.material.index_medium = 1.5;
        // }
    }
};
#endif
