#ifndef SCENE_H
#define SCENE_H

#include "Mesh.h"
#include <GL/gl.h>
#include <cfloat>
#include <cstdlib>
#include <limits>
#include <locale>
#include <vector>
#include <string>
#include "Camera.h"
#include "Material.h"
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
    bool drawABBS = false;

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

        if(drawABBS) {
            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
            glColor3f(1.0f, 0.0f, 0.0f); // Couleur rouge pour les AABB
            glDisable(GL_LIGHTING);
            glEnable(GL_LINE_SMOOTH);

            kdTree.drawAABBs();

            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
            glDisable(GL_LINE_SMOOTH);
            glEnable(GL_LIGHTING);
        }

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

        // Traitement des meshes avec kd-tree

        RayTriangleIntersection triangleIntersection;
        if(!kdTree.isEmpty()) {
            if (kdTree.intersect(ray, triangleIntersection)) {
                if (triangleIntersection.intersectionExists && triangleIntersection.t < closestIntersection) {
                    closestIntersection = triangleIntersection.t;
                    result.intersectionExists = true;
                    result.t = triangleIntersection.t;
                    result.objectIndex = triangleIntersection.tIndex;
                    result.typeOfIntersectedObject = MeshType_Mesh;
                    result.rayMeshIntersection = triangleIntersection;
                }
            }
        }

        // // Méthode avant kd-tree
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

    // Calcule la direction réfractée du rayon considéré
    Vec3 computeRefractedDirection(const Vec3 &incident, const Vec3 &normal, const float &airRefractionIndice, const float &ior) {
        float cosi = std::clamp(-1.0f, 1.0f, Vec3::dot(incident, normal));
        float etai = airRefractionIndice, etat = ior;
        Vec3 n = normal;
        if (cosi < 0) { cosi = -cosi; } else { std::swap(etai, etat); n= -1*normal; }
        float eta = etai / etat;
        float k = 1 - eta * eta * (1 - cosi * cosi);
        return k < 0 ? Vec3(0.0f,0.0f,0.0f) : eta * incident + (eta * cosi - sqrtf(k)) * n;
    }

    float computeFresnelEffect(const Vec3 &I, const Vec3 &N, const float &airRefractionIndice, const float &ior) {
        float cosi = std::clamp(-1.0f, 1.0f, Vec3::dot(I, N));
        float etai = airRefractionIndice, etat = ior;
        if (cosi > 0) { std::swap(etai, etat); }
        // Calcul de sin(theta_t)^2 via la loi de Snell-Descartes
        float sint = etai / etat * sqrtf(std::max(0.f, 1 - cosi * cosi));
        // Gestion totale de la réflexion interne
        if (sint >= 1) {
            return 1.0f;
        }
        else {
            float cost = sqrtf(std::max(0.f, 1 - sint * sint));
            cosi = fabsf(cosi);
            float Rs = ((etat * cosi) - (etai * cost)) / ((etat * cosi) + (etai * cost));
            float Rp = ((etai * cosi) - (etat * cost)) / ((etai * cosi) + (etat * cost));
            return (Rs * Rs + Rp * Rp) / 2;
        }
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
        }
        
        if (m.type == Material_Glass && NRemainingBounces > 0) {

            float airRefractionIndice = 1.00029f;
            float materialRefractionIndice = m.index_medium;

            Vec3 refractedDirection = computeRefractedDirection(rayDirection, normalToIntersectionPoint, airRefractionIndice, materialRefractionIndice);
            refractedDirection.normalize();

            Vec3 refractedColor = Vec3(0.0f, 0.0f, 0.0f);
            Vec3 reflectedColor = Vec3(0.0f, 0.0f, 0.0f);
            float reflectance = computeFresnelEffect(rayDirection, normalToIntersectionPoint, airRefractionIndice, m.index_medium);
            bool totalInternalReflection = refractedDirection == Vec3(0.0f, 0.0f, 0.0f);

            if(!totalInternalReflection) {
                Vec3 refractedRayOrigin = intersectionPoint;
                refractedColor = m.transparency * rayTraceRecursive(Ray(refractedRayOrigin, refractedDirection), NRemainingBounces - 1);
            }
                
            reflectedColor = rayTraceRecursive(Ray(intersectionPoint, reflectedDirection), NRemainingBounces - 1);

            color = reflectance * reflectedColor + (1 - reflectance) * refractedColor;

        }
        
        if (m.type == Material_Diffuse_Blinn_Phong && NRemainingBounces > 0) {
            color = ambient + diffuse + specular;
            //color = Vec3(1.0f,0.0f,0.0f);
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
            //s.setupTextures("img/sphereTextures/s2.ppm", "img/normalMaps/n1.ppm");
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
            s.material.transparency = 0.7;
            s.material.index_medium = 1.05;
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

    void setup_dof_scene(){
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();

        {
            lights.resize( lights.size() + 1 );
            Light & light = lights[lights.size() - 1];
            light.pos = Vec3(0,0,4);
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

        {
            spheres.resize( spheres.size() + 1 );
            Sphere & s = spheres[spheres.size() - 1];
            s.m_radius = 0.25f;
            s.m_center = Vec3(1.0, -1.75, 0.5);
            s.build_arrays();
            s.material.type = Material_Diffuse_Blinn_Phong;
            s.material.diffuse_material = Vec3( 1.,0.,0. );
            s.material.specular_material = Vec3( 1.,0.,0. );
            s.material.shininess = 5;
            s.material.transparency = 1.0;
            s.material.index_medium = 1.05;
        }


         {
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0.0f,0.0f,1.5f));
            //s.scale(Vec3(1, 1, 0.5));
            s.build_arrays();
            s.material.diffuse_material = Vec3( 1.0,1.0,0.0 );
            s.material.specular_material = Vec3( 0.0,1.0,0.0 );
            s.material.shininess = 5;
        }

        { //MIRRORED Sphere
            spheres.resize( spheres.size() + 1 );
            Sphere & s = spheres[spheres.size() - 1];
            s.m_center = Vec3(-1.0, -1.25, -0.5);
            s.m_radius = 0.75f;
            s.build_arrays();
            s.material.type = Material_Diffuse_Blinn_Phong;
            s.material.diffuse_material = Vec3( 1.,1.,1. );
            s.material.specular_material = Vec3(  1.,1.,1. );
            s.material.shininess = 10;
            s.material.transparency = 0.;
            s.material.index_medium = 0.;
        };
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
        {
            meshes.resize(meshes.size() + 1);
            Mesh &mesh = meshes[meshes.size() - 1];
            mesh.loadOFF("models/suzanne.off");
            mesh.translate(Vec3(0.0f,-1.0f,-1.5f));
            mesh.build_arrays();
            mesh.meshIndex = 0;
            mesh.material.type = Material_Diffuse_Blinn_Phong;
            mesh.material.transparency = 1.0;
            mesh.material.index_medium = 10.0;
            mesh.material.shininess = 10;
            mesh.material.ambient_material = Vec3(0.3f,0.1f,0.5f);
            mesh.material.diffuse_material = Vec3(0.3f,0.1f,0.5f);
            mesh.material.specular_material = Vec3(0.3f,0.1f,0.5f);
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
            light.pos = Vec3(-1,0.2,4);
            light.radius = 2.5f;
            light.powerCorrection = 2.f;
            light.ambientPower = 1.0f;
            light.diffusePower = 1.0f;
            light.specularPower = 0.6f;
            light.type = LightType_Spherical;
            light.material = Vec3(1,1,1);
            light.isInCamSpace = false;
        }
        // Back Wall
        {
            squares.resize(squares.size() + 1);
            Square &s = squares[squares.size() - 1];
            s.setQuad(Vec3(-10.0, -1.0, -10.0), Vec3(1.0, 0.0, 0.0), Vec3(0.0, 1.0, 0.0), 40.0, 40.0);
            s.translate(Vec3(0.0f,-1.0f,0.0f));
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
            s.rotate_x(180);
            s.translate(Vec3(0.0f,-3.0f,0.0f));
            s.build_arrays();
            s.material.diffuse_material = Vec3(0.75, 0.75, 0.75);
            s.material.specular_material = Vec3(0.0, 0.0, 0.0);
            s.material.shininess = 0;
        }
        {
            meshes.resize(meshes.size() + 1);
            Mesh &mesh = meshes[meshes.size() - 1];
            mesh.loadOFF("models/blob.off");
            mesh.build_arrays();
            // mesh.recomputeNormals();
            // mesh.centerAndScaleToUnit();
            mesh.material.type = Material_Mirror;
            mesh.material.ambient_material = Vec3(0.4f,0.2f,0.1f);
            mesh.material.diffuse_material = Vec3(0.5f,0.5f,0.5f);
            mesh.material.specular_material = Vec3(0.1f,0.2f,0.2f);
        }
    }
};
#endif
