#include "Mesh.h"
#include <iostream>
#include <fstream>

void Mesh::loadOFF (const std::string & filename) {
    std::ifstream in (filename.c_str ());
    if (!in)
        exit (EXIT_FAILURE);
    std::string offString;
    unsigned int sizeV, sizeT, tmp;
    in >> offString >> sizeV >> sizeT >> tmp;
    vertices.resize (sizeV);
    triangles.resize (sizeT);
    for (unsigned int i = 0; i < sizeV; i++)
        in >> vertices[i].position;
    int s;
    for (unsigned int i = 0; i < sizeT; i++) {
        in >> s;
        for (unsigned int j = 0; j < 3; j++)
            in >> triangles[i].v[j];
    }
    in.close ();
}

// void Mesh::loadOFF(const std::string &filename) {
//     std::ifstream in(filename.c_str());
//     if (!in) {
//         std::cerr << "Cannot open file: " << filename << std::endl;
//         exit(EXIT_FAILURE);
//     }
    
//     std::string offString;
//     unsigned int sizeV, sizeT, tmp;
//     in >> offString;
//     if (offString != "OFF") {
//         std::cerr << "File is not a valid OFF file: " << filename << std::endl;
//         exit(EXIT_FAILURE);
//     }
    
//     in >> sizeV >> sizeT >> tmp;
//     if (!in) {
//         std::cerr << "Reading sizes failed!" << std::endl;
//         exit(EXIT_FAILURE);
//     }

//     vertices.resize(sizeV);
//     for (unsigned int i = 0; i < sizeV; i++) {
//         if (!(in >> vertices[i].position)) {
//             std::cerr << "Reading vertex " << i << " failed!" << std::endl;
//             exit(EXIT_FAILURE);
//         }
//     }

//     int s;
//     for (unsigned int i = 0; i < sizeT; i++) {
//         in >> s;
//         if (s != 3) {
//             std::cerr << "Face " << i << " is not a triangle!" << std::endl;
//             exit(EXIT_FAILURE);
//         }
//         for (unsigned int j = 0; j < 3; j++) {
//             if (!(in >> triangles[i].v[j])) {
//                 std::cerr << "Reading triangle " << i << " failed!" << std::endl;
//                 exit(EXIT_FAILURE);
//             }
//         }
//     }

//     in.close();
//     if (!in) {
//         std::cerr << "File was not closed properly!" << std::endl;
//         exit(EXIT_FAILURE);
//     }
// }


void Mesh::recomputeNormals () {
    for (unsigned int i = 0; i < vertices.size (); i++)
        vertices[i].normal = Vec3 (0.0, 0.0, 0.0);
    for (unsigned int i = 0; i < triangles.size (); i++) {
        Vec3 e01 = vertices[triangles[i].v[1]].position -  vertices[triangles[i].v[0]].position;
        Vec3 e02 = vertices[triangles[i].v[2]].position -  vertices[triangles[i].v[0]].position;
        Vec3 n = Vec3::cross (e01, e02);
        n.normalize ();
        for (unsigned int j = 0; j < 3; j++)
            vertices[triangles[i].v[j]].normal += n;
    }
    for (unsigned int i = 0; i < vertices.size (); i++)
        vertices[i].normal.normalize ();
}

void Mesh::centerAndScaleToUnit () {
    Vec3 c(0,0,0);
    for  (unsigned int i = 0; i < vertices.size (); i++)
        c += vertices[i].position;
    c /= vertices.size ();
    float maxD = (vertices[0].position - c).length();
    for (unsigned int i = 0; i < vertices.size (); i++){
        float m = (vertices[i].position - c).length();
        if (m > maxD)
            maxD = m;
    }
    for  (unsigned int i = 0; i < vertices.size (); i++)
        vertices[i].position = (vertices[i].position - c) / maxD;
}
