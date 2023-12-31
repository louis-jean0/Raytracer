#include "Mesh.h"
#include "Triangle.h"
#include "KdNode.h"
#include <cstddef>
#include <stdexcept>
#include <vector>
#include <algorithm>
#define MAX_DEPTH 500
#define MIN_TRIANGLES 1

class KdTree {
    public:
        KdTree() : root(nullptr) {}

        ~KdTree() {
            deleteTree(root);
        }

        bool isEmpty() {
            return root == nullptr;
        }

        void drawAABBs() {
            drawNodeAABB(root);
        }

        void build(const std::vector<Mesh> &meshes) {
            if(root != nullptr) {
                deleteTree(root);
            }

            std::vector<Triangle> allTriangles;

            for (const auto& mesh : meshes) {
                auto meshTriangles = mesh.getTriangles();
                allTriangles.insert(allTriangles.end(), meshTriangles.begin(), meshTriangles.end());
            }

            AABB globalBox = AABB::computeGlobalBox(allTriangles);

            root = buildRecursiveKdTree(allTriangles, globalBox, 0);

        }

        void sortTrianglesByAxis(std::vector<Triangle> &triangles, int axis) {
            if (triangles.size() <= 1) return;
            // Lambda fonction pour trier les centres des triangles selon l'axe considéré
            std::sort(triangles.begin(), triangles.end(), [axis](const Triangle a, const Triangle b) -> bool {
                Vec3 centerA = a.getCenter();
                Vec3 centerB = b.getCenter();
                return centerA[axis] < centerB[axis];
            });
        }

        float computeMedian(const std::vector<Triangle> &triangles, int axis) {
            if (triangles.empty()) {
                return 0.0f;
            }

            std::vector<float> values;
            values.resize(triangles.size());

            for (const Triangle &tri : triangles) {
                Vec3 center = tri.getCenter();
                values.push_back(center[axis]);
            }

            size_t n = values.size() / 2;
            std::nth_element(values.begin(), values.begin() + n, values.end()); // Pour trouver la médiane sans avoir à tout trier
            float median = values[n];

            // Si le nombre de triangles est pair, trouvez les deux médianes centrales et prenez la moyenne.
            if (values.size() % 2 == 0) {
                std::nth_element(values.begin(), values.begin() + n - 1, values.begin() + n);
                median = (median + values[n - 1]) * 0.5f;
            }
            return median;
        }

        void divideTrianglesAtMedian(const std::vector<Triangle>& triangles, int axis, float median, std::vector<Triangle>& leftChildTriangles, std::vector<Triangle>& rightChildTriangles) {
            for (const Triangle& tri : triangles) {
                Vec3 center = tri.getCenter();
                if (center[axis] <= median) {
                    leftChildTriangles.push_back(tri);
                } else {
                    rightChildTriangles.push_back(tri);
                }
            }
        }

        bool intersect(const Ray &ray, RayTriangleIntersection &intersection) {
            return traverseTree(root, ray, intersection);
        }

        bool intersectTriangles(KdNode* node, const Ray& ray, RayTriangleIntersection& intersection) {
            bool hit = false;
            float closestDistance = FLT_MAX;

            for (const Triangle &tri : node->getTriangles()) {
                RayTriangleIntersection tempIntersection = tri.getIntersection(ray);
                if (tempIntersection.intersectionExists && tempIntersection.t < closestDistance) {
                    closestDistance = tempIntersection.t;
                    intersection = tempIntersection;
                    hit = true;
                }
            }
            return hit;
        }

    private:
        KdNode* root;

        void drawNodeAABB(KdNode* node) {
            if(node == nullptr) return;

            node->getBox().drawBox();
            drawNodeAABB(node->getLeftChild());
            drawNodeAABB(node->getRightChild());
        }
        
        void deleteTree(KdNode* node) {
            if(node != nullptr) {
                deleteTree(node->getLeftChild());
                deleteTree(node->getRightChild());
                delete node;
            }
        }

        KdNode* buildRecursiveKdTree(const std::vector<Triangle> &triangles, const AABB& box, int depth) {
            KdNode* node = new KdNode();
            node->setBox(box);
            std::cout << "Profondeur : " << depth << ", Nombre de triangles : " << triangles.size() << std::endl;

            if(depth >= MAX_DEPTH || triangles.size() <= MIN_TRIANGLES) {
                node->setLeaf(true);
                node->setTriangles(triangles);
                return node;
            }

            int currentAxis = depth % 3; // Pour alterner selon les trois axes
            
            // Triez les triangles en fonction de l'axe actuel et calculez la médiane
            std::vector<Triangle> sortedTriangles = triangles;
            sortTrianglesByAxis(sortedTriangles, currentAxis);
            float median = computeMedian(sortedTriangles, currentAxis);

            size_t medianIndex = sortedTriangles.size() / 2;
            Triangle medianTriangle = sortedTriangles[medianIndex];

            std::vector<Triangle> leftChildTriangles(sortedTriangles.begin(), sortedTriangles.begin() + medianIndex);
            std::vector<Triangle> rightChildTriangles(sortedTriangles.begin() + medianIndex, sortedTriangles.end());

            // Divisez les triangles en fonction de la médiane
            // std::vector<Triangle> leftChildTriangles;
            // std::vector<Triangle> rightChildTriangles;
            // divideTrianglesAtMedian(triangles, median, currentAxis, leftChildTriangles, rightChildTriangles);

            // Calculez les boîtes englobantes pour les enfants gauche et droit
            AABB leftChildBox = AABB::computeGlobalBox(leftChildTriangles);
            AABB rightChildBox = AABB::computeGlobalBox(rightChildTriangles);

            // Construisez récursivement les sous-arbres pour les enfants gauche et droit
            node->setLeftChild(buildRecursiveKdTree(leftChildTriangles, leftChildBox, depth + 1));
            node->setRightChild(buildRecursiveKdTree(rightChildTriangles, rightChildBox, depth + 1));
            // std::cout<<"Noeud gauche : "<<node->getLeftChild()->getTriangles().size()<<std::endl;
            // std::cout<<"Noeud droit : "<<node->getRightChild()->getTriangles().size()<<std::endl;

            // Mise à jour des AABB du parent
            if (node->getLeftChild() && node->getRightChild()) {
                // Englobe les boîtes des enfants gauche et droit
                node->setBox(AABB::combine(node->getLeftChild()->getBox(), node->getRightChild()->getBox()));
            } else if (node->getLeftChild()) {
                // Si seulement un enfant gauche, utilisez son AABB
                node->setBox(node->getLeftChild()->getBox());
            } else if (node->getRightChild()) {
                // Si seulement un enfant droit, utilisez son AABB
                node->setBox(node->getRightChild()->getBox());
            }

            return node;
        }


        bool traverseTree(KdNode* node, const Ray &ray, RayTriangleIntersection &intersection) {
            if(node == nullptr) {
                return false;
            }

            if(node->getIsLeaf()) {
                return intersectTriangles(node, ray, intersection);
            }

            if(!node->getBox().intersect(ray)) {
                return false;
            }

            bool hitLeft = traverseTree(node->getLeftChild(), ray, intersection);
            bool hitRight = traverseTree(node->getRightChild(), ray, intersection);

            return hitLeft || hitRight;
        }
        
};