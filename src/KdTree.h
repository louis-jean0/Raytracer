#include "Mesh.h"
#include "Triangle.h"
#include "KdNode.h"
#include <cstddef>
#include <stdexcept>
#include <vector>
#include <algorithm>
#define MAX_DEPTH 6
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

            if(meshes.empty()) return;

            std::vector<Triangle> allTriangles;
            std::vector<Vec3> verticesNormals;

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
                    std::vector<Vec3> verticesNormals = tri.getVerticesNormals();
                    Vec3 interpolatedNormal = tempIntersection.w0 * verticesNormals[0]
                                            + tempIntersection.w1 * verticesNormals[1]
                                            + tempIntersection.w2 * verticesNormals[2];
                    interpolatedNormal.normalize();
                    tempIntersection.normal = interpolatedNormal;
                    tempIntersection.tIndex = tri.getMeshIndex();
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
            //std::cout << "Profondeur : " << depth << ", Nombre de triangles : " << triangles.size() << std::endl;

            if(depth >= MAX_DEPTH || triangles.size() <= MIN_TRIANGLES) {
                node->setLeaf(true);
                node->setTriangles(triangles);
                return node;
            }

            int currentAxis = depth % 3; // Pour alterner selon les trois axes
            
            // Tri des triangles
            std::vector<Triangle> sortedTriangles = triangles;
            sortTrianglesByAxis(sortedTriangles, currentAxis);
            size_t medianIndex = sortedTriangles.size() / 2;
            Triangle medianTriangle = sortedTriangles[medianIndex];

            std::vector<Triangle> leftChildTriangles(sortedTriangles.begin(), sortedTriangles.begin() + medianIndex);
            std::vector<Triangle> rightChildTriangles(sortedTriangles.begin() + medianIndex, sortedTriangles.end());

            // Boîtes englobantes pour les enfants gauche et droit
            AABB leftChildBox = AABB::computeGlobalBox(leftChildTriangles);
            AABB rightChildBox = AABB::computeGlobalBox(rightChildTriangles);

            // Construction récursive des sous-arbres pour les enfants gauche et droit
            node->setLeftChild(buildRecursiveKdTree(leftChildTriangles, leftChildBox, depth + 1));
            node->setRightChild(buildRecursiveKdTree(rightChildTriangles, rightChildBox, depth + 1));
            // std::cout<<"Noeud gauche : "<<node->getLeftChild()->getTriangles().size()<<std::endl;
            // std::cout<<"Noeud droit : "<<node->getRightChild()->getTriangles().size()<<std::endl;

            // Mise à jour des AABB du parent
            if (node->getLeftChild() && node->getRightChild()) {
                // Englobe les boîtes des enfants gauche et droit
                node->setBox(AABB::combine(node->getLeftChild()->getBox(), node->getRightChild()->getBox()));
            } else if (node->getLeftChild()) {
                node->setBox(node->getLeftChild()->getBox());
            } else if (node->getRightChild()) {
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

            RayTriangleIntersection leftIntersection,rightIntersection;
            bool hitLeft = traverseTree(node->getLeftChild(), ray, leftIntersection);
            bool hitRight = traverseTree(node->getRightChild(), ray, rightIntersection);

            if(hitLeft && hitRight) {
                if(leftIntersection.t < rightIntersection.t) {
                    intersection = leftIntersection;
                } else {
                    intersection = rightIntersection;
                }
                return true;
            } else if(hitLeft) {
                intersection = leftIntersection;
                return true;
            } else if(hitRight) {
                intersection = rightIntersection;
                return true;
            }

            return false;
        }
};