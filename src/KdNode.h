#include "AABB.h"
#include "Mesh.h"
#include "Triangle.h"
#include <vector>

class KdNode {
    public:
        KdNode() : isLeaf(false), leftChild(nullptr), rightChild(nullptr) {}
        
        ~KdNode() {} // Destruction gérée dans KdTree

        AABB getBox() {
            return this->box;
        }

        void setBox(AABB box) {
            this->box = box;
        }

        KdNode* getLeftChild() {
            return this->leftChild;
        }

        void setLeftChild(KdNode* leftChild) {
            this->leftChild = leftChild;
        }

        KdNode* getRightChild() {
            return this->rightChild;
        }

        void setRightChild(KdNode* rightChild) {
            this->rightChild = rightChild;
        }

        void setLeaf(bool isLeaf) {
            this->isLeaf = isLeaf;
        }

        bool getIsLeaf() {
            return this->isLeaf;
        }

        void setTriangles(const std::vector<Triangle> &triangles) {
            this->triangles = triangles;
        }

        std::vector<Triangle> &getTriangles() {
            return this->triangles;
        }

        const std::vector<Triangle>& getTriangles() const {
            return this->triangles;
        }

    private:
        AABB box;
        KdNode* leftChild;
        KdNode* rightChild;
        bool isLeaf;
        std::vector<Triangle> triangles;
};