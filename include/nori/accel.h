/*
Group members:

Name			    ID
Han Wang		    47143106
Jiawen Ye			54035818
Yingchen Zhou		55879895


*/
#pragma once

#include <nori/mesh.h>
#include <map>

NORI_NAMESPACE_BEGIN

/**
 * \brief Acceleration data structure for ray intersection queries
 *
 * The current implementation falls back to a brute force loop
 * through the geometry.
 */

struct Node {
    BoundingBox3f myB;
    Node* leftChild;
    Node* rightChild;
    bool isLeaf = false;
    std::vector<Node*> children;
    std::vector<uint32_t> triangles;
};

class Accel {
public:
    /**
     * \brief Register a triangle mesh for inclusion in the acceleration
     * data structure
     *
     * This function can only be used before \ref build() is called
     */
    void addMesh(Mesh *mesh);
    //bool  Accel::cmpI(uint32_t a, uint32_t b);
    //bool  Accel::cmpII(uint32_t a, uint32_t b);
    //bool  Accel::cmpIII(uint32_t a, uint32_t b);
    void remove(uint32_t);

    Node* nodeBuild(BoundingBox3f, std::vector<uint32_t>, int);
    /// Build the acceleration data structure (currently a no-op)
    void build();

    /// Return an axis-aligned box that bounds the scene
    const BoundingBox3f &getBoundingBox() const { return m_bbox; }

    /**
     * \brief Intersect a ray against all triangles stored in the scene and
     * return detailed intersection information
     *
     * \param ray
     *    A 3-dimensional ray data structure with minimum/maximum extent
     *    information
     *
     * \param its
     *    A detailed intersection record, which will be filled by the
     *    intersection query
     *
     * \param shadowRay
     *    \c true if this is a shadow ray query, i.e. a query that only aims to
     *    find out whether the ray is blocked or not without returning detailed
     *    intersection information.
     *
     * \return \c true if an intersection was found
     */
    bool rayIntersect(const Ray3f &ray, Intersection &its, bool shadowRay) const;

private:
    Mesh         *m_mesh = nullptr; ///< Mesh (only a single one for now)
    BoundingBox3f m_bbox;           ///< Bounding box of the entire scene
    Node* root;
    const int MAX_DEPTH = 100;
    
};

NORI_NAMESPACE_END
