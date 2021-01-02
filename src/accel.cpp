/*
Group members:

Name			    ID
Han Wang		    47143106
Jiawen Ye			54035818
Yingchen Zhou		55879895


*/

#include <nori/accel.h>
#include <Eigen/Geometry>
#include <stack>
#include <queue>
#include <utility>

NORI_NAMESPACE_BEGIN

void Accel::addMesh(Mesh* mesh) {
    if (m_mesh)
        throw NoriException("Accel: only a single mesh is supported!");
    m_mesh = mesh;
    m_bbox = m_mesh->getBoundingBox();
}

struct cmpI {
    bool operator()(std::pair<uint32_t,Point3f> a, std::pair<uint32_t,Point3f> b) {
        return a.second.z() > b.second.z();
    }
};

struct cmpII {
    bool operator()(std::pair<uint32_t,Point3f> a, std::pair<uint32_t,Point3f> b) {
        return a.second.x() > b.second.x();
    }
};

struct cmpIII {
    bool operator()(std::pair<uint32_t,Point3f> a, std::pair<uint32_t,Point3f> b) {
        return a.second.y() > b.second.y();
    }
};

Node* Accel::nodeBuild(BoundingBox3f bbox, std::vector<uint32_t> triangles, int depth) {
    // std::cout<<"In build"<<std::endl;
    // we suppose Kt=1,ki=2;
    if (triangles.empty())
        return nullptr;
    
    if (triangles.size() < 5) {
        // std::cout << "ttttttttttttttttttttttt" << std::endl;
        Node* leafNode = new Node();
        leafNode->myB = bbox;
        leafNode->leftChild = nullptr; //=======
        leafNode->rightChild = nullptr; //=======
        leafNode->isLeaf = true;
        leafNode->triangles = triangles;
        return leafNode;
    }
    
    // Divide Bouding Box
    std::vector<BoundingBox3f> bboxB;
    std::vector<BoundingBox3f> left, right;

    /* min.x = min(corner.x, midpoint.x);
       min.y = min(corner.y, midpoint.y);
       min.z = min(corner.z, midpoint.z);*/

    // Point3f midpoint = bbox.getCenter();

    std::vector<std::vector<uint32_t>> triangleList(2, std::vector<uint32_t>{});

    if (depth % 3 == 0) { // divide according to x======= z half 
        
        std::priority_queue<std::pair<uint32_t, Point3f>, std::vector<std::pair<uint32_t, Point3f>>, cmpI >  cmpx;
        for (uint32_t tri : triangles) {
            cmpx.push(std::make_pair(tri, m_mesh->getCentroid(tri)));
        }
        int cutIndex = cmpx.size()/2, i = 0;
        // int cutIndex = 1, i = 0; //float cutZ = cmpx[cmpx.size() / 2]).z();
        while (!cmpx.empty()) {
            if (i < cutIndex) {
                triangleList[0].push_back(cmpx.top().first);
            }
            else
                triangleList[1].push_back(cmpx.top().first);
            cmpx.pop();
            i++;
        }
    }
    else if (depth % 3 == 1) {// divide according to y======= x half
        std::priority_queue<std::pair<uint32_t, Point3f>, std::vector<std::pair<uint32_t, Point3f>>, cmpII >  cmpy;
        for (uint32_t tri : triangles) {
            cmpy.push(std::make_pair(tri, m_mesh->getCentroid(tri)));
        }
        // int cutIndex = 1, i = 0; //float cutZ = cmpx[cmpx.size() / 2]).z();
        int cutIndex = cmpy.size()/2, i = 0;
        while (!cmpy.empty()) {
            if (i < cutIndex) {
                triangleList[0].push_back(cmpy.top().first);
            }
            else
                triangleList[1].push_back(cmpy.top().first);
            cmpy.pop();
            i++;
        }
    }
    else {// divide according to z=======  y half
        std::priority_queue<std::pair<uint32_t, Point3f>, std::vector<std::pair<uint32_t, Point3f>>, cmpIII >  cmpz;
        for (uint32_t tri : triangles) {
            cmpz.push(std::make_pair(tri, m_mesh->getCentroid(tri)));
        }
        // int cutIndex = 1, i = 0; //float cutZ = cmpx[cmpx.size() / 2]).z();
        int cutIndex = cmpz.size()/2 , i = 0;
        while (!cmpz.empty()) {
            if (i < cutIndex) {
                triangleList[0].push_back(cmpz.top().first);
            }
            else
                triangleList[1].push_back(cmpz.top().first);
            cmpz.pop();
            i++;
        }
    }

    BoundingBox3f newLeft = BoundingBox3f();

    for (uint32_t box : triangleList[0]) {
        newLeft.expandBy(m_mesh->getBoundingBox(box));
    }

    BoundingBox3f newRight = BoundingBox3f();
    for (uint32_t box : triangleList[1]) {
        newRight.expandBy(m_mesh->getBoundingBox(box));
    }

    // we suppose Kt=1,ki=2;
    float x1 = bbox.max.x() - bbox.min.x();
    float y1 = bbox.max.y() - bbox.min.y();
    float z1 = bbox.max.z() - bbox.min.z();
    float Sav = 2 * x1 * y1 + 2 * y1 * z1 + 2 * x1 * z1;

    float x2 = newLeft.max.x() - newLeft.min.x();
    float y2 = newLeft.max.y() - newLeft.min.y();
    float z2 = newLeft.max.z() - newLeft.min.z();
    float Savl = 2 * x2 * y2 + 2 * y2 * z2 + 2 * x2 * z2;

    float x3 = newRight.max.x() - newRight.min.x();
    float y3 = newRight.max.y() - newRight.min.y();
    float z3 = newRight.max.z() - newRight.min.z();
    float Sav2 = 2 * x3 * y3 + 2 * y3 * z3+ 2 * x3 * z3;
    float KT = 4, KI = 2;
    float Cvp = KT + KI * (Savl / Sav * float(triangleList[0].size()) + Sav2 / Sav * float(triangleList[1].size()));

    // std::cout <<"CVP"<< Cvp << std::endl;
    // std::cout << "total is " << triangles.size() * 2 << std::endl;
    //td::cout << "o.size: " << Sav << std::endl;
    //std::cout << "1.size: " << Sav2 << std::endl;
    if (Cvp > float(triangles.size() * 2)) {
        // std::cout << "xxxxxxxxxxxxxxxxxxxxxx" << std::endl;
        Node* leafNode = new Node();
        leafNode->myB = bbox;
        leafNode->leftChild = nullptr; //=======
        leafNode->rightChild = nullptr; //=======
        leafNode->isLeaf = true;
        leafNode->triangles = triangles;
        return leafNode;
    }
    
    //std::cout << "yyyyyyyyyyyyyyyyyyyyyyyyy" << std::endl;
    Node* node = new Node();
    node->myB = bbox;
    //std::cout<< triangleList[0].size()
    node->leftChild = nodeBuild(newLeft, triangleList[0], depth + 1);
    node->rightChild = nodeBuild(newRight, triangleList[1], depth + 1);
    //node->isLeaf = (node->leftChild == nullptr && node->rightChild == nullptr);
    return node;
}

void Accel::build() {
    /* Nothing to do here for now */

    std::vector<uint32_t> triangleLists;

    for (uint32_t idx = 0; idx < m_mesh->getTriangleCount(); ++idx) {
        triangleLists.push_back(idx);
    }
    this->root = this->nodeBuild(m_bbox, triangleLists, 0);
    if(this->root->rightChild == nullptr && this->root->leftChild == nullptr){
        this->root->isLeaf = true;
        this->root->triangles = triangleLists;
    }
}

// struct cmp {
//     bool operator()(std::pair<int, float>a, std::pair<int, float>b){
//         return a.second > b.second;
//     }
// };

bool cmp(std::pair<int, float>a, std::pair<int, float>b) {
    return a.second > b.second;
};

// struct sortedNode{
// };

bool Accel::rayIntersect(const Ray3f& ray_, Intersection& its, bool shadowRay) const {
    bool foundIntersection = false;  // Was an intersection found so far?
    uint32_t f = (uint32_t)-1;      // Triangle index of the closest intersection

    Ray3f ray(ray_); /// Make a copy of the ray (we will need to update its '.maxt' value)

    std::stack<std::pair<Node*, int>> dfsStack;

    dfsStack.push(std::make_pair(this->root, 0));


    while(!dfsStack.empty()){

        Node* current = dfsStack.top().first;

        dfsStack.pop();

        // bool inter = false;
        // std::cout<<"1"<<std::endl;
        if(current->myB.rayIntersect(ray)){
            if (current->isLeaf) {  
                for (uint32_t idx = 0; idx < current->triangles.size(); ++idx) {
                    float u, v, t;

                    if (m_mesh->rayIntersect(current->triangles[idx], ray, u, v, t)) {
                        /* An intersection was found! Can terminate
                        immediately if this is a shadow ray query */
                        if (shadowRay){
                            return true;
                        }   
                        ray.maxt = its.t = t;
                        its.uv = Point2f(u, v);
                        its.mesh = m_mesh;
                        f = current->triangles[idx];
                        foundIntersection = true;

                        while(!dfsStack.empty() && dfsStack.top().second != 0){
                            dfsStack.pop();
                        }
                    }
                }
            }

            if (!current->isLeaf) { 

                std::pair<int, float> temp[2];
                std::pair<int, float> first = std::make_pair(-1, (float)999999999);
                for(int i = 0; i < 2; i++){
                    temp[i] = first;
                }
                if(current->leftChild != nullptr){
                    float nearT, farT; 
                    if(current->leftChild->myB.rayIntersect(ray, nearT, farT)){
                        // std::cout<<"Child "<< i << " NearT = "<< nearT<<std::endl;
                        temp[0] = std::make_pair(0, nearT);  
                    }
                }
                if(current->rightChild != nullptr){
                    float nearT, farT; 
                    if(current->rightChild->myB.rayIntersect(ray, nearT, farT)){
                        // std::cout<<"Child "<< i << " NearT = "<< nearT<<std::endl;
                        temp[1] = std::make_pair(1, nearT);  
                    }
                }

                std::sort(temp, temp + 2, cmp);
                // int index =0;
                for(int i = 0; i < 2; i++){
                    float nearT, farT;
                    if(temp[i].first == 0 && current->leftChild->myB.rayIntersect(ray, nearT, farT)){
                        dfsStack.push(std::make_pair(current->leftChild, i));
                        // dfsStack.push(std::make_pair(current->children[temp[i].first], index));
                        // index++;
                    }
                    if(temp[i].first == 1 && current->rightChild->myB.rayIntersect(ray, nearT, farT))
                        dfsStack.push(std::make_pair(current->rightChild, i));
                }
            }
        }
        // ============================================================================================================================================
    }


    if (foundIntersection) {
        /* At this point, we now know that there is an intersection,
           and we know the triangle index of the closest such intersection.

           The following computes a number of additional properties which
           characterize the intersection (normals, texture coordinates, etc..)
        */

        /* Find the barycentric coordinates */
        Vector3f bary;
        bary << 1 - its.uv.sum(), its.uv;

        /* References to all relevant mesh buffers */
        const Mesh* mesh = its.mesh;
        const MatrixXf& V = mesh->getVertexPositions();
        const MatrixXf& N = mesh->getVertexNormals();
        const MatrixXf& UV = mesh->getVertexTexCoords();
        const MatrixXu& F = mesh->getIndices();

        /* Vertex indices of the triangle */
        uint32_t idx0 = F(0, f), idx1 = F(1, f), idx2 = F(2, f);

        Point3f p0 = V.col(idx0), p1 = V.col(idx1), p2 = V.col(idx2);

        /* Compute the intersection positon accurately
           using barycentric coordinates */
        its.p = bary.x() * p0 + bary.y() * p1 + bary.z() * p2;

        /* Compute proper texture coordinates if provided by the mesh */
        if (UV.size() > 0)
            its.uv = bary.x() * UV.col(idx0) +
            bary.y() * UV.col(idx1) +
            bary.z() * UV.col(idx2);

        /* Compute the geometry frame */
        its.geoFrame = Frame((p1 - p0).cross(p2 - p0).normalized());

        if (N.size() > 0) {
            /* Compute the shading frame. Note that for simplicity,
               the current implementation doesn't attempt to provide
               tangents that are continuous across the surface. That
               means that this code will need to be modified to be able
               use anisotropic BRDFs, which need tangent continuity */

            its.shFrame = Frame(
                (bary.x() * N.col(idx0) +
                    bary.y() * N.col(idx1) +
                    bary.z() * N.col(idx2)).normalized());
        }
        else {
            its.shFrame = its.geoFrame;
        }
    }

    return foundIntersection;
}

NORI_NAMESPACE_END

