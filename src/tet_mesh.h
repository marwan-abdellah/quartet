#ifndef TET_MESH_H
#define TET_MESH_H

#include "vec.h"
#include <vector>

/**
 * @brief The Tet struct
 * Tetrahedron
 */
struct Tet
{
    /**
     * @brief Tet
     */
    Tet();

    /**
     * @brief Tet
     * @param v0
     * @param v1
     * @param v2
     * @param v3
     */
    Tet(const Vec3f& v0, const Vec3f& v1, const Vec3f& v2, const Vec3f& v3);

    /**
     * @brief computeOrientation
     * @param sixSignedVolume
     * @return
     */
    int computeOrientation(double &sixSignedVolume) const;

    /**
     * @brief computeOrientation
     * @return
     */
    int computeOrientation() const;

    /**
     * @brief computeVolume
     * @return
     */
    float computeVolume() const;

    /**
     * @brief computeAspectRatio
     * @return
     */
    float computeAspectRatio() const;

    /**
     * @brief computeDihedralAngles
     * @param angles
     */
    void computeDihedralAngles(std::vector<float>& angles) const;

    /**
     * @brief computeFaceAngles
     * @param angles
     */
    void computeFaceAngles(std::vector<float>& angles) const;

    /**
     * @brief computeMaxEdgeLength
     * @return
     */
    float computeMaxEdgeLength() const;

public:

    /**
     * @brief operator []
     * @param i
     * @return
     */
    Vec3f& operator[] (int i)
    {
        return v[i];
    }

    /**
     * @brief operator []
     * @param i
     * @return
     */
    const Vec3f& operator[] (int i) const
    {
        return v[i];
    }

public:

    /**
     * @brief v
     * The four vertices of the tetrahedron.
     */
    Vec3f v[4];
};

/**
 * @brief The TetMesh struct
 * Data structure representing a tetrahedral mesh.
 */
struct TetMesh
{
public:

    /**
     * @brief TetMesh
     * Default constructor
     */
    TetMesh()
        : _dirty(false)
    {
        /// EMPTY CONSTRUCTOR
    }

    /**
     * @brief TetMesh
     * Copy constructor
     * @param mesh
     * Input tetrahedral mesh.
     */
    TetMesh(const TetMesh& mesh)
        : _vertices(mesh._vertices)
        , _tets(mesh._tets)
        , _dirty(true)
    {
        /// EMPTY CONSTRUCTOR
    }

    /**
     * @brief TetMesh
     * Construct from vertices and tet indices.
     * @param verts
     * @param tets
     */
    TetMesh(const std::vector<Vec3f>& verts, 
            const std::vector<Vec4i>& tets)
        : _vertices(verts)
        , _tets(tets)
        , _dirty(true)
    {
         /// EMPTY CONSTRUCTOR
    }

    /**
     * @brief getVertices
     * Returns the position of the vertex stored at index i.
     * @param i
     * @return
     */
    Vec3f& getVertices(int i)
    {
        return _vertices[i];
    }

    /**
     * @brief getVertices
     * Returns the position of the vertex stored at index i.
     * @param i
     * @return
     */
    const Vec3f& getVertices(int i) const
    {
        return _vertices[i];
    }

    /**
     * @brief getTets
     * Returns the indices of the tetrahedron stored at index i.
     * @param i
     * @return
     */
    Vec4i& getTets(int i)
    {
        _dirty = true;
        return _tets[i];
    }

    /**
     * @brief getTets
     * @param i
     * @return
     */
    const Vec4i& getTets(int i) const
    {
        return _tets[i];
    }

    /**
     * @brief getListVertices
     * Returns the list of this tet mesh's vertices.
     * @return
     */
    std::vector<Vec3f>& getListVertices()
    {
        _dirty = true;
        return _vertices;
    }

    /**
     * @brief getListVertices
     * Returns the list of this tet mesh's vertices.
     * @return
     */
    const std::vector<Vec3f>& getListVertices() const
    {
        return _vertices;
    }

    /**
     * @brief getListTets
     * Returns the list of this tet mesh's tetrahedron indices.
     * @return
     */
    std::vector< Vec4i >& getListTets()
    {
        _dirty = true;
        return _tets;
    }

    /**
     * @brief getListTets
     * Returns the list of this tet mesh's tetrahedron indices.
     * @return
     */
    const std::vector< Vec4i >& getListTets() const
    {
        return _tets;
    }

    /**
     * @brief getNumberVertices
     * Returns the number of vertices in the tet mesh.
     * @return
     */
    uint64_t getNumberVertices() const
    {
        return _vertices.size();
    }

    /**
     * @brief getNumberTets
     * Returns the number of tetrahedra in the mesh.
     * @return
     */
    uint64_t getNumberTets() const
    {
        return _tets.size();
    }

    /**
     * @brief getTet
     * Returns the tetrahedron stored at index t.
     * @param i
     * @return
     */
    Tet getTet(uint64_t i) const;

    /**
     * @brief getAdjacentVertices
     * Fills the given list with all indices of all the vertices that are
     * adjacent to the given vertex.
     * @param vertex
     * @param _vertices
     */
    void getAdjacentVertices(int vertex, std::vector<int>& _vertices) const;

    /**
     * @brief getIncidentTets
     * Fills the given list with all the tets incident to the given vertex.
     * @param vertex
     * @param _tets
     */
    void getIncidentTets(int vertex, std::vector<Tet>& getListTets) const;

    /**
     * @brief getIncidentTetIndices
     * Fills the given list with the indices of tets incident to vertex.
     * @param vertex
     * @param tetIndices
     */
    void getIncidentTetIndices(int vertex, std::vector<int>& tetIndices) const;

    /**
     * @brief compactMesh
     * Clean up the mesh, getting rid of unused vertices.
     * Tetrahedral indices adjusted accordingly.
     */
    void compactMesh();

    /**
     * @brief getBoundary
     * Extract boundary vertices and triangles from the tet mesh.
     * Vertices will be returned in boundaryVerts as indices into the
     * tet mesh vertex list _vertices.
     *  Triangles are stored in boundaryTris.
     * These are also indices into _vertices, NOT as indices into boundaryVerts.
     * Assumes that the tet mesh is well-formed and closed.
     * @param boundaryVerts
     * @param boundaryTris
     */
    void getBoundary(std::vector<int>& boundaryVerts,
                     std::vector<Vec3i>& boundaryTris) const;

    /**
     * @brief getBoundary
     * Extract boundary vertices and triangles from the tet mesh.
     * Vertices will be copied into boundaryVerts.
     * Triangles are stored in boundaryTris.
     * These are indices into boundaryVerts.
     * Assumes that the tet mesh is well-formed and closed.
     * @param boundaryVerts
     * @param boundaryTris
     */
    void getBoundary(std::vector<Vec3f>& boundaryVerts,
                     std::vector<Vec3i>& boundaryTris) const;

    // Write this TetMesh to file with the given name.
    // The output file will be of the .TET format.
    bool writeToFile(const char* filename) const;
    bool writeNodeFile(const char* filename) const;
    bool writeEdgeFile(const char* filename) const;
    bool writeEleFile(const char* filename) const;
    bool writeFaceFile(const char* filename) const;
    bool writeSmeshFile(const char* filename) const;
    bool writeToGMshFile(const char* filename) const;

    // Compute information about this TetMesh (dihedral angles, tet volumes,
    // etc.) and write it to the file with the given name.
    bool writeInfoToFile(const char* filename) const;

private:

    /**
     * @brief _buildIncidenceMap
     * Re-create the incidence map from the current _vertices and _tets lists.
     * @note This is const because we need to be able to call it from
     * getIncidentTets(), which is const.
     * _incidenceMap and _dirty are mutable as a result.
     */
    void _buildIncidenceMap() const;

private:

    /**
     * @brief _vertices
     * A list of all the vertices of the tetrahedral mesh.
     */
    std::vector< Vec3f > _vertices;

    /**
     * @brief _tets
     * A list of all the tets of the tetrahedral mesh.
     */
    std::vector< Vec4i > _tets;

    /**
     * @brief _incidenceMap
     * Stores the indices of the tetrahedra incident to each vertex.
     */
    mutable std::vector< std::vector< int > > _incidenceMap;

    /**
     * @brief _dirty
     * Indicates whether the incidence map may need to be rebuilt.
     */
    mutable bool _dirty;
};

#endif
