#ifndef SDF_H
#define SDF_H

#include "array3.h"
#include "vec.h"

/**
 * @brief The SDF struct
 * The Signed Distance Function (or Field) struct.
 * The SDF class encapsulates a scalar field that stores the shortest distance
 * to an implicit surface, that is therefore found at the zero-isosurface.
 * The field has a negative sign inside the surface and a positive sign outside.
 * The scalar field phi is stored as a three-dimensional array of values, and
 * trilinear interpolation is provided for arbitrary spatial lookups.
 */
struct SDF
{
    Array3f phi;
    const Vec3f origin;
    const float dx;

    // Constructor
    // origin is the origin of the 3D array.
    // dx is the grid spacing in 3D space
    // ni*nj*nk are the dimensions of the grid.
    SDF(const Vec3f& origin,
        const float dx,
        int ni, int nj, int nk);

    /**
     * @brief operator ()
     * Evaluates the scalar field at the given point in the space p using
     * tri-linear interpolation.
     * @param p
     * A three-dimensional point the space.
     * @return
     * The SDF value at the given point.
     */
    float operator() (const Vec3f& p) const;

    /**
     * @brief computeGradient
     * Compute the gradient vector of the field at the given point p.
     * @param p
     * A three-dimensional point the space.
     * @return
     * The gradient vector of the field at the given point p.
     */
    Vec3f computeGradient(const Vec3f& p) const;

    // Compute the "normal" vector of the field at the given point.
    // This is equivalent to the normalized gradient vector.

    /**
     * @brief computeNormal
     * Compute the "normal" vector of the field at the given point p.
     * @param p
     * A three-dimensional point the space.
     * @return
     * The the "normal" vector of the field at the given point p.
     */
    Vec3f computeNormal(const Vec3f& p) const;

    /**
     * @brief projectToIsosurface
     * Project the given point x onto the zero-isosurface of the field.
     * @param p
     * @note
     * Note that this is only accurate if x is already near the isosurface.
     * @return
     * The given point x onto the zero-isosurface of the field.
     */
    Vec3f projectToIsosurface(const Vec3f& p) const;

private:
    const float over_dx; // 1/dx
};

#endif
