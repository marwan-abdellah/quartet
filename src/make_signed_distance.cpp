#include "make_signed_distance.h"
#include "geometry_queries.h"
#include <omp.h>

// Comment out to use a faster, less accurate point-tri distance measure.
// #define EXACT_POINT_TRI_DISTANCE


static void
check_neighbour(const std::vector<Vec3i> &tri,
                const std::vector<Vec3f> &x,
                Array3f &phi,
                Array3i &closest_tri,
                const Vec3f &gx,
                int i0, int j0, int k0,
                int i1, int j1, int k1)
{
    if(closest_tri(i1,j1,k1) >= 0) {
        int p, q, r; assign(tri[closest_tri(i1,j1,k1)], p, q, r);
	
#ifdef EXACT_POINT_TRI_DISTANCE
		// This is the more accurate, exact distance-to-tri routine
        float d = point_triangle_distance(gx, x[p], x[q], x[r]);
#else
		// This is a fast approximation, using the triangle centroid
        float d = dist(gx, 0.3333333333333333333f*(x[p]+x[q]+x[r]));
#endif
		
        if(d < phi(i0,j0,k0)) {
            phi(i0,j0,k0) = d;
            closest_tri(i0,j0,k0) = closest_tri(i1,j1,k1);
        }
    }
}


static void
sweep(const std::vector<Vec3i> &tri,
      const std::vector<Vec3f> &x,
      SDF &sdf,
      Array3i &closest_tri,
      int di, int dj, int dk)
{
    // figure out loop bounds from direction (di,dj,dk)
    int i0, i1;
    if (di>0) { i0=1; i1=sdf.phi.ni; }
    else { i0=sdf.phi.ni-2; i1=-1; }
    int j0, j1;
    if (dj>0) { j0=1; j1=sdf.phi.nj; }
    else { j0=sdf.phi.nj-2; j1=-1; }
    int k0, k1;
    if (dk>0) { k0=1; k1=sdf.phi.nk; }
    else { k0=sdf.phi.nk-2; k1=-1; }

    // sweep through the grid propagating closest triangle information
    int i, j, k;
    for(k=k0; k!=k1; k+=dk)
        for(j=j0; j!=j1; j+=dj)
            for(i=i0; i!=i1; i+=di)
            {
                Vec3f gx(i*sdf.dx+sdf.origin[0],
                        j*sdf.dx+sdf.origin[1],
                        k*sdf.dx+sdf.origin[2]);

                check_neighbour(tri,x, sdf.phi, closest_tri, gx, i,j,k, i-di,j,   k);
                check_neighbour(tri,x, sdf.phi, closest_tri, gx, i,j,k, i,   j-dj,k);
                check_neighbour(tri,x, sdf.phi, closest_tri, gx, i,j,k, i-di,j-dj,k);
                check_neighbour(tri,x, sdf.phi, closest_tri, gx, i,j,k, i,   j,   k-dk);
                check_neighbour(tri,x, sdf.phi, closest_tri, gx, i,j,k, i-di,j,   k-dk);
                check_neighbour(tri,x, sdf.phi, closest_tri, gx, i,j,k, i,   j-dj,k-dk);
                check_neighbour(tri,x, sdf.phi, closest_tri, gx, i,j,k, i-di,j-dj,k-dk);
            }
}


void
calculateSignedDistanceFunction(const std::vector<Vec3i> &tri,
                     const std::vector<Vec3f> &x,
                     SDF &sdf)
{
    // Prepare an array to store (or approximate) the closest triangle to
    // the grid points
    Array3i closestTriangle(sdf.phi.ni, sdf.phi.nj, sdf.phi.nk, -1);

    // Prepare an array to count the intersections of triangles of the
    // iso-surface with grid edges
    // intersectionCount(i,j,k) is the number of tri intersections in
    // (i-1,i]x{j}x{k}
    Array3i intersectionCount(sdf.phi.ni, sdf.phi.nj, sdf.phi.nk, 0);

    // Pre-compute double-precision dx and 1/dx (or the multiplicative inverse)
    double dxDouble = static_cast< double >(sdf.dx);
    double dxDoubleInverse = 1.0 / dxDouble;
    float dxInverse = static_cast< float >(dxDoubleInverse);

    // Pre-compute double-precision origin
    double oX = static_cast< double >(sdf.origin[0]);
    double oY = static_cast< double >(sdf.origin[1]);
    double oZ = static_cast< double >(sdf.origin[2]);

    // initialize correct distances near the mesh and count intersections
    std::cout<<"  rasterizing triangles"<<std::endl;
    int i, j, k;

    #pragma omp parallel for
    for(size_t t = 0; t < tri.size( ); ++t)
    {
        int p, q, r;
        assign(tri[t], p, q, r);

        // get padded bounding box of triangle in grid coordinates
        int i0=(int)((min(x[p][0],x[q][0],x[r][0]) - oX)*dxDoubleInverse - 0.5),
            i1=(int)((max(x[p][0],x[q][0],x[r][0]) - oX)*dxDoubleInverse + 1.5),
            j0=(int)((min(x[p][1],x[q][1],x[r][1]) - oY)*dxDoubleInverse - 0.5),
            j1=(int)((max(x[p][1],x[q][1],x[r][1]) - oY)*dxDoubleInverse + 1.5),
            k0=(int)((min(x[p][2],x[q][2],x[r][2]) - oZ)*dxDoubleInverse - 0.5),
            k1=(int)((max(x[p][2],x[q][2],x[r][2]) - oZ)*dxDoubleInverse + 1.5);

        // clamp appropriately to lie within the grid
        if (i0<0) i0=0; else if (i0>sdf.phi.ni-1) i0=sdf.phi.ni-1;
        if (i1<0) i1=0; else if (i1>sdf.phi.ni-1) i1=sdf.phi.ni-1;
        if (j0<0) j0=0; else if (j0>sdf.phi.nj-1) j0=sdf.phi.nj-1;
        if (j1<0) j1=0; else if (j1>sdf.phi.nj-1) j1=sdf.phi.nj-1;
        if (k0<0) k0=0; else if (k0>sdf.phi.nk-1) k0=sdf.phi.nk-1;
        if (k1<0) k1=0; else if (k1>sdf.phi.nk-1) k1=sdf.phi.nk-1;

        // do distances and intersection counts in this chunk
        for(k=k0; k<=k1; ++k) for(j=j0; j<=j1; ++j) for(i=i0; i<=i1; ++i) {
            Vec3f gx(i*sdf.dx+sdf.origin[0], 
                     j*sdf.dx+sdf.origin[1], 
                     k*sdf.dx+sdf.origin[2]);
            float d = point_triangle_distance(gx, x[p], x[q], x[r]);
            if (d < sdf.phi(i,j,k)) {
                sdf.phi(i,j,k) = d;
                closestTriangle(i,j,k) = t;
            }
        }

        // and do intersection counts
        for(int k=k0; k<=k1; ++k)
            for(int j=j0; j<=j1; ++j)
            {
                double gy=j*dxDouble+sdf.origin[1], gz=k*dxDouble+sdf.origin[2];
                double a, b, c;
                if (point_in_triangle_2d(gy, gz,
                                         x[p][1], x[p][2],
                                         x[q][1], x[q][2],
                                         x[r][1], x[r][2], a, b, c))
                {

                    // estimate x coordinate of intersection
                    double xint = a*x[p][0] + b*x[q][0] + c*x[r][0];

                    // convert to a grid coordinate
                    int i_interval=(int)(std::ceil((xint-sdf.origin[0])*dxInverse));

                    // intersection is in (i_interval-1,i_interval]
                    if (i_interval < 0)
                        ++intersectionCount(0, j, k);
                    else if (i_interval < sdf.phi.ni)
                        ++intersectionCount(i_interval, j, k);
                    // ignore intersections beyond the +x side of the grid
                }
        }
    }

    // fill in the rest of the distances with fast sweeping
    // (a single sweep is adequate to get rough approximation everywhere)
    std::cout<<"  sweeping to rest of grid"<<std::endl;
    for (unsigned int pass=0; pass<1; ++pass){
        sweep(tri, x, sdf, closestTriangle, +1, +1, +1);
        sweep(tri, x, sdf, closestTriangle, -1, -1, -1);
        sweep(tri, x, sdf, closestTriangle, +1, +1, -1);
        sweep(tri, x, sdf, closestTriangle, -1, -1, +1);
        sweep(tri, x, sdf, closestTriangle, +1, -1, +1);
        sweep(tri, x, sdf, closestTriangle, -1, +1, -1);
        sweep(tri, x, sdf, closestTriangle, +1, -1, -1);
        sweep(tri, x, sdf, closestTriangle, -1, +1, +1);
    }

    // figure out signs (inside/outside) from parity of intersection counts
    std::cout<<"  calculating signs"<<std::endl;
    for(int k=0; k<sdf.phi.nk; ++k)
    {
        for(int j=0; j<sdf.phi.nj; ++j)
        {
            int total_count = 0;
            for (int i=0; i<sdf.phi.ni; ++i)
            {
                total_count += intersectionCount(i,j,k);

                if (total_count%2 == 1)
                {
                    // if parity of intersections so far is odd,
                    sdf.phi(i,j,k) = -sdf.phi(i,j,k); // we are inside the mesh
                }
            }
        }
    }
}

