#ifndef _LIBUMESH_H_INCLUDED_
#define _LIBUMESH_H_INCLUDED_

// LibUMesh: A library for geometry operations on uniform
//           meshes defined by lattices.
//
// We aim to generate the most well-centered meshes with proper
// circumcentric duals and positive-definite Hodge stars. Thus,
// for orthogonal lattices, we compute not simplicial meshes,
// but hexahedral meshes. The meshes are implicitly defined
// by local incidence information for a single unit cell, so
// these meshes are suitable for both truncated and periodic
// domains.
//
// For the uninitiated:
//  We want to make meshes where the the elements are as
// isotropic in shape as possible. This means that triangles
// should have angles as close to 60 degrees as possible and
// tetrahedra should have dihedral angles that are not too
// big or small.
//  We define the dual mesh of the mesh to be the following
// (in 3D):
//   Each tet in the original mesh is a vertex in the dual
//  mesh, whose position is the circumcenter of the tet.
//   Each face in the original mesh is an edge in the dual
//  mesh connecting dual vertices corresponding to adjacent
//  tets.
//   Each edge in the original mesh is a face in the dual mesh
//  perpendicular to the original edge.
//   Each vertex in the original mesh is a polyhedral region
//  around the vertex.
// The fancy name for this is a simplicial complex.

typedef struct UMesh2_struct{
	// The uniform mesh is defined by two lattice vectors,
	// but there may be a third for non-orthogonal lattices.
	// n_edges is the number of edge vectors in one unit cell.
	int n_edges;
	// The edge vectors in a unit cell, in pairs of xy coordinates.
	// The first two are always the lattice vectors, the third being
	// the extra vector for non-orthgonal lattices.
	double Lr[3*2];
	// For an orthogonal lattice, we would have:
	//
	//                  |                 |
	//                  +-----------------+----
	//                 ^|                 |
	//                 ||                 |
	//   (Lr[2],Lr[3]) ||                 |
	//                 ||                 |
	//                  +-----------------+----
	//                   ---------------->
	//                     (Lr[0],Lr[1])
	// For nonorthgonal lattices:
	//
	//
	//                  +-----------+
	//                ^/ \^        /
	// (Lr[2],Lr[3]) //   \\ _____/____ (Lr[4],Lr[5])
	//              //     \\    /
	//             //       \\  /
	//            //         \\/
	//            +-----------+
	//             ---------->
	//             (Lr[0],Lr[1])
	// etc.
	
	
	// inc01 is the incidence matrix between edges and vertices
	//        from to
	//         u v u v
	// u vec [ 0 0 1 0 ]
	// v vec [ 0 0 0 1 ]
	// w vec [ 0 0 1 1 ] w = u+v in this case
	// w could also be [ 1 0 0 1 ], w = from u to v = v-u
	// where u is (Lr[0],Lr[1]), v and w are the next two.
	signed char inc01[3*(2+2)];
	
	// There can be only 1 face for orthogonal lattices.
	// There are 2 faces for non-orthgonal lattices.
	int n_faces;
	
	// inc12 is the incidence matrix between edges and faces
	//    4 x {which u-off v-off}
	// face0 [ 1 0 0   2 1 0   -1 0 1   -2 0 0 ] for a square lattice
	// face1 [ 0 ]
	// or for a non-orthogonal lattice:
	//  [ 1 0 0    2 1 0   -3 1 1   0 0 0 ] for type 1
	//  [ 3 0 0   -1 0 1   -2 0 0   0 0 0 ]
	//signed char inc12[2*4*3];
	unsigned int inc12[2];
#define LibUMesh2_inc12_getedge(UI,EI) (((UI)>>((EI)<<3))&0xFF)
#define LibUMesh2_inc12_which_MASK 0x1C
#define LibUMesh2_inc12_uoff_MASK  0x02
#define LibUMesh2_inc12_voff_MASK  0x01
#define LibUMesh2_inc12_which(UC) ((int)(((UC)&LibUMesh2_inc12_which_MASK)>>2)-3)
#define LibUMesh2_inc12_uoff(UC) ((int)((UC)&LibUMesh2_inc12_uoff_MASK))
#define LibUMesh2_inc12_voff(UC) ((int)((UC)&LibUMesh2_inc12_voff_MASK))
	
	
#define LibUMesh_d_ind_MASK  0x0F
#define LibUMesh_d_uoff_MASK 0x10
#define LibUMesh_d_voff_MASK 0x20
#define LibUMesh_d_woff_MASK 0x40
#define LibUMesh_d_sign_MASK 0x80

	// d10 tells which vertices are incident on each edge.
	// Each row is one edge (there can be up to 3 rows).
	// There are two columns; the first is the from vertex (-1),
	// the second is the to vertex (+1).
	// Each entry stores the element index (for vertices, always zero)
	// in the lowest 4 bits. The next bit give give the u-offset of
	// the vertex, and the next higher bit gives the v-offset.
	// The offsets are encoded as {0,+1} -> {0,1}. The highest
	// bit gives the sign of the vertex (0 = +1, 1 = -1). For d10
	// the first column is all -1, the second all +1.
	unsigned char d10[3*2]; // #edges by 2 matrix, stride 3
	
	// d21 tells which edges are incident on each face.
	unsigned char d21[2*4]; // #faces by #edges/face matrix, stride 2
	unsigned char *d[2];
	unsigned char ld[2];
	unsigned char drows[2];
	unsigned char dcols[2];
	
	// Diagonal hodge stars for vertices, edges, and faces
	double star0, star1[3], star2[2];
	
	// type = 0: orthogonal, n_edges = 2, n_faces = 1
	// type = 1: u.v < 0, n_edges = 3, n_faces = 2
	// type = 2: u.v > 0, n_edges = 3, n_faces = 2
	int type;
} UMesh2;

// Purpose: Fills in the mesh information corresponding to the
//          implicitly defined mesh of a lattice defined by
//          lattice vectors u,v.
int LibUMesh2_Create(const double u[2], const double v[2], UMesh2 *mesh);

// Purpose: Computes the Voronoi neighborhood of an vertex.
// Description:
//   The Voronoi neighborhood of a vertex is the set of points
//  nearer (by Euclidean distance) to the vertex than to any other
//  vertex. It is described by a convex polygonal region. This
//  function returns the polygon in boundary representation as a
//  set of ordered vertex positions.
//   There can be up to 6 vertices involved.
// Arguments:
//  p - OUT - Returns half the vertices of Vornoi polygon, in positive
//            orientation order (right hand rule). Only half are
//            computed since the other half are simply the negations
//            of the returned vertices.
// Returns value:
//  Number of vertices returned (2 or 3).
//  Or, -n if the n-th argument is invalid.
int LubUMesh2_Neighborhood0(const UMesh2 *mesh, double p[6]);

// Purpose: Computes the Voronoi neighborhood of an edge.
// Description:
//   The Voronoi neighborhood of an edge is the convex hull of
//  the edge itself with the line segment connecting the circum-
//  centers of the faces incident to the edge.
// Arguments:
//  which - IN - Index of edge whose neighborhood to compute.
//               Can be 0 to n_edges-1.
//  p    - OUT - Returns the two adjacent circumcenters, ordered so
//               that the vector from the first to the second is in
//               the direction of a 90 degree CW rotation of the
//               edge vector (positive orientation)
// Return value:
//  0 on success.
//  Or, -n if the n-th argument is invalid.
int LibUMesh2_Neighborhood1(const UMesh2 *mesh, int which, double p[4]);


typedef struct UMesh3_struct{
	// There are 3 possible dot products between the lattice vectors:
	//   a.b   b.c   c.a
	// We classify primarily based on which of these products are zero.
	// First, when all 3 are zero (all orthogonal)
	//   type = 0: all vectors orthogonal (cubic, tetragonal, orthorhombic)
	// Next, when 2 are zero (1 special orthogonal axis)
	//   type = 1,2: one vector is orthogonal to the other 2,
	//               same as types 1,2 for 2D meshes (hexagonal, monoclinic)
	// At this point, there are just triclinic and rhombohedral lattices.
	// It is difficult to classify these, so they are all lumped together
	//   type = 3: everything else
	int type;
	
	// 3 lattice vectors, 3 face divisions, and one long diagonal
	double Lr[7*3];
	int n_edges, n_faces, n_tets;
	
	// Let the primary lattice vectors be a,b,c
	// Let the primary face subdivision vectors be d,e,f,
	//  where d is on the face opposite a, joining tips of b and c.
	// Let the diagonal vector be g.
	
	// inc01 is the incidence matrix between vertices and edges
	//      from     to
	//      a b c   a b c
	//  a [ 0 0 0   1     ]
	//  b [ 0 0 0     1   ]
	//  c [ 0 0 0       1 ]
	//  d [               ] (for a cubic lattice there are only 3 edges)
	//  e [               ]
	//  f [               ]
	//  g [               ]
	//signed char inc01[7*(3*3)];
	// We pack the above into bits:
	unsigned char inc01[7];
#define LibUMesh3_inc01_from_a_MASK(UC) 0x20
#define LibUMesh3_inc01_from_b_MASK(UC) 0x10
#define LibUMesh3_inc01_from_c_MASK(UC) 0x08
#define LibUMesh3_inc01_to_a_MASK(UC)   0x04
#define LibUMesh3_inc01_to_b_MASK(UC)   0x02
#define LibUMesh3_inc01_to_c_MASK(UC)   0x01

#define LibUMesh3_inc01_from_a(UC) ((int)((UC)&LibUMesh3_inc01_from_a_MASK))
#define LibUMesh3_inc01_from_b(UC) ((int)((UC)&LibUMesh3_inc01_from_b_MASK))
#define LibUMesh3_inc01_from_c(UC) ((int)((UC)&LibUMesh3_inc01_from_c_MASK))
#define LibUMesh3_inc01_to_a(UC)   ((int)((UC)&LibUMesh3_inc01_to_a_MASK))
#define LibUMesh3_inc01_to_b(UC)   ((int)((UC)&LibUMesh3_inc01_to_b_MASK))
#define LibUMesh3_inc01_to_c(UC)   ((int)((UC)&LibUMesh3_inc01_to_c_MASK))
	
	// inc12 gives incidence information about faces and edges
	//    4 x {which, a-off, b-off, c-off}
	// face0 [   4      0      1      1    ... ]
	// which = 1-based index of which of a,b,c,d,e,f,g (0 for none, negative for negative orientation)
	// a,b,c-off = (0-1) integer multiples of a,b,c to add to the vector
	//signed char inc12[12*4*4];
	// We pack the above into bits: top 5 bits for which+7, bottom 3 for a,b,c-off:
	// lowest order 8 bits for first edge, etc.
	unsigned int inc12[12];
#define LibUMesh3_inc12_getedge(UI,EI) (((UI)>>((EI)<<3))&0xFF)
#define LibUMesh3_inc12_which_MASK 0xF8
#define LibUMesh3_inc12_aoff_MASK  0x04
#define LibUMesh3_inc12_boff_MASK  0x02
#define LibUMesh3_inc12_coff_MASK  0x01
#define LibUMesh3_inc12_which(UC) ((int)(((UC)&LibUMesh3_inc12_which_MASK)>>3)-7)
#define LibUMesh3_inc12_aoff(UC) ((int)((UC)&LibUMesh3_inc12_aoff_MASK))
#define LibUMesh3_inc12_boff(UC) ((int)((UC)&LibUMesh3_inc12_boff_MASK))
#define LibUMesh3_inc12_coff(UC) ((int)((UC)&LibUMesh3_inc12_coff_MASK))
	
	// inc23 gives incidence information about faces and tets/prisms/cubes
	//    6 x {which, a-off, b-off, c-off}
	// tet0 [ ... ]
	// which = 1-based index of which face (0 for none, negative for negative orientation)
	//int inc23[6*6*4];
	// We pack the above into bits: top 5 bits for which+12, bottom 3 for a,b,c-off:
	unsigned char inc23[6*6];
#define LibUMesh3_inc23_which_MASK 0xF8
#define LibUMesh3_inc23_aoff_MASK  0x04
#define LibUMesh3_inc23_boff_MASK  0x02
#define LibUMesh3_inc23_coff_MASK  0x01
#define LibUMesh3_inc23_which(UC) ((int)(((UC)&LibUMesh3_inc23_which_MASK)>>3)-12)
#define LibUMesh3_inc23_aoff(UC) ((int)((UC)&LibUMesh3_inc23_aoff_MASK))
#define LibUMesh3_inc23_boff(UC) ((int)((UC)&LibUMesh3_inc23_boff_MASK))
#define LibUMesh3_inc23_coff(UC) ((int)((UC)&LibUMesh3_inc23_coff_MASK))
	
	// Each incidence relationship is {-1,0,+1}, plus offset
	// information {-1,0,+1} in each direction, so we need 3*2 total
	// bits per incidence matrix entry.
	// For each {-1,0,+1}, we store them as {2,0,1}
	// Lowest 2 bits are the incidence relationship,
	// next 2 are u-offset, next 2 are v-offset
	unsigned char d10[7*2]; // #edges by 2 matrix, stride 7
	unsigned char d21[12*4]; // #faces by #edges/face matrix, stride 12
	unsigned char d32[6*6]; // #tets by #faces/tet matrix, stride 6
	unsigned char *d[3];
	unsigned char ld[3]; // strides
	
	// Diagonal hodge stars for vertices, edges, faces, and tets.
	// Each corresponds to the row of the incidence info.
	double star0, star1[7], star2[12], star3[6];
} UMesh3;


// Purpose: Fills in the mesh information corresponding to the
//          implicitly defined mesh of a lattice defined by
//          lattice vectors a,b,c.
int LibUMesh3_Create(
	const double a[3],
	const double b[3],
	const double c[3],
	UMesh3 *mesh
);

// Purpose: Computes the Voronoi neighborhood of a vertex
// Description:
//   The Voronoi neighborhood of a vertex is the set of points
//  nearer (by Euclidean distance) to the vertex than to any other
//  vertex. It is described by a convex polyhedral region with flat
//  polygonal faces. This function returns the polyhedron in boundary
//  representation as a set of vertex positions and face vertex lists.
//   There can be up to 24? vertices (truncated octahedron) involved.
// Arguments:
//  v - OUT - Returns vertices of polyhedron.
//            Must be length at least 24*3.
//  f - OUT - Returns the indices (into v) of vertices comprising each face.
//            Each face is allocated 6 entries, and there can be up to 14
//            faces, so f should be at least length 6*14. Unused entries
//            are filled with -1.
// Returns value:
//  Number of faces.
//  Or, -n if the n-th argument is invalid.
int LibUMesh3_Neighborhood0(const UMesh3 *mesh, double *v, int *f);

// Purpose: Computes the Voronoi neighborhood of an edge.
// Description:
//   The Voronoi neighborhood of an edge is the union of tets.
//  Each tet is the convex hull of two line segments. One segment
//  is always the edge itself. The other edge connects circumcenters
//  of adjacent tets that are incident to the edge.
// Arguments:
//  which - IN - Index of edge whose neighborhood to compute.
//               Can be 0 to n_edges-1.
//  v    - OUT - Length 3*6. There can be up to 6 incident tets,
//               Hence 6 circumcenters in the 1-ring. The apex
//               vertices of the ends of the edge are implicit
//               in which edge is chosen, and hence not returned.
//               The vertices are returned in positive orientation
//               order (right hand rule).
// Return value:
//  The number of vertices.
//  Or, -n if the n-th argument is invalid.
int LibUMesh3_Neighborhood1(const UMesh3 *mesh, int which, double *v);

// Purpose: Computes the Voronoi neighborhood of a face.
// Description:
//   The Voronoi neighborhood of a face is the convex hull of the
//  face itself with the line segment joining the circumcenters of
//  the tets/prisms/cubes incident to the face.
// Arguments:
//  which - IN - Index of face whose neighborhood to compute.
//               Can be 0 to n_faces-1.
//  v    - OUT - The two vertices (in positive orientation order) of
//               the adjacent circumcenters.
// Return value:
//  0 on success.
//  Or, -n if the n-th argument is invalid.
int LibUMesh3_Neighborhood2(const UMesh3 *mesh, int which, double v[6]);

#endif // _LIBUMESH_H_INCLUDED_
