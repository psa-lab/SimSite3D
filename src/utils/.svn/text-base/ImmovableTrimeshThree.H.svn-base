#ifndef IMMOVABLE_TRIMESH_THREE_HEADER_INCLUDED
#define IMMOVABLE_TRIMESH_THREE_HEADER_INCLUDED

#include <SimpleTrimeshThree.H>
#include <FaceBins.H>

namespace SimSite3D{
namespace geometry{

#if 0
class ImmovableFace : public FaceAttrib{

  //! Default constructor
  ImmovableFace() : FaceAttrib()
  {
  }

  //! Construct using 3 vertices
  /*!
   * Assumption: The vertices are in counter clockwise order so that 
   * barycentric method to determine if the projection of a point lies inside 
   * the face places the centroid of the face inside.  Besides, the right
   * hand rule has been typically followed up to this point
   *
   * NEED TO THINK ABOUT THIS MORE -- after rotation we could mess up 
   * the cross product?  If it flips to be clockwise with respect to
   * the right hand rule we should be ok since the resulting vector of the
   * cross products will still face the same direction for all 3 vectors
   */
  ImmovableFace(vert_attrib_vi V0, vert_attrib_vi V1, vert_attrib_vi V2)
    : FaceAttrib(V0, V1, V2)
  {
    A_V[2] = V2;
    my_float_t U[3], V[3];
    vector(3, A_V[0]->pos, A_V[1]->pos, U);
    vector(3, A_V[2]->pos, A_V[1]->pos, V);
    cross(U, V, A_normal);
    A_area = normalize(A_normal);
  }

  //! Copy constructor
  ImmovableFace(const ImmovableFace& F) : FaceAttrib(F)
  {
    do_copy(F);
  }

  //! Assignment operator
  const ImmovableFace&
  operator=(const ImmovableFace& F)
  {
    if(&F != this){
      do_copy(F);
      FaceAttrib::operator= (F);
    }
    return *this;
  }

  virtual ~ImmovableFace()
  {
  }



};
#endif


class ImmovableTrimeshThree : public SimpleTrimeshThree{
public:
  //! Default constructor
  ImmovableTrimeshThree(const std::string fname,
                        const my_float_t max_pt2surf_dist = 1.5);

  //! Do nothing explicit dstr
  ~ImmovableTrimeshThree()
  { ; } 

  bool
  compare(const my_float_t *other_vert_begin, const size_t num_other_vert, 
          my_float_t *closest_pts, my_float_t *dists, my_float_t *Ns_for_cp = 0,
          const my_float_t max_dist=1.5) const;

  FaceBins::bin_vci
  get_bin(const my_float_t *pt) const
  { return A_grid.get_bin(pt); }

private:
  FaceBins A_grid;
  my_float_t A_max_pt2surf_dist;  //!< Maximum distance for pt 2 surface distance to be reliable

};

} // end namespace geometry
} // end namespace SimSite3D
#endif
