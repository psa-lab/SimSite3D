#ifndef FACE_BINS_HEADER_INCLUDED
#define FACE_BINS_HEADER_INCLUDED

#include <FaceAttrib.H>

namespace ASCbase{
namespace geometry{

//! Make this a simple class to get it workiing -- then template it for generic
//! associated datatype(s) -- keep all functions here for now because if 
//! we use templates gcc requires them in header anyhow
class FaceBins{
public:
  typedef std::vector< std::vector<FaceAttrib::vci> > bin_vec;
  typedef bin_vec::iterator bin_vi;
  typedef bin_vec::const_iterator bin_vci;

  // The code would probably be faster if we changed this to 
  // std::map<bin_vci, std::vector<const my_float*> >
  typedef std::multimap<bin_vci, const my_float_t*> bin2vert_mmap;
  typedef bin2vert_mmap::const_iterator bin2vert_mmci;

  //! Default constructor -- sometimes the array cannot be constructed 
  //! at class construction 
  FaceBins()
  { init(); }

  FaceBins(VertAttrib::vci verts_begin, VertAttrib::vci verts_end,
           FaceAttrib::vci faces_begin, FaceAttrib::vci faces_end,
           const my_float_t bin_width, const my_float_t max_dist)
  {
    init();
    setup_grid(verts_begin, verts_end, faces_begin, faces_end, bin_width, 
               max_dist);
  }
  
  ~FaceBins() {;}

  void
  setup_grid(VertAttrib::vci verts_begin, VertAttrib::vci verts_end,
             FaceAttrib::vci faces_begin, FaceAttrib::vci faces_end,
             const my_float_t bin_width, const my_float_t max_dist)
  {
    if(A_bins_are_setup){
      std::cerr << "Bins are already setup - cannot set them up again\n";
      return;
    }
    A_bin_width = bin_width;
    A_bin_width_inv = 1.0 / bin_width;
    A_max_dist = max_dist;

    my_float_t min_pos[3], max_pos[3];
    get_min_max_coords(verts_begin, verts_end, min_pos, max_pos);
    determine_space(min_pos, max_pos, A_lower_corner, A_upper_corner, 
                    A_num_bins);
    populate_grid(faces_begin, faces_end, A_max_dist, &A_bins);
    A_bins_are_setup = true;
  }

  void 
  get_min_max_coords(VertAttrib::vci verts_begin, VertAttrib::vci verts_end,
                     my_float_t *min_pos, my_float_t *max_pos)
  {
    std::copy(verts_begin->pos, verts_begin->pos + 3, min_pos);
    std::copy(verts_begin->pos, verts_begin->pos + 3, max_pos);
    for(VertAttrib::vci V = verts_begin; V < verts_end; ++V)
      for(int i = 0; i < 3; ++i){
        if(max_pos[i] < V->pos[i]) max_pos[i] = V->pos[i];
        if(V->pos[i] < min_pos[i]) min_pos[i] = V->pos[i];
      }
  }

  void
  determine_space(const my_float_t *min_pos, const my_float_t *max_pos,
                  my_float_t *min_corner, my_float_t *max_corner, int *num_bins)
  {
    for(int i = 0; i < 3; ++i){
      min_corner[i] = 
        A_bin_width * std::floor(min_pos[i] / A_bin_width) - 2*A_bin_width;
      max_corner[i] = 
        A_bin_width * std::floor(max_pos[i] / A_bin_width) + 2*A_bin_width;
      num_bins[i] = 
        static_cast<int>((max_corner[i] - min_corner[i]) / A_bin_width);
    }
  }

  void
  populate_grid(FaceAttrib::vci faces_begin, FaceAttrib::vci faces_end,
                const my_float_t max_dist, bin_vec *bins)
  {
    // First dimension
    bin_vec dim1_bins;
    dim1_bins.resize(A_num_bins[0]);
    for(FaceAttrib::vci F = faces_begin; F != faces_end; ++F){
      VertAttrib::vi V[3];
      F->get_vertices(V);
      my_float_t min_pt = A_lower_corner[0] - max_dist; 
      my_float_t max_pt = A_lower_corner[0] + A_bin_width + max_dist; 
      for(bin_vi bin = dim1_bins.begin(); bin < dim1_bins.end(); ++bin){
        for(int zz = 0; zz < 3; ++zz){
          if(min_pt <= V[zz]->pos[0] && V[zz]->pos[0] <= max_pt){
            bin->push_back(F);
            break;
          }
        }
        min_pt += A_bin_width; 
        max_pt += A_bin_width; 
      }
    }

    // Second dimension
    bin_vec dim2_bins;
    subdivide_bins(dim1_bins, &dim2_bins, max_dist, 1);

    // Third dimension
    subdivide_bins(dim2_bins, &A_bins, max_dist, 2);
  }

  bin_vci
  get_bin(const my_float_t *pt) const
  {
    // First check if point is outside of bin ranges
    for(int i = 0; i < 3; ++i)
      if(pt[i] < A_lower_corner[i] || A_upper_corner[i] < pt[i])
        return A_bins.end();

    int idx[3];
    for(int i = 0; i < 3; ++i)
      idx[i] = static_cast<int>(std::floor((pt[i] - A_lower_corner[i]) * 
                                           A_bin_width_inv)); 
    return A_bins.begin() + (idx[0] * A_num_bins[1] * A_num_bins[2] +
                             idx[1] * A_num_bins[2] + idx[2]);
  }

  void
  get_bins(const my_float_t *vert_begin, const uint num_vert, 
           bin2vert_mmap *pts_in_bins) const
  {
    const my_float_t *vert_end = vert_begin + 3*num_vert;
    for(const my_float_t *v = vert_begin; v < vert_end; v += 3){
      pts_in_bins->insert(std::pair<bin_vci, const my_float_t*>(get_bin(v), v));
    }  
  }

  //! Constant iterator to first bin in the vector
  inline bin_vci
  bins_begin() const
  { return A_bins.begin(); }

  //! Constant iterator to one past the last bin in the vector
  inline bin_vci
  bins_end() const
  { return A_bins.end(); }
                

private:
  void
  init()
  {
    A_max_dist = my_float_max;
    A_bin_width = my_float_max;
    std::fill(A_upper_corner, A_upper_corner + 3, 0.0);
    std::fill(A_lower_corner, A_lower_corner + 3, 0.0);
    std::fill(A_num_bins, A_num_bins + 3, 0);
    A_bins_are_setup = false;
  }

  void
  subdivide_bins(const bin_vec &bins_in, bin_vec *bins_out, 
                 const my_float_t max_dist, const size_t dim)
  {
    typedef std::vector< FaceAttrib::vci >::const_iterator face_vec_vci;
    bins_out->resize(bins_in.size() * A_num_bins[dim]);

    // Loop over the "large" bins
    bin_vi S_begin = bins_out->begin();
    for(bin_vci L_bin = bins_in.begin(); L_bin < bins_in.end(); ++L_bin){
      for(face_vec_vci F_i = L_bin->begin(); F_i < L_bin->end(); ++F_i){
        VertAttrib::vi V[3];
        (*F_i)->get_vertices(V);
        my_float_t min_pt = A_lower_corner[dim] - max_dist; 
        my_float_t max_pt = A_lower_corner[dim] + A_bin_width + max_dist;  
        // Store in the correct "small" bins
        for(bin_vi S_bin = S_begin; S_bin < S_begin + A_num_bins[dim]; ++S_bin){
          for(int zz = 0; zz < 3; ++zz){
            if(min_pt <= V[zz]->pos[dim] && V[zz]->pos[dim] <= max_pt){
              S_bin->push_back(*F_i);
              break;
            }
          }
          min_pt += A_bin_width; 
          max_pt += A_bin_width; 
        } 
      }
      S_begin += A_num_bins[dim];
    }
  }

  my_float_t A_max_dist;  //!< Maximum distance for grid to be valid
  my_float_t A_bin_width; //!< Side length of a bin
  my_float_t A_bin_width_inv; //!< multiplicative inverse (1 over) side length
  my_float_t A_upper_corner[3];
  my_float_t A_lower_corner[3];
  int A_num_bins[3];
  bin_vec A_bins;
//  bin_vec A_bins_for_moved_pts;
  bool A_bins_are_setup;
  //std::vector<my_float_t *> A_bins;
};
} // end namespace geometry
} // end namespace ASCbase
#endif
