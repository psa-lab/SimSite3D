#ifndef DISTANCE_ARRAY_HEADER_INCLUDED
#define DISTANCE_ARRAY_HEADER_INCLUDED


namespace SimSite3D{
namespace geometry{

//! Make this a simple class to get it workiing -- then template it for generic
//! associated datatype(s) -- keep all functions here for now because if 
//! we use templates gcc requires them in header anyhow
class DistanceArray{
public:
  typedef std::vector< std::vector<const my_float_t *> > bin_vec;
  typedef bin_vec::const_iterator bin_vci;

  typedef std::multimap<bin_vci, const my_float_t*> bin2vert_mmap;
  typedef bin2vert_mmap::const_iterator bin2vert_mmci;

  //! Default constructor -- sometimes the array cannot be constructed 
  //! at class construction 
  DistanceArray()
  { init(); }

  DistanceArray(const my_float_t *positions, const uint num_pos,
                const my_float_t bin_width, const my_float_t max_dist)
  {
    init();
    setup_grid(positions, num_pos, bin_width, max_dist);
  }
  
  ~DistanceArray() {;}

  void
  setup_grid(const my_float_t *positions, const uint num_pos, 
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
    get_min_max_coords(positions, num_pos, min_pos, max_pos);
    determine_space(min_pos, max_pos, A_lower_corner, A_upper_corner, 
                    A_num_bins);
    populate_grid(positions, num_pos, A_max_dist, &A_bins);
//    A_bins_for_moved_pts->resize(A_num_bins[0]*A_num_bins[1]*A_num_bins[2]);
    A_bins_are_setup = true;
  }

  void 
  get_min_max_coords(const my_float_t *positions, const uint num_pos,
                     my_float_t *min_pos, my_float_t *max_pos)
  {
    std::copy(positions, positions + 3, min_pos);
    std::copy(positions, positions + 3, max_pos);
    const my_float_t *pos_end = positions + 3*num_pos;
    for(const my_float_t *pos = positions; pos < pos_end; pos += 3)
      for(int i = 0; i < 3; ++i){
        if(max_pos[i] < pos[i]) max_pos[i] = pos[i];
        if(pos[i] < min_pos[i]) min_pos[i] = pos[i];
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
  populate_grid(const my_float_t *positions, const uint num_pos,
                const my_float_t max_dist,
                std::vector< std::vector<const my_float_t*> > *bins)
  {
    // Would this be faster using dimension binning?
    //
    bins->resize(A_num_bins[0]*A_num_bins[1]*A_num_bins[2]);

    std::vector< std::vector<const my_float_t*> >::iterator cbin = bins->begin();
    //std::vector<my_float_t*>::iterator cbin = bins->begin();
    const my_float_t *pos_end = positions + 3*num_pos;
    my_float_t min_pt[3], max_pt[3];
    min_pt[0] = A_lower_corner[0]; 
    max_pt[0] = A_lower_corner[0] + A_bin_width; 
    for(int i = 0; i < A_num_bins[0]; ++i){
      min_pt[1] = A_lower_corner[1]; 
      max_pt[1] = A_lower_corner[1] + A_bin_width; 
      for(int j = 0; j < A_num_bins[1]; ++j){
        min_pt[2] = A_lower_corner[2]; 
        max_pt[2] = A_lower_corner[2] + A_bin_width; 
        for(int k = 0; k < A_num_bins[2]; ++k){
	  cbin->clear();
          for(const my_float_t *p = positions; p < pos_end; p += 3){
            if(min_pt[0] - max_dist <= p[0] &&
	       p[0] <= max_pt[0] + max_dist &&
               min_pt[1] - max_dist <= p[1] &&
	       p[1] <= max_pt[1] + max_dist &&
               min_pt[2] - max_dist <= p[2] &&
	       p[2] <= max_pt[2] + max_dist)
	    cbin->push_back(p);
	  }
          min_pt[2] += A_bin_width; 
          max_pt[2] += A_bin_width; 
	  ++cbin;
        }
        min_pt[1] += A_bin_width; 
        max_pt[1] += A_bin_width; 
      }
      min_pt[0] += A_bin_width; 
      max_pt[0] += A_bin_width; 
    }
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
} // end namespace SimSite3D
#endif
