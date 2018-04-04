#ifndef DISTANCE_ARRAY_HEADER_INCLUDED
#define DISTANCE_ARRAY_HEADER_INCLUDED


namespace ASCbase{
namespace geometry{

//! Make this a simple class to get it working -- then template it for generic
//! associated datatype(s) -- keep all functions here for now because if 
//! we use templates gcc requires them in header anyhow
/*!
 *  The goal of this method is to make things easier to update and hopefully at
 * the same time a bit faster to setup and access 
 */
class DistanceArray2{
public:
  typedef std::vector< std::vector<const my_float_t *> > bin_vec;
  typedef bin_vec::iterator bin_vi;
  typedef bin_vec::const_iterator bin_vci;

  typedef std::multimap<bin_vci, const my_float_t*> bin2vert_mmap;
  typedef bin2vert_mmap::const_iterator bin2vert_mmci;

  //! Default constructor -- sometimes the array cannot be constructed 
  //! at class construction 
  DistanceArray2()
  { init(); }

  DistanceArray2(const my_float_t *positions, const uint num_pos,
                 const my_float_t bin_width)
  {
    init();
    setup_grid(positions, num_pos, bin_width);
  }
  
  ~DistanceArray2() {;}

  void
  setup_grid(const my_float_t *positions, const uint num_pos, 
             const my_float_t bin_width)
  {
    if(A_bins_are_setup){
      std::cerr << "Bins are already setup - cannot set them up again\n";
      return;
    }
    A_bin_width = bin_width + 0.01;
    A_bin_width_inv = 1.0 / bin_width;
    //A_max_dist = max_dist;

    my_float_t min_pos[3], max_pos[3];
    get_min_max_coords(positions, num_pos, min_pos, max_pos);
    determine_space(min_pos, max_pos, A_lower_corner, A_upper_corner, 
                    A_num_bins);
    populate_grid(positions, num_pos, A_bin_width, &A_bins);
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
        A_bin_width * std::ceil(max_pos[i] / A_bin_width) + 2*A_bin_width;
      num_bins[i] = 
        static_cast<int>(std::ceil((max_corner[i] - min_corner[i]) 
                                   / A_bin_width));
    }
  }

  void
  populate_grid(const my_float_t *positions, const uint num_pos,
                const my_float_t bin_width,
                std::vector< std::vector<const my_float_t*> > *bins)
  {
    std::cout << A_num_bins[0] << " X " << A_num_bins[1] << " X " 
              << A_num_bins[2] << std::endl;

    bins->resize(A_num_bins[0]*A_num_bins[1]*A_num_bins[2]);
    A_current_bin_idz.resize(A_num_bins[0]*A_num_bins[1]*A_num_bins[2]);

    std::cout << "number of bins: " << A_bins.size() << std::endl;
    const my_float_t *pos_end = positions + 3*num_pos;
    for(const my_float_t *p = positions; p < pos_end; p += 3){
      int idx[3];
      for(int i = 0; i < 3; ++i)
        idx[i] = static_cast<int>(std::floor((p[i] - A_lower_corner[i]) * 
                                             A_bin_width_inv)); 

//      std::cout << "bin idz " << idx[0] << " " << idx[1] << " " 
//                << idx[2] << std::endl;

      // add p to the 27 bins (the computed index and +/- 1 for each coordinate)
      A_current_bin_idz.push_back(idx[0] * A_num_bins[1] * A_num_bins[2] +
                                  idx[1] * A_num_bins[2] + idx[2]);

      for(int i = -1; i < 2; ++i){
        bin_vi cbin = A_bins.begin();
        cbin += (idx[0] + i) * A_num_bins[1] * A_num_bins[2];
        cbin += (idx[1] - 1) * A_num_bins[2];
        for(int j = -1; j < 2; ++j){
          for(int k = -1; k < 2; ++k){
# if 0
            std::cout << "try this: " 
                      << ((idx[0] + i) * A_num_bins[1] * A_num_bins[2] +
                           (idx[1] + j) * A_num_bins[2] + (idx[2] + k))
                      << " "
                      << (A_bins.begin() + 
                          ((idx[0] + i) * A_num_bins[1] * A_num_bins[2] +
                           (idx[1] + j) * A_num_bins[2] + (idx[2] + k))) - A_bins.begin() << "\n";
            std::cout << "pushing back to bin # " << (cbin+k) - A_bins.begin()
                      << " " << (idx[0] + i) << " " << (idx[1] + j)
                      << " " << idx[2] + k << "\n";
#endif
            (cbin + idx[2] + k)->push_back(p);
          }
          cbin += A_num_bins[2];
        }
      }
    }

#if 0
    int cnt = 0;
    for(bin_vci bin = A_bins.begin(); bin < A_bins.end(); ++bin, ++cnt){
      if(bin->size() > 0) std::cout << "bin " << cnt << " has " << bin->size() 
                                    << " entries\n"; 
    }
#endif

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
#if 0
    std::cout << "pt is in bin idz " << idx[0] << " " << idx[1] << " " << idx[2] << "\n";
    std::cout << "pt is in bin # " 
              << idx[0] * A_num_bins[1] * A_num_bins[2] + idx[1] * A_num_bins[2] + idx[2] << "\n";
#endif
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
    //A_max_dist = my_float_max;
    A_bin_width = my_float_max;
    std::fill(A_upper_corner, A_upper_corner + 3, 0.0);
    std::fill(A_lower_corner, A_lower_corner + 3, 0.0);
    std::fill(A_num_bins, A_num_bins + 3, 0);
    A_bins_are_setup = false;
  }

//  my_float_t A_max_dist;  //!< Maximum distance for grid to be valid
  my_float_t A_bin_width; //!< Side length of a bin
  my_float_t A_bin_width_inv; //!< multiplicative inverse (1 over) side length
  my_float_t A_upper_corner[3];
  my_float_t A_lower_corner[3];
  int A_num_bins[3];
  bin_vec A_bins;
  std::vector<int> A_current_bin_idz;
  bin_vec A_bins_for_moved_pts;
  bool A_bins_are_setup;
};
} // end namespace geometry
} // end namespace ASCbase
#endif
