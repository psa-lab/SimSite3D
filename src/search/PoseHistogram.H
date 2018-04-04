#include <math_basics.H>
#include <vector>
#include <list>
#include <cmath>
#include <iostream>

namespace pose_clustering{

  typedef std::vector<my_float_t> pose_t;
  typedef std::vector<pose_t> pose_vec;

  class PoseHistogram{
  public:
    PoseHistogram(uint bins_in, uint min_poses_in)
    {
      nbins = bins_in;
      min_poses = min_poses_in;
    }
  
    ~PoseHistogram() 
    {
      std::list<pose_vec*>::iterator li;
      for(li = bins.begin(); li != bins.end(); ++li) delete *li;
      bins.clear();
    }

    // Expects [Tx Ty Tz theta0 theta1 s]
    void histogram(const pose_vec& poses_in)
    {
      pose_vec poses(poses_in.begin(), poses_in.end());
      get_space_bounds(poses);

      // Again need to check if theta1 and s are circular
      circular.resize(6);
      circular[0] = circular[1] = circular[2] = circular[5] = false;
      circular[4] = false;
      circular[3] = true;

      subdivide_current_bin(&poses, 0);
    }

    std::list<pose_vec*>::const_iterator begin() const
    { return bins.begin(); }

    std::list<pose_vec*>::const_iterator end() const
    { return bins.end(); }
     
    size_t size() const { return bins.size(); }
  
  private:
    void get_space_bounds(const pose_vec& poses)
    {
      bnds.resize(12);
      std::fill(bnds.begin(), bnds.end(), my_float_max);
      for(uint i = 0; i < 6; ++i) bnds[i+6] = -bnds[i];
      pose_vec::const_iterator pose;
      for(pose = poses.begin(); pose != poses.end(); ++pose){
        pose_t tmp = *pose;
        for(uint i = 0; i < 3; ++i){
          if(tmp[i] < bnds[i]) bnds[i] = tmp[i];
          else if(tmp[i] > bnds[i+6]) bnds[i+6] = tmp[i];
        }
      }

      // theta0 is in [-pi,pi]
      bnds[3] = -M_PI - 0.0001;
      bnds[9] = M_PI + 0.0001;
      // theta2 is in [-pi/2, pi/2]
      bnds[4] = -M_PI_2 - 0.0001;
      bnds[10] = M_PI_2 + 0.0001;
      // s is in [0,1]
      bnds[5] = -0.0001;
      bnds[11] = 1.0001;
      
      // Kludge to enlarge the bounds enough so that the extreme points of the
      // bounded regions are not included in the poses.
      for(uint i = 0; i < 3; ++i){
        bnds[i] -= 0.0001;
        bnds[i+6] += 0.0001;
      }

      step_sizes.resize(6);
      for(uint i = 0; i < 6; ++i) step_sizes[i] = (bnds[i+6] - bnds[i]) / nbins;
    }

    // depends on implementation of push_back whether it is more efficient
    // to store intermediaries in vectors or lists
    void subdivide_current_bin(pose_vec *poses_in, const uint level)
    {
      // Init the bins for this level
      std::vector<pose_vec*> level_bins(nbins);
      for(uint i = 0; i < level_bins.size(); ++i) level_bins[i] = new pose_vec;

      //std::cout << "level is: " << level << std::endl;
      //std::cout << "lower bound is: " << bnds[level] << std::endl;
      //std::cout << "step size is: " << step_sizes[level] << std::endl;
      pose_vec::iterator pose;
      for(pose = poses_in->begin(); pose != poses_in->end(); ++pose){
        //for(uint i = 0; i < 6; ++i) std::cout << (*pose)[i] << " ";
        //std::cout << std::endl; 
        my_float_t val = ((*pose)[level] - bnds[level]) / step_sizes[level];
        //std::cout << "val is: " << val << std::endl;
        int bin = static_cast<int>(val);
	//std::cout << "bin: " << bin << "   nbins: " << nbins << std::endl;
        level_bins[bin]->push_back(*pose);
        
        // check for overlap
        my_float_t rem = val - bin;
        if(rem <= 0.2){
          if(bin) level_bins[bin - 1]->push_back(*pose);
          else if(circular[level]) level_bins[nbins - 1]->push_back(*pose);
        }else if(rem >= 0.8){
          if(bin + 1 < nbins) level_bins[bin + 1]->push_back(*pose);
          else if(circular[level]) level_bins[0]->push_back(*pose);
        }
      }
      poses_in->clear();
      if(level) delete poses_in;

      std::vector<pose_vec*>::iterator vi;
      for(vi = level_bins.begin(); vi < level_bins.end(); ++vi){
        //std::cout << "number of poses in bin: " << (*vi)->size() << std::endl;
        if((*vi)->size() < min_poses){
          (*vi)->clear();
          delete *vi;
          *vi = NULL;
        }else if(level < 5) subdivide_current_bin(*vi, level + 1);
        // Depends on how we wish to represent the bins -- for now we are trying
        // to get pose clustering working and to evaluate it -- we wish to 
        // have access to all the poses in each bin for stats/analysis
        else bins.push_back(*vi);
      }
    }

    std::vector<bool> circular;
    std::vector<my_float_t> step_sizes;
    std::vector<my_float_t> bnds;
    // Likely to be faster if use a vector of pose_vecs -- location should
    // trump cost of copying data.
    std::list<pose_vec*> bins;
    uint nbins;
    uint min_poses;
  };
}
