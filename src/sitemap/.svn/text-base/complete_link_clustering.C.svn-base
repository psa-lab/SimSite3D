/******************************************************************************
 * Copyright (c) 2006,2007, Michigan State University (MSU) Board of Trustees.
 *   All rights reserved.
 *
 * This file is part of the ASCbase Software project.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *
 * Authors: Jeffrey Van Voorst, vanvoor4@msu.edu
 *          Leslie Kuhn, Ph.D., KuhnL@msu.edu 
 *****************************************************************************/

/*
 * $Source: /psa/share/repository/pfizer_proj/src/gen_points/complete_link_clustering.C,v $
 * $Revision: 1.6 $
 * $Author: vanvoor4 $
 * $Date: 2008-05-15 17:43:48 $
 *
 * $Log: not supported by cvs2svn $
 * Revision 1.5  2008/03/31 17:53:21  vanvoor4
 * Half hearted/baked attempt to get the complete link
 * clustering method to use an octree and other methods to speed
 * it up.
 *
 * Revision 1.4  2007/08/21 18:09:10  vanvoor4
 * Changed the point type from hbond_point to hphob_point
 *
 * Revision 1.3  2007/02/07 15:32:30  vanvoor4
 * Header files changed.  Added ASCbase namespace and they point_t conflict
 * required the necessary changes here as well.
 *
 * Revision 1.2  2006/10/20 13:19:58  vanvoor4
 * vectorified
 *
 * Revision 1.1  2006/03/30 18:39:11  vanvoor4
 * initial checkin
 *
 *
 */

#include <cstdio>
#include <cstdlib>
#include <defs.H>
#include <basics.H>
#include <complete_link_clustering.H>
#include <SimpleOctree.H>

using namespace ASCbase;

typedef enum my_stuff{
  NOT_SET,
  POINT_IN,
  POINT_OUT
}flag_type;

void complete_link_clustering(hphob_point_vec& hphob_pts, 
                              hphob_point_list* points, double threshold)
{
  double  **proximity_matrix;
  size_t     **cluster;
  size_t     *number_of_cluster_points;
  flag_type *marker;
  size_t     min_i = 0; 
  size_t     min_j = 0;

  hphob_point_vi pts_begin = hphob_pts.begin();
  hphob_point_vi pts_end = hphob_pts.end();

  size_t number_of_points = pts_end - pts_begin; 
  proximity_matrix = new double*[number_of_points];
  cluster = new size_t*[number_of_points];
  number_of_cluster_points = new size_t[number_of_points];
  std::fill(number_of_cluster_points, 
            number_of_cluster_points + number_of_points, 1);
  marker = new flag_type[number_of_points];

  // initialization 
  for(size_t i = 0; i < number_of_points; ++i){
    proximity_matrix[i] = new double[number_of_points];
    std::fill(proximity_matrix[i], proximity_matrix[i] + number_of_points, 
              std::numeric_limits<double>::infinity());
    proximity_matrix[i][i] = 0.0;
    cluster[i] = new size_t[number_of_points];
    cluster[i][0] = i;
    marker[i] = POINT_IN;
  }
  double min_dist = 999999.9;

  // Generate an upper diagonal pairwise distance matrix
  SimpleOctree<hphob_point_vci> Otree(threshold, 10); 
  Otree.build(pts_begin, pts_end); 
  for(hphob_point_vi pt_i = pts_begin; pt_i < pts_end; ++pt_i){
    const std::vector<hphob_point_vci> *tmp_pts =
      Otree.near_by_points(pt_i->pos, threshold);
    const size_t ii = pt_i - pts_begin;
    for(size_t idx = 0; idx < tmp_pts->size(); ++idx){
      const size_t jj = (*tmp_pts)[idx] - pts_begin;
      if(jj <= ii) continue;

      double my_dist = dist((*tmp_pts)[idx]->pos, pt_i->pos);
      if(my_dist <= threshold){
        proximity_matrix[ii][jj] = my_dist;
        if(my_dist < min_dist){
          min_dist = my_dist;
          min_i = ii;
          min_j = jj;
        }
      }
    }
  }

  int count = 0;
  while(min_dist < threshold){      
    ++count; 
    for(size_t i = 0; i < number_of_cluster_points[min_j]; ++i){
      cluster[min_i][number_of_cluster_points[min_i]] = cluster[min_j][i];
      number_of_cluster_points[min_i]++;
    }
    marker[min_j] = POINT_OUT;

    for(size_t i = 0; i < min_i; i++)
      if(marker[i] != POINT_OUT)
        if(proximity_matrix[i][min_i] < proximity_matrix[i][min_j] && 
           ((proximity_matrix[i][min_i] - proximity_matrix[i][min_j]) > MIN_DOUBLE || \
		 (proximity_matrix[i][min_i] - proximity_matrix[i][min_j]) < -MIN_DOUBLE) )
	    proximity_matrix[i][min_i] = proximity_matrix[i][min_j];

      for (size_t i = min_i + 1; i < min_j; i++ )
	if ( marker[i] != POINT_OUT )
	  if ( proximity_matrix[min_i][i] < proximity_matrix[i][min_j] && \
	       ( (proximity_matrix[min_i][i] - proximity_matrix[i][min_j]) > MIN_DOUBLE || \
		 (proximity_matrix[min_i][i] - proximity_matrix[i][min_j]) < -MIN_DOUBLE) )
	    proximity_matrix[min_i][i] = proximity_matrix[i][min_j];

      for (int j = min_j + 1; j < number_of_points; j++ )
	if ( marker[j] != POINT_OUT )
	    if ( proximity_matrix[min_i][j] < proximity_matrix[min_j][j] && \
		 ( (proximity_matrix[min_i][j] - proximity_matrix[min_j][j]) > MIN_DOUBLE || \
		   (proximity_matrix[min_i][j] - proximity_matrix[min_j][j]) < -MIN_DOUBLE) )
	    proximity_matrix[min_i][j] = proximity_matrix[min_j][j];

      min_dist = 999999.9;

      for(size_t i = 0; i < number_of_points; i++ )
        if(marker[i] == POINT_IN){
	  for(int j = i + 1; j < number_of_points; j++ )
	    if( marker[j] == POINT_IN)
	      if ( proximity_matrix[i][j] < min_dist && \
		   ( (proximity_matrix[i][j] - min_dist) > MIN_DOUBLE || \
		     (proximity_matrix[i][j] - min_dist) < -MIN_DOUBLE) )
		{
		  min_dist = proximity_matrix[i][j];
		  min_i = i;
		  min_j = j;
		}
	      }
    }
 
  for(size_t i = 0; i < number_of_points; i++)
    if(marker[i] != POINT_OUT){
      hphob_point_t tmp_pt;    
      memset(tmp_pt.pos, 0, 3 * my_float_size);
      for(size_t j = 0; j < number_of_cluster_points[i]; j++ )
        for(uint k = 0; k < 3; k++ ) 
          tmp_pt.pos[k] += (pts_begin + cluster[i][j])->pos[k];
      for(uint k = 0; k < 3; k++ ) 
        tmp_pt.pos[k] /= number_of_cluster_points[i];
      points->push_back(tmp_pt);
    }

  for(size_t i = 0; i < number_of_points; i++ ){
    delete[] ( cluster[i] );
    delete[] ( proximity_matrix[i] );
  }
  delete[] ( cluster );
  delete[] ( number_of_cluster_points );
  delete[] ( proximity_matrix );
  delete[] ( marker );
}
