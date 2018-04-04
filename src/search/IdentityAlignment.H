/******************************************************************************
 * Copyright (c) 2010, Michigan State University (MSU) Board of Trustees.
 *   All rights reserved.
 *
 * This file is part of the SimSite3D Software project.
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

#ifndef IDENTITY_ALIGNMENT_HEADER_FILE
#define IDENTITY_ALIGNMENT_HEADER_FILE

#include <ScoreMapBase.H>
//#include <Sitemap.H>
#include <ModelSitemap.H>

namespace SimSite3D{

//! Used when user wants to score the sites as given -- generally not
//! useful in the case of a standard screening dataset
class IdentityAlignment{
public:

  //! Cstor 
  IdentityAlignment()
  { ; }

  //! basic destruction
  virtual ~IdentityAlignment()
  { ; }

  //! Required alignment method
  template< typename align_T > bool
  align(ModelSitemap *query_site, Sitemap& dset_site, 
        std::vector<align_T> *alignments)
  {
    align_T my_align;
    std::fill(my_align.R, my_align.R + 9, 0.0);
    my_align.R[0] = my_align.R[4] = my_align.R[8] = 1.0;
    std::fill(my_align.T, my_align.T + 3, 0.0);
    my_float_t Q[] = {1.0, 0.0, 0.0, 0.0};
    my_align.Q = Quaternion(Q, 4);
    alignments->push_back(my_align);
    return true;
  }

};
} // end namespace SimSite3D
#endif
