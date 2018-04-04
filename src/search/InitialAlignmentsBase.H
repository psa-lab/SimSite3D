/******************************************************************************
 * Copyright (c) 2010, Michigan State University (MSU) Board of Trustees.
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

#ifndef INITIAL_ALIGNMENTS_BASE_HEADER_FILE
#define INITIAL_ALIGNMENTS_BASE_HEADER_FILE

#include <ScoreMapBase.H>
#include <Sitemap.H>

namespace ASCbase{

//! The way to get around/handle template function specialization
/*! This is as per Herb Sutter, Peter Dimov, and Dave Abrahams 
 * who likely know the best way to handle C++ intricacies because they are
 * writing articles about C++ details and are
 * instrumental in building libraries such as Boost, etc.
 */
template < typename align_T >
struct align_impl;


//! Base class for initial alignment(s) of site maps / binding sites
class InitialAlignmentsBase{
public:

  //! Cstor 
  InitialAlignmentsBase()
  { ; }

  //! basic destruction
  virtual ~InitialAlignmentsBase()
  { ; }

  // Don't touch this
  template < typename align_T > bool
//  align(ModelSitemap *query_site, Sitemap& dset_site, 
  align(Sitemap& dset_site, std::vector<align_T> *alignments)
  {
    align_impl<align_T>::align(dset_site, alignments);
  }

/*
  //! Required alignment method
  template< typename align_T > bool
  align(Sitemap& dset_site, std::vector<align_T> *alignments)
  {
    std::cerr << "This function is supposed to be overloaded!!\n";
    return false;
  }
*/

//private:
};

//! specialize this template to get your alignment method
template < typename align_T >
struct align_impl
{
  // Overload this function
  static bool align(Sitemap& dset_site, std::vector<align_T> *alignments);
};





#if 0
bool
align_impl<rigid_align_t>::align(Sitemap& dset_site, std::vector<rigid_align_t> *alignments)
{
  std::cout << "test1\n";
}
#endif



} // end namespace ASCbase
#endif

