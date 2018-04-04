/******************************************************************************
 * Copyright (c) 2006,2007, Michigan State University (MSU) Board of Trustees.
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

/*
 * $Source: /psa/share/repository/pfizer_proj/src/search/PsuedoLigRmsd.C,v $
 * $Revision: 1.3 $
 * $Author: vanvoor4 $
 * $Date: 2008-01-04 15:42:42 $
 * 
 * $Log: not supported by cvs2svn $
 * Revision 1.2  2007/11/01 16:35:41  vanvoor4
 * Now pass in the parameters obj instead of those needed
 *
 * Revision 1.1  2007/10/11 16:21:15  vanvoor4
 * Initial checkin
 *
 * 
 */

#include <PsuedoLigRmsd.H>
#include <math_basics.H>
#include <iostream>

using namespace SimSite3D;

const std::string PsuedoLigRmsd::_fname = "PsuedoLigRmsd.C";

PsuedoLigRmsd::PsuedoLigRmsd(ModelSitemap* model_in, 
                             const SearchParameters& params)
 : ScoreRigidAlignments<std::less<my_float_t> >(model_in, params)
{
  model_lig = new mol2File(model->ligand_file_name());
}

PsuedoLigRmsd::~PsuedoLigRmsd()
{
  if(model_lig) delete model_lig;
  model_lig = 0;
}

my_float_t 
PsuedoLigRmsd::score(const DbaseSitemap& search, rigid_align_vi scores)
{
  mol2File search_lig(search.ligand_file_name());
  search_lig.transform(scores->R, scores->T);
  scores->score = simple_rmsd(*model_lig, search_lig); 
  return(scores->score);
}
