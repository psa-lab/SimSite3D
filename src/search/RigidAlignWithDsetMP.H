/******************************************************************************
 * Copyright (c) 2011, Jeffrey Van Voorst
 * This file is part of the SimSite3D software project.
 *
 * Authors: Jeffrey Van Voorst, vanvoor4@msu.edu
 *          Leslie Kuhn, Ph.D., KuhnL@msu.edu 
 *
 *  SimSite3D is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  SimSite3D is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License along
 *  with this program; if not, write to the Free Software Foundation, Inc.,
 *  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *  You may also visit http://www.gnu.org/licenses/gpl-2.0.html.
 * 
 * $Date: $: Date of last commit
 * $Author: $: Author of last commit
 * $Rev: $: svn revision of last commit
 *
 * svn file: $Id: $
 * file location: $URL: $
 *****************************************************************************/
#ifndef RIGID_ALIGN_WITH_DATASET_MATCH_PRINTS_HEADER_INCLUDED
#define RIGID_ALIGN_WITH_DATASET_MATCH_PRINTS_HEADER_INCLUDED

#include <ScoreMapBase.H>

namespace SimSite3D{

class RigidAlignWithDsetMP : public rigid_align_t{
public:
  RigidAlignWithDsetMP() : rigid_align_t()
  { ; }

  RigidAlignWithDsetMP(const RigidAlignWithDsetMP& other) : rigid_align_t(other)
  {
    do_copy(other);
  }

  const RigidAlignWithDsetMP&
  operator=(const RigidAlignWithDsetMP& other)
  {
    if(&other != this){
      do_copy(other);
      rigid_align_t::operator=(other);
    }
    return *this;
  }

  virtual ~RigidAlignWithDsetMP()
  { ; }

  virtual void
  write_score_fields(std::ostream& out, const uint orient_num,
                     const bool wrote_ligs, const std::string& ext_SF_id_in,
                     const std::string& struct_id,
                     const std::string& lig_id) const
  {
    // First write the base class's score line
    rigid_align_t::write_score_fields(out, orient_num, wrote_ligs, ext_SF_id_in,
                                      struct_id, lig_id);
    std::vector<bool>::const_iterator mp;
    for(mp = dataset_mp.begin(); mp < dataset_mp.end(); ++mp)
      out << (*mp ? "1" : "0");
    out << "|";
  }

  virtual void
  write_score_fields(std::ostream& out) const
  {
    rigid_align_t::write_score_fields(out);
    std::vector<bool>::const_iterator mp;
    for(mp = dataset_mp.begin(); mp < dataset_mp.end(); ++mp)
      out << (*mp ? "1" : "0");
    out << "|";
  }

  virtual void
  get_score_field_labels(std::vector<std::string> *fields,
                         const bool normalize_score) const
  {
    rigid_align_t::get_score_field_labels(fields, normalize_score);
    fields->push_back("Dataset site point matchprint");
  }
 

  std::vector<bool> dataset_mp; 
private:

  void
  do_copy(const RigidAlignWithDsetMP& other)
  {
    dataset_mp = other.dataset_mp; 
  }

};
}
#endif
