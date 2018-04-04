#include <ScoreMapBase.H>

using namespace ASCbase;

void
rigid_align_t::write_score_fields(std::ostream& out, const uint orient_num,
                                  const bool wrote_ligs,
                                  const std::string& ext_SF_id_in,
                                  const std::string& struct_id, 
                                  const std::string& lig_id) const
{
  // struct id or ligand frag fname
  if(wrote_ligs)
    out << lig_id << "_" << std::setfill('0') << std::right
        << std::setw(5) << orient_num << "_f.mol2" << "|";
  else out << struct_id << "|";
  //else if(!include_struct_id_field) out << struct_id << "|";
  //else out << "|";

  // Write the different score fields
  write_score_fields(out);
  
  if(ext_SF_id_in.length())
    for(uint i = 0; i < ext_scores.size(); ++i) out << ext_scores[i] << "|";
}

void
rigid_align_t::write_score_fields(std::ostream &out) const
{
  // Score is likely accurate to at best 2 or 3 decimal places 
  out << std::fixed << std::setprecision(3) << std::left << std::setfill(' ')
      << score << "|";

  // Transformation (rotation followed by translation) 
  out << std::fixed << std::setprecision(14);
  // Our rotation matrices are to be muliplied from the left by the 
  // positions, but the output matrices are for post multiplication by 
  // vectors
  std::string follower = "        |";
  for(size_t i = 0; i < 3; ++i)
    for(size_t j = 0; j < 3; ++j){
      out << R[3*j + i] << follower[3*i + j];
    }
  for(uint i = 0; i < 3; ++i) out << T[i] << follower[6 + i];


  // Sitemap points match string
  std::vector<bool>::const_iterator p;
  for(p = match_print.begin(); p < match_print.end(); ++p){
    if(*p) out << "1";
    else out << "0";
  }
  out << "|";

  // Ligand fragment atoms string
  std::vector<bool>::const_iterator f;
  for(f = frag_atoms_flags.begin(); f < frag_atoms_flags.end(); ++f){
    if(*f) out << "1";
    else out << "0";
  }
  out << "|";

  // Orientation/Alignment features
  if(terms.size()){
    out << std::fixed << std::setprecision(7) << terms[0];
    for(size_t i = 1; i < terms.size(); ++i) out << " " << terms[i];
  }
  out << "|" << std::fixed << std::setprecision(7);
}

void
rigid_align_t::get_score_field_labels(std::vector<std::string> *fields, 
                                      const bool normalize_score) const
{
  std::string tmp;
  if(normalize_score)
    tmp = "Normalized ASCbase alignment score of target to query";
  else tmp = "Raw ASCbase alignment score of target to query";
  fields->push_back(tmp);

  fields->push_back("Rotation matrix to align target to query");
  fields->push_back("Translation vector to move target to query");
  tmp = "Match print of the query's sitemap points satisfied by ";
  tmp += "sitemap points\n#     in the database hit";
  fields->push_back(tmp);
  tmp = "Ligand fragment binary string:  1 or 0 in nth position ";
  tmp += "implies that the\n#     nth mol2 ligand atom is or is not in ";
  tmp += "the mol2 ligand fragment (resp.)";
  fields->push_back(tmp);
  fields->push_back("scoring function terms/alignment features");
}
