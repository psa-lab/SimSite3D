#include <SitemapPointsFile.H>
#include <cstring>
#include <sstream>
#include <iomanip>
#include <map>

using namespace SimSite3D;

int main(int argc, char** argv)
{
  if(argc != 4){
    std::cerr << "Usage:  " << *argv << " <hits>.out <starting_dset_dir> "
              << "<ref_dset_dir> " << std::endl;
    std::cerr << "Where <hits>.out is the SimSite3D results file\n" 
              << "      <starting_dset_dir> is the dir holding the initial site map orientations (before alignments)\n"
              << "      <ref_dset_dir> is the dir holding the reference site map orientations (i.e. 0.0 RMSD orientations\n" 
              << std::endl;
    std::cerr << "The site maps need to be identical in terms of filenames and all data in the files except for the positions of the points\n";
    return -1;
  }

  std::map<std::string,bool> ade_ids;
  ade_ids["1G55_SAH"] = true;
  ade_ids["1MSK_SAM"] = true;
  ade_ids["1MRJ_ADN"] = true;
  ade_ids["1F3L_SAH"] = true;
  ade_ids["1MXT_FAE"] = true;
  ade_ids["1KA1_A3P"] = true;
  ade_ids["1B37_FAD"] = true;
  ade_ids["1BX4_ADN"] = true;
  ade_ids["2DPM_SAM"] = true;
  ade_ids["1KJQ_ADP"] = true;
  ade_ids["1C1D_NAD"] = true;
  ade_ids["1DL5_SAH"] = true;
  ade_ids["1F20_NAP"] = true;
  ade_ids["1KOL_NAD"] = true;
  ade_ids["1T2D_NAD"] = true;
  ade_ids["1J09_ATP"] = true;
  ade_ids["1KPF_AMP"] = true;
  ade_ids["1RLZ_NAD"] = true;
  ade_ids["1XVA_SAM"] = true;
  ade_ids["1HP1_ATP"] = true;
  ade_ids["1O2D_NAP"] = true;
  ade_ids["1KRH_FAD"] = true;
  ade_ids["1AF7_SAH"] = true;

  std::ifstream in;
  if(!open_ifstream(in, argv[1])) return -1;

  const std::string init_site_dir = argv[2];
  const std::string ref_site_dir = argv[3];

  my_float_t R[9]; 
  my_float_t T[3];  
  std::string prev_struct_id;
  SitemapPointsFile* ref_site = 0;
  SitemapPointsFile* move_site = 0;
  Quaternion Q_id;
  my_float_t ref_centroid[3];
  int num = 0;
  for(std::string line; std::getline(in, line); ){
    if(line[0] == '#') continue;

    // parse the line into columns
    ++num;
    std::vector<std::string> toks;
    string_tok(line, &toks, '|');

    // Get the struct_id so we know the name of the sitemap file
    std::string struct_id;
    if(toks[0].find("init_") != std::string::npos){
      std::string tmp = toks[0].substr(5);
      toks[0] = tmp;
    }else if(toks[0].find("preIK_") != std::string::npos){ 
      std::string tmp = toks[0].substr(6);
      toks[0] = tmp;
    }
  
    // becuase I am too lazy to handle norm db another way
    std::map<std::string, bool>::const_iterator found = 
      ade_ids.find(toks[0]);
    if(found == ade_ids.end()){
      std::cout << "-1|-1|-1|\n";
      continue;
    }

    // Remove the ".mol2" extension if it exists, and append an underscore
    if(toks[0].find("mol2") == std::string::npos)
      struct_id = toks[0] + "_";
    else{
      std::vector<std::string> more_toks;
      string_tok(toks[0], &more_toks, '_');
      for(size_t i = 0; i < more_toks.size() - 2; ++i)
        struct_id += more_toks[i] + "_";
    }

    if(prev_struct_id != struct_id){
      prev_struct_id = struct_id;
      std::string init_fname = init_site_dir;
      init_fname += std::string("/") + struct_id + "s.pdb";
      std::string ref_fname = ref_site_dir;
      ref_fname += std::string("/") + struct_id + "s.pdb";

      if(ref_site) delete ref_site;
      ref_site = new SitemapPointsFile(ref_fname);
      if(move_site) delete move_site;
      move_site = new SitemapPointsFile(init_fname);
      ref_site->centroid_3D(ref_centroid);
    }

    // Get the transformation
    std::vector<std::string> more_toks;
    string_tok(toks[2], &more_toks, ' ');
    std::vector<std::string>::iterator t;
    int cnt = 0;
    for(t = more_toks.begin();t < more_toks.end(); ++t, ++cnt)
      my_strtof(*t,  R + cnt);

    cnt = 0;
    more_toks.clear();
    string_tok(toks[3], &more_toks, ' ');
    for(t = more_toks.begin();t < more_toks.end(); ++t, ++cnt)
      my_strtof(*t,  T + cnt);

    move_site->revert();
    my_float_t R_transpose[9];
    for(size_t i = 0; i < 3; ++i)
      for(size_t j = 0; j < 3; ++j) R_transpose[3*j + i] = R[3*i + j];
    move_site->transform(R_transpose,T);
    my_float_t moved_centroid[3];
    move_site->centroid_3D(moved_centroid);
    Quaternion Q(R,9);
    std::cout << simple_rmsd(*ref_site, *move_site) << "|"
              << dist(ref_centroid, moved_centroid) << "|"
              << std::acos(1.0 - distance(Q, Q_id, 1.0)) << "|" << "\n";
  }
}
