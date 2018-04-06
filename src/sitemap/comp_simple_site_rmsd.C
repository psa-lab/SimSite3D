#include <SitemapPointsFile.H>
#include <cstring>
#include <sstream>
#include <iomanip>

using namespace SimSite3D;

int main(int argc, char** argv)
{
  if(argc != 4){
    std::cerr << "Usage:  " << *argv << " <input_SIMSITE3D_test>.out <dbase_site_dir>" 
              << " <tag | NONE>"
              << std::endl;
    return -1;
  }

  std::ifstream in;
  if(!open_ifstream(in, argv[1])) return -1;

  my_float_t R[9]; 
  my_float_t T[3]; 
  std::string prev_struct_id;
  SitemapPointsFile* ref_site = 0;
  SitemapPointsFile* move_site = 0;
  Quaternion Q_id;
  my_float_t centroid[3];
  int num = 0;
  for(std::string line; std::getline(in, line); ){
    if(line[0] == '#') continue;
  
    ++num;
    std::vector<std::string> toks;
    string_tok(line, &toks, '|');

    // The output of SimSite3D has changed between v3.0 and v3.3 --
    // if the pocket size is 0 ligands are not written and the struct_id
    // is (by default) printed in column 1
    std::string struct_id;
    if(toks[0].find("init_") != std::string::npos){
      std::string tmp = toks[0].substr(5);
      toks[0] = tmp;
    }else if(toks[0].find("preIK_") != std::string::npos){
      std::string tmp = toks[0].substr(6);
      toks[0] = tmp;
    }


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
      std::string site_fname = argv[2];
      if(strcmp(argv[3], "NONE"))
        site_fname += std::string("/") + struct_id + argv[3] + "_s.pdb";
      else  site_fname += std::string("/") + struct_id + "s.pdb";
      if(ref_site) delete ref_site;
      
      ref_site = new SitemapPointsFile(site_fname);
      if(move_site) delete move_site;
      move_site = new SitemapPointsFile(site_fname);
      ref_site->centroid_3D(centroid);
    } 

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
              << dist(centroid, moved_centroid) << "|" 
              << std::acos(1.0 - distance(Q, Q_id, 1.0)) << "|" << "\n";

//    ofname << std::setfill('0') << std::setw(5) << << "_ADE_l.mol2";
//    move_lig->write(ofname.str());
  }
}
