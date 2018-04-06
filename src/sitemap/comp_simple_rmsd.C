#include <SitemapPointsFile.H>
#include <sstream>
#include <iomanip>

using namespace SimSite3D;

int main(int argc, char** argv)
{
  if(argc != 4){
    std::cerr << "Usage:  " << *argv << " <input_SIMSITE3D_test>.out <ligs_dir>" 
              << " <tag | NONE>"
              << std::endl;
    return -1;
  }

  std::ifstream in;
  if(!open_ifstream(in, argv[1])) return -1;

  my_float_t R[9]; 
  my_float_t T[3]; 
  std::string prev_lig;
  SitemapPointsFile* ref_lig = 0;
  SitemapPointsFile* move_lig = 0;
  my_float_t centroid[3];
  int lig_num = 0;
  for(std::string line; std::getline(in, line); ){
    if(line[0] == '#') continue;
  
    ++lig_num;
    std::vector<std::string> toks;
    string_tok(line, &toks, '|');

    // load ref lig
    std::ostringstream ofname;
    //std::vector<std::string> more_toks;
    //string_tok(toks[0], &more_toks, '/');
    std::string struct_id = toks[0].substr(0, toks[0].length() - 12);
    

    /*
    ofname << lig_fname.substr(0, lig_fname.length() - 10);
    lig_fname = argv[2] + std::string("/") + 
      lig_fname.substr(0, lig_fname.length() - 10) + "ade_l.pdb";
    if(lig_fname != prev_lig){
      prev_lig = lig_fname;
      if(ref_lig) delete ref_lig;
      ref_lig = new SitemapPointsFile(lig_fname); 
      if(move_lig) delete move_lig;
      move_lig = new SitemapPointsFile(lig_fname); 
      ref_lig->centroid_3D(centroid);
    }

    more_toks.clear();
    string_tok(toks[7], &more_toks, ' ');
    std::vector<std::string>::iterator t;
    int cnt = 0;
    for(t = more_toks.begin();t < more_toks.end(); ++t, ++cnt)
      my_strtof(*t,  R + cnt);

    cnt = 0;
    more_toks.clear();
    string_tok(toks[8], &more_toks, ' ');
    for(t = more_toks.begin();t < more_toks.end(); ++t, ++cnt)
      my_strtof(*t,  T + cnt);

    move_lig->revert();
    move_lig->transform(R,T);
    my_float_t moved_centroid[3];
    move_lig->centroid_3D(moved_centroid);
    std::cout << simple_rmsd(*ref_lig, *move_lig) << "|" 
              << dist(centroid, moved_centroid) << "|" << std::endl;;
    */

//    ofname << std::setfill('0') << std::setw(5) << lig_num << "_ADE_l.mol2";
//    move_lig->write(ofname.str());
  }
}
