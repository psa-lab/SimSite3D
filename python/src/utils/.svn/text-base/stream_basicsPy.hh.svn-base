#include <fstream>
#include <iostream>

#ifndef STL_OFSTREAM_HEADER_FILE_INCLUDED
#define STL_OFSTREAM_HEADER_FILE_INCLUDED

class stl_ofstream{
public:
  stl_ofstream() { ; }
  
  ~stl_ofstream() { A_ofile.close(); }

  bool
  open(std::string fname) 
  { 
    A_ofile.open(fname.c_str(), std::ios_base::out|std::ios_base::trunc);
    if(A_ofile.fail()) return false; 
    return true;
  }

  void
  close()
  {
    A_ofile.close();
  }

  void 
  write(std::string str)
  {
    A_ofile << str;
  }

  std::ofstream&
  get() { return A_ofile; }

private:
  stl_ofstream(stl_ofstream&) {;}
  std::ofstream A_ofile;
};

#endif
