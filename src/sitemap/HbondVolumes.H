#ifndef HBOND_VOLUMES_HEADER_FILE_INCLUDED
#define HBOND_VOLUMES_HEADER_FILE_INCLUDED

#include <hbond_points.H>



namespace ASCbase{

  class hbond_volume_t : public hbond_ideal_pt_base{
  public: 
    hbond_volume_t(const atom_vci hbond_atom, const atom_vci C_nbr_atom,
                   const atom_vci second_nbr_atom, 
                   const BoundingVolume &site_vol, const int cap_number,
                   const bool include_metals = false,
                   const alloc_t a = ALLOC_POSIITON)
      : hbond_ideal_pt_base(hbond_atom, C_nbr_atom, second_nbr_atom, site_vol,
                            cap_number, include_metals, a)
    {

    }



  };

  // Concentrate on the searching for now
  class HbondVolumes : public HbondBase{
  public:

  private:
  };


} // end namespace ASCbase

#endif
