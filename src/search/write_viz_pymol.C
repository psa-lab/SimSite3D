/*
 * $Source: /psa/share/repository/pfizer_proj/src/search/write_viz_pymol.C,v $
 * $Revision: 1.1 $
 * $Author: vanvoor4 $
 * $Date: 2006-04-11 20:36:00 $
 * 
 * $Log: not supported by cvs2svn $
 *
 */
 
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <err_handle.H>
#include <utils.H>
#include <write_viz_pymol.H>
#include <errno.h>

static const std::string accept_color = "  COLOR, 0.0, 1.0, 0.0";
static const std::string donor_color = "  COLOR, 0.0, 0.0, 1.0";
static const std::string donept_color = "  COLOR, 1.0, 1.0, 1.0";
static const std::string hydrophob_color = "  COLOR, 1.0, 0.0, 0.0";

static const double template_sp_sz = 0.2;
static const double search_sp_sz = 0.1;

int dump_points(FILE* fp, int& curr_idx, int max_idx, atom_pt atoms);
int dump_points(FILE* fp, int& curr_idx, int max_idx, 
		interaction_pt* atoms_ptr);

bool write_viz_pymol(global_data_pt global, char *filename, double score,
		     std::vector<int>& pairs)
{
  FILE *fp = fopen(filename, "w");
  if(fp == NULL){
    int my_err = errno;
    std::string msg = std::string("Unable to open ") + filename;
    err_error("write_viz_pymol", msg, my_err);
    return false;
  }

  fprintf(fp, "# Score                  : %5.3f\n", score);
  fprintf(fp, "# conformer              : %s\n", global->ligand_file_name); 

  /* Output geometric measures of docking */
  fprintf(fp, "# dist matrix error      : %5.3f\n", global->dist_matrix_err );
  fprintf(fp, "# rms error              : %5.3f\n", global->rms_err );
  fprintf(fp, "# matched triangle:\n" );
  int* ligand_matches = global->ligand_matches;
  int* template_matches = global->template_matches;
  atom_pt search_atoms = global->ligand->atoms;
  for(int i = 0; i < 3; i++ ) { 
    int index = ligand_matches[i];
    int template_act = global->template_interactions[template_matches[i]].act;
    fprintf(fp, "# Search point %3d  %-3s  %-5s to template point "
	    "%2d %s\n", search_atoms[index].atom_number, 
	    search_atoms[index].name, search_atoms[index].type_str, 
	    template_matches[i], 
	    action_type_to_string( (interactionType) template_act).c_str());
  }

  fprintf(fp, "from pymol.cgo import *    # get constants\n");
  fprintf(fp, "from pymol import cmd\n\n");
  fprintf(fp, "\nmy_spheres = [\n");
  interaction_pt* temp_pts = global->template_points;
  int act = temp_pts[0]->act;
  int max_i = global->number_of_template_points;
  for(int i = 0; i < max_i && act != -1; ){
    if(i != 0) fprintf(fp, ",\n");
    if(act == ACCEPTOR) fprintf(fp, "%s", accept_color.c_str());
    else if(act == DONOR) fprintf(fp, "%s", donor_color.c_str());
    else if(act == DONEPTOR)
      fprintf(fp, "%s", donept_color.c_str());
    else fprintf(fp, "%s", hydrophob_color.c_str());

    act = dump_points(fp, i, max_i, temp_pts);
  } 
  fprintf(fp, ",\n");
  atom_pt srch_pts = global->ligand->atoms;
  act = srch_pts[0].act;
  max_i = global->ligand->number_of_atoms - 
    global->ligand->number_of_added_hydrogens;
  for(int i = 0; i < max_i && act != -1; ){
    if(i != 0) fprintf(fp, ",\n");
    if(act == ACCEPTOR) fprintf(fp, "%s", accept_color.c_str());
    else if(act == DONOR) fprintf(fp, "%s", donor_color.c_str());
    else if(act == DONEPTOR)
      fprintf(fp, "%s", donept_color.c_str());
    else fprintf(fp, "%s", hydrophob_color.c_str());

    act = dump_points(fp, i, max_i, srch_pts);
  } 
  fprintf(fp, "\n  ]\ncmd.load_cgo(my_spheres, \'Template Points\')\n\n");


  fprintf(fp, "\nmy_lines = [\n  BEGIN, LINES,\n  COLOR, 1.0, 0.0, 0.0");
  for(int i = 0; i < pairs.size(); i+=2){
    float* temp_pos = temp_pts[pairs[i]]->pos;
    float* srch_pos = srch_pts[pairs[i+1]].pos;
    fprintf(fp, ",\n\n  VERTEX, %7.3f, %7.3f, %7.3f,\n", temp_pos[X],
	    temp_pos[Y], temp_pos[Z]);
    fprintf(fp, "  VERTEX, %7.3f, %7.3f, %7.3f", srch_pos[X], srch_pos[Y],
	    srch_pos[Z]);
  }
  fprintf(fp, ",\n  END\n  ]\ncmd.load_cgo(my_lines, \'Pairs\')\n");
  
  fclose(fp);
}



int dump_points(FILE* fp, int& curr_idx, int max_idx, atom_pt atoms)
{
  int curr_act = atoms[curr_idx].act;
  for(int idx = curr_idx ; idx < max_idx; idx++){
    if(atoms[idx].act != curr_act){
      curr_act = atoms[idx].act;
      curr_idx = idx;
      return atoms[idx].act;
    }
    fprintf(fp, ",\n  SPHERE, %7.3f, %7.3f, %7.3f, %7.2f", atoms[idx].pos[X],
	    atoms[idx].pos[Y], atoms[idx].pos[Z], search_sp_sz);
  }
  curr_idx = max_idx;
  return -1;
}

int dump_points(FILE* fp, int& curr_idx, int max_idx, 
		interaction_pt* atoms_ptr)
{
  int curr_act = atoms_ptr[curr_idx]->act;
  for(int idx = curr_idx ; idx < max_idx; idx++){
    if(atoms_ptr[idx]->act != curr_act){
      curr_act = atoms_ptr[idx]->act;
      curr_idx = idx;
      return atoms_ptr[idx]->act;
    }
    fprintf(fp, ",\n  SPHERE, %7.3f, %7.3f, %7.3f, %7.2f", 
	    atoms_ptr[idx]->pos[X], atoms_ptr[idx]->pos[Y], 
	    atoms_ptr[idx]->pos[Z], template_sp_sz);
  }
  curr_idx = max_idx;
  return -1;
}
