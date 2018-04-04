#include <string>
#include <iostream>
#include <mat_ops.H>
#include <SimpleTrimesh.H>
#include <TransformableTrimesh.H>
#include <ImmovableTrimesh.H>
#include <ImmovableTrimeshThree.H>
#include <DistanceArray.H>
#include <Timer.H>

using namespace SimSite3D;
using SimSite3D::geometry::DistanceArray;
using SimSite3D::geometry::SimpleTrimeshTwo;
using SimSite3D::geometry::TransformableTrimesh;
using SimSite3D::geometry::ImmovableTrimesh;
using SimSite3D::geometry::ImmovableTrimeshThree;

void 
get_extents(const SimpleTrimeshTwo &mesh, my_float_t *min_c, my_float_t *max_c);

void
check_all_triangles(const SimpleTrimesh &mesh, const my_float_t *pos,
                    my_float_t *cp, my_float_t *dist);

bool
test_vertices(const SimpleTrimesh &mesh);

bool
test_vertices_using_mess(const ImmovableTrimesh &mesh);

bool
check_correspondences(const my_float_t *cp0, const my_float_t *cp1, 
                     const my_float_t *ref_dists, const int num_pts);

int 
main(int argc, char** argv)
{
  const int num_runs = 100; 
  std::cout << "Number of non-brute force runs: " << num_runs << std::endl;
  // Farthest away from a surface
  const my_float_t radius = 1.5;

  if(argc != 3){
    std::cerr << "Usage: " << argv[0] << " <query_file_name>.vert <db_file>.vert\n";
    return -1;
  }

  std::string q_fname = argv[1];
  if(q_fname.rfind(".vert") != std::string::npos){
    std::string tmp = q_fname.substr(0, q_fname.length() - 5);
    q_fname = tmp;
  }

  std::string db_fname = argv[2];
  if(db_fname.rfind(".vert") != std::string::npos){
    std::string tmp = db_fname.substr(0, db_fname.length() - 5);
    db_fname = tmp;
  }

  Timer my_timer;
  my_timer.start();
  SimpleTrimesh old_q_mesh(q_fname);
  old_q_mesh.build_vert_octree();
  my_float_t timer_val[3];
  my_timer.get(timer_val, timer_val + 1, timer_val + 2);
  std::cout << "Finished" << std::endl;
  std::cout.precision(2);
  std::cout << "Initialization of SimpleTrimesh + build_vert_octree time: " 
            << timer_val[2] << "\n";

  my_float_t start_time[3];
  my_timer.get(start_time, start_time + 1, start_time + 2);
  SimpleTrimesh old_db_mesh(db_fname);
  old_db_mesh.build_vert_face_map();
  my_float_t end_time[3];
  my_timer.get(end_time, end_time + 1, end_time + 2);
  std::cout << "Finished" << std::endl;
  std::cout.precision(2);
  std::cout << "Initialization of SimpleTrimesh + build_vert_face_map time: " 
            << end_time[2] - start_time[2] << "\n";

  TransformableTrimesh new_q_mesh(q_fname);
  std::cout << "nubmer of vertices in new_q_mesh:" << new_q_mesh.number_of_vertices() << std::endl;
  TransformableTrimesh new_xformable_db_mesh(db_fname);

  my_timer.get(start_time, start_time + 1, start_time + 2);
  ImmovableTrimesh new_db_mesh(db_fname);
  std::cout << "Finished" << std::endl;
  std::cout.precision(2);
  std::cout << "Initialization of ImmovableTrimesh time: " 
            << end_time[2] - start_time[2] << "\n";

  my_timer.get(start_time, start_time + 1, start_time + 2);
  ImmovableTrimeshThree newest_db_mesh(db_fname);
  std::cout << "Finished" << std::endl;
  std::cout.precision(2);
  std::cout << "Initialization of ImmovableTrimeshThree time: " 
            << end_time[2] - start_time[2] << "\n";




  // Test for 2 different meshes

  // Method 0 all triangles
  std::cout << "Testing closest point using all triangles ... " << std::flush;
  int N = old_q_mesh.number_of_vertices();
  my_float_t *closest_points_0 = new my_float_t[4*N];
  my_float_t *dists_0 = closest_points_0 + 3*N;
  std::fill(dists_0, dists_0 + N, my_float_max);
  my_float_t *cp = closest_points_0;
  my_float_t *d = dists_0;
  my_timer.get(start_time, start_time + 1, start_time + 2);
//  for(int zz = 0; zz < num_runs; ++zz){
    cp = closest_points_0;
    d = dists_0;
    const my_float_t *v_end = 
      old_q_mesh.vertices_begin() + 3*old_q_mesh.number_of_vertices();
    for(const my_float_t *v = old_q_mesh.vertices_begin(); v < v_end; v += 3)
    {
      check_all_triangles(old_db_mesh, v, cp, d);
      ++d;
      cp += 3;
    }
//  }
  my_timer.get(end_time, end_time + 1, end_time + 2);
  std::cout << "Finished" << std::endl;
  std::cout.precision(2);
  std::cout << "Checking all triangles total time: " 
            << end_time[2] - start_time[2]<< "\n";

  // Method 1 SimpleTrimesh
  my_float_t prev_timer_val[3];
  std::cout << "Testing closest point using SimpleTrimesh ... " << std::flush;
  my_timer.get(prev_timer_val, prev_timer_val + 1, prev_timer_val + 2);
  for(int zz = 0; zz < num_runs; ++zz){
    size_t num_points;
    my_float_t RMSE;
    SimpleTrimesh::correspond_map cmap;
    old_q_mesh.compare(old_db_mesh, 1.5, &num_points, &RMSE, &cmap);
  }
  my_timer.get(timer_val, timer_val + 1, timer_val + 2);
  std::cout << "Finished" << std::endl;
  std::cout << "SimpleTrimesh total time: " 
            << timer_val[2] - prev_timer_val[2] << "\n";

  // Method 2 ImmovableTrimesh -- closest_point
  std::cout << "Testing closest point using ImmovableTrimesh ... " << std::flush;
  my_timer.get(prev_timer_val, prev_timer_val + 1, prev_timer_val + 2);
  my_float_t *closest_points_2 = new my_float_t[4*N];
  my_float_t *dists_2 = closest_points_2 + 3*N;
  std::fill(dists_2, dists_2 + N, my_float_max);
  for(int zz = 0; zz < num_runs; ++zz){
    d = dists_2;
    cp = closest_points_2;

    v_end = new_q_mesh.vertices_begin() + 3*new_q_mesh.number_of_vertices();
    for(const my_float_t *v = new_q_mesh.vertices_begin(); v < v_end; v += 3)
    {
      new_db_mesh.closest_point(v, &(*d), cp, 1.5);
      ++d;
      cp += 3;
    }
  }
  my_timer.get(timer_val, timer_val + 1, timer_val + 2);
  std::cout << "Finished" << std::endl;
  std::cout << "ImmovableTrimesh total time: " 
            << timer_val[2] - prev_timer_val[2] << "\n";

  // Method 3 TransformableTrimesh -- compare 
  std::cout << "Testing TransformableTrimesh::compare ... " << std::flush;
  my_timer.get(prev_timer_val, prev_timer_val + 1, prev_timer_val + 2);
  my_float_t *closest_points_3 = new my_float_t[4*N];
  my_float_t *dists_3 = closest_points_3 + 3*N;
  new_xformable_db_mesh.compare(new_q_mesh.vertices_begin(), 
                                new_q_mesh.number_of_vertices(), 
                                closest_points_3, dists_3, 0, 0, 0, 1.5);
  my_timer.get(timer_val, timer_val + 1, timer_val + 2);
  std::cout << "Finished" << std::endl;
  std::cout << "TransformableTrimesh::compare total time: " 
            << timer_val[2] - prev_timer_val[2] << "\n";

  // Method 4 ImmovableTrimesh -- compare
  std::cout << "Testing ImmovableTrimesh::compare ... " << std::flush;
  my_timer.get(prev_timer_val, prev_timer_val + 1, prev_timer_val + 2);
  my_float_t *closest_points_4 = new my_float_t[4*N];
  my_float_t *dists_4 = closest_points_4 + 3*N;
  for(int zz = 0; zz < num_runs; ++zz)
    new_db_mesh.compare(new_q_mesh.vertices_begin(), 
                        new_q_mesh.number_of_vertices(), 
                        closest_points_4, dists_4, 0, 0, 0, 1.5);
  my_timer.get(timer_val, timer_val + 1, timer_val + 2);
  std::cout << "Finished" << std::endl;
  std::cout << "ImmovableTrimesh::compare total time: " 
            << timer_val[2] - prev_timer_val[2] << "\n";

  // Method 5 ImmovableTrimeshThree -- compare
  std::cout << "Testing ImmovableTrimeshThree::compare ... " << std::flush;
  my_timer.get(prev_timer_val, prev_timer_val + 1, prev_timer_val + 2);
  my_float_t *closest_points_5 = new my_float_t[4*N];
  my_float_t *dists_5 = closest_points_0 + 3*N;
  for(int zz = 0; zz < num_runs; ++zz)
    newest_db_mesh.compare(new_q_mesh.vertices_begin(), 
                           new_q_mesh.number_of_vertices(), 
                           closest_points_5, dists_5, 0, 1.5);
  my_timer.get(timer_val, timer_val + 1, timer_val + 2);
  std::cout << "Finished" << std::endl;
  std::cout << "ImmovableTrimesh::compare total time: " 
            << timer_val[2] - prev_timer_val[2] << "\n";

  std::cout.precision(3);
  std::cout << "V0: " << new_q_mesh.vertices_begin()[0] << " "
            << new_q_mesh.vertices_begin()[1] << " "
            << new_q_mesh.vertices_begin()[2] << "\n";
  std::cout << "Good Correspond: " << closest_points_0[0] << " "
            << closest_points_0[1] << " "
            << closest_points_0[2] << "\n";
  std::cout << "Bad (new) Correspond: " << closest_points_5[0] << " "
            << closest_points_5[1] << " "
            << closest_points_5[2] << "\n";


  // Check corrspondences
//  std::cout << "Correspondences between 0 and 2\n";
//  check_correspondences(closest_points_0, closest_points_2, 
//                       dists_0, old_q_mesh.number_of_vertices());
//  std::cout << "\nCorrespondences between 0 and 3\n";
//  check_correspondences(closest_points_0, closest_points_3, 
//                       dists_0, old_q_mesh.number_of_vertices());
  std::cout << "\nCorrespondences between 0 and 5\n";
  check_correspondences(closest_points_0, closest_points_5, 
                        dists_0, old_q_mesh.number_of_vertices());


  delete [] closest_points_0; 
  delete [] closest_points_2;
//  delete [] closest_points_3;
  delete [] closest_points_4;
  delete [] closest_points_5;
}

void 
get_extents(const SimpleTrimeshTwo &mesh, my_float_t *min_c, my_float_t *max_c)
{
  size_t N = mesh.number_of_vertices();
  const my_float_t *v = mesh.vertices_begin();
  std::copy(v, v+3, min_c);
  std::copy(v, v+3, max_c);
  for(size_t i = 0; i < N; ++i, v += 3)
    for(int j = 0; j < 3; ++j){
      if(max_c[j] < v[j]) max_c[j] = v[j]; 
      if(min_c[j] < v[j]) min_c[j] = v[j]; 
    }
}

// Brute force method: check each and every triangle in the mesh for the 
// point on the mesh that is closest to the given position (pos).
void
check_all_triangles(const SimpleTrimesh &mesh, const my_float_t *pos,
                    my_float_t *cp, my_float_t *dist)
{
  *dist = my_float_max;

  my_float_t abs_saved_d = my_float_max;
  for(SimpleTrimesh::face_vci f = mesh.faces_begin(); f < mesh.faces_end(); ++f)
  {
    my_float_t d, tmp_cp[3];
    // Check if corresponding point is on a face
    bool rv = corresponding_point(pos, f->vertices[0], f->vertices[1],
                                  f->vertices[2], d, &d, tmp_cp);
    if(rv){
      my_float_t abs_d = d;
      if(abs_d < 0) abs_d *= -1.0;
      if(abs_d < abs_saved_d){
        abs_saved_d = abs_d;
        *dist = d;
        std::copy(tmp_cp, tmp_cp + 3, cp);
      }
    }else{
      for(uint j = 0; j < 3; ++j){
        rv = corresponding_point(pos, f->vertices[j], f->vertices[(j+1) % 3], 
                                 &d, tmp_cp);
        my_float_t abs_d = d;
        if(abs_d < 0) abs_d *= -1.0;
        if(abs_d < abs_saved_d){
          abs_saved_d = abs_d;
          *dist = d;
          std::copy(tmp_cp, tmp_cp + 3, cp);
        }
      }
    }
  }
}

bool
test_vertices(const SimpleTrimesh &mesh)
{
  const my_float_t dist_tol = 0.001;

  const my_float_t *v_end = 
    mesh.vertices_begin() + 3 * mesh.number_of_vertices();
  for(const my_float_t *v = mesh.vertices_begin(); v < v_end; v += 3)
  {
    my_float_t cp[3], dist;
    check_all_triangles(mesh, v, cp, &dist);
    if(dist < 0.0) dist *= -1.0;
    if(dist > dist_tol){
      std::cerr << "check_all_triangles failed to find correct vertex in "
                << "test_vertices\n";
      return false;
    }
  }

  return true;  
}

bool
test_vertices_using_mess(const ImmovableTrimesh &mesh)
{
  const my_float_t dist_tol = 0.001;

  const my_float_t *v_end = 
    mesh.vertices_begin() + 3 * mesh.number_of_vertices();
  for(const my_float_t *v = mesh.vertices_begin(); v < v_end; v += 3)
  {
    my_float_t cp[3], dist;
    mesh.closest_point(v, &dist, cp);
    if(dist < 0.0) dist *= -1.0;
    if(dist > dist_tol){
      std::cerr << "check_all_triangles failed to find correct vertex for "
                << "vertex " << (v - mesh.vertices_begin())/3 << " in "
                << "test_vertices_using_mess\n";
      return false;
    }
  }

  return true;  
}

bool
check_correspondences(const my_float_t *cp0, const my_float_t *cp1, 
                     const my_float_t *ref_dists, const int num_pts)
{
  std::cout.precision(6);
  const my_float_t *p0 = cp0;
  const my_float_t *p1 = cp1;
  int cnt1 = 0, cnt2 = 0;
  for(int i = 0; i < num_pts; ++i, p0 += 3, p1 += 3){
#if 0
    std::cout << "P0[0]: " << p0[0] << " " << p0[1] << " " << p0[2] << "\n";
    std::cout << "P1[0]: " << p1[0] << " " << p1[1] << " " << p1[2] << "\n";
    std::cout << "distance is: " 
              << std::sqrt((p0[0]-p1[0])*(p0[0]-p1[0])+
                           (p0[1]-p1[1])*(p0[1]-p1[1])+
                           (p0[2]-p1[2])*(p0[2]-p1[2])) << "\n";
#endif

    my_float_t diff = dist(p0, p1);
    std::cout << i << " "; 
    if(-1.5 <= ref_dists[i] && ref_dists[i] <= 1.5){
       ++cnt1;
    }else{

    }
    std::cout << ref_dists[i];
    std::cout << "\n";

    
    std::cout << i << " Distance was: " << ref_dists[i] 
              << "     diff: " << diff << std::endl;
#if 0
    if(diff > 1e-07){
      std::cout << i << " Distance was: " << ref_dists[i] 
                << "     diff: " << diff << std::endl;
    }else{
      std::cout << i << " " << dist(p0, p1) << std::endl;
      if(dist(p0, p1) <= 1.5) ++cnt2;
    }
#endif
  }

  std::cout << "Count 1 is " << cnt1 << "\n";
//  std::cout << "Count 2 is " << cnt2 << "\n";
}
