import ASCbasePy

class RandomAlignments:

  def __init__(self):
    pass

  def align(self, dbase_site, aligns, centroid, N=100, q_width=0.05, 
            half_side_len=0.25):
    """
    self -- calling class
    dbase_site -- (only need the point positions at this time -- need hphob point interface )
    aligns -- an ASCbasePy.search.vector_less__rigid_align_t__greater_() 
      instance
    centroid -- centroid of sitemap or ligand points -- use only ligand atoms
      in the query pocket
    N -- number of random alignments
    q_width -- width of a quaternion parameter bin (percent written as decimal)
    half_side_len -- half of the side length of a transformation coordinates
      cube

    Returns the parameters used to generate the orientations

    Assumptions:
      db_sites are aligned to query as best as possible using structure, etc
      The aligns vector is not explicitly cleared -- this allows for calling
      this method multiple times or adding the identity alignment before calling
      this function.

    """
#    aligns.clear()
    (orientations, orient_params) = \
      ASCbasePy.utils.gen_orientations(N, centroid, q_width=q_width,
                                       half_side_len=half_side_len, 
                                       overlap = 0.0)

    for i in range(N):
      (q,t) = orientations[i]
      x = ASCbasePy.search.rigid_align_t()
      R = q.get_ortho_rot_mat().reshape((9,))
      for j in range(9): x.R[j] = R[j]
      for j in range(3): x.T[j] = t[j]
      aligns.append(x)

    return orient_params
