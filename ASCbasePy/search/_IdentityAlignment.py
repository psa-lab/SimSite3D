import ASCbasePy

class IdentityAlignment:

  def __init__(self):
    pass

  def align(self, dbase_site, aligns):
    """
    This function does not explicitly clear the aligns vector
    """
    x = ASCbasePy.search.rigid_align_t()
    x.R[0] = 1.0
    x.R[1] = 0.0
    x.R[2] = 0.0
    x.R[3] = 0.0
    x.R[4] = 1.0
    x.R[5] = 0.0
    x.R[6] = 0.0
    x.R[7] = 0.0
    x.R[8] = 1.0
    x.T[0] = 0.0
    x.T[1] = 0.0
    x.T[2] = 0.0
    aligns.append(x)
