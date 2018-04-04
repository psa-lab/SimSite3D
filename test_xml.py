#<hydrogen_bond_clouds>
#  <atom>
#    <chainID></chainID>
#    <resSeq></resSeq>
#    <iCode></iCode>
#    <name></name>
#    <cloud>
#      <number></number>
#      <position></position>
#      <neighbor>
#        <position></position>
#      </neighbor>
#    </cloud>
#  </atom>
#</hydrogen_bond_clouds>

import xml.dom.minidom

class test_xml:

  def __init__(self, csv_fname):
    self.dom = xml.dom.minidom.parse(csv_fname)
    tmp = self.dom.getElementsByTagName("hydrophobic_points")
    self.hphob_pts_el = tmp[0]

    # create the clouds section
    clouds = self.dom.createElement("hydrogen_bond_clouds")
    clouds = self.dom.documentElement.insertBefore(clouds, self.hphob_pts_el)
    self.clouds = clouds

#    # create a template cloud node
#    cloud = self.dom.createElement("cloud")
#    items = ["chainID", "resSeq", "iCode", "name"]
#    for item in items: self.add_text_node(cloud, item)
#
#    # create a template point node
#    point = self.dom.createElement("point")
#    self.add_text_node(point, "number")
#    self.add_text_node(point, "position")
#    cloud.appendChild(point)
#
#    # Neighbors will need to added later
#    self.cloud = cloud

  def add_text_node(self, parent, name, text):
    cnode = self.dom.createElement(name)
    txt = self.dom.createTextNode(text)
    cnode.appendChild(txt)
    parent.appendChild(cnode)

  def add_clouds(self, atom, ideal_pts, nbrs_list):

    # create an atom node
    atom_node = self.dom.createElement("atom")
    self.add_text_node(atom_node, "chainID", atom.chainID)
    self.add_text_node(atom_node, "resSeq", "%d" % (atom.resSeq))
    self.add_text_node(atom_node, "iCode", atom.iCode)
    self.add_text_node(atom_node, "name", atom.name)

    for i in range(len(ideal_pts)):
      ideal_pt = ideal_pts[i]
      if(len(ideal_pt) == 0): continue

      # create a cloud node
      pt_node = self.dom.createElement("cloud")
      self.add_text_node(pt_node, "number", "%d" % (i + 1))
      self.add_text_node(pt_node, "position",
                         "%f %f %f" % (ideal_pt[0], ideal_pt[1], ideal_pt[2]))
      atom_node.appendChild(pt_node)
  
      # Add the neighbors' positions
      for nbr in nbrs_list[i]:
        nbr_node = self.dom.createElement("neighbor")
        self.add_text_node(nbr_node, "position", 
                           "%f %f %f" % (nbr.position[0], nbr.position[1], nbr.position[2]))
        pt_node.appendChild(nbr_node)

    # Add the node
    self.clouds.appendChild(atom_node)
    

  def toxml(self, xml_file):
    #print self.dom.toprettyxml(indent="  ", newl="")
    #print self.dom.toprettyxml(indent="  ")
    print self.dom.writexml(xml_file)

  def cleanup(self):
    self.dom.unlink()


if __name__ == "__main__":

  csv_fname = "/home/vanvoor4/data/new_sampling/pterins/dbase/2qx0_ph2_s.csv"
  prot_fname = "/home/vanvoor4/data/new_sampling/pterins/proteins/2qx0_ph2_p.pdb" 
  lig_fname = "/home/vanvoor4/data/new_sampling/pterins/ligands/2qx0_ph2_l.mol2"
  
  my_klas = test_xml(csv_fname)
  


