from xml.dom import minidom

############################## 
class score_info:

  def __init__(self, mean = 0.0, stdev = -1.0, N = -1):
    self.mean = mean
    self.stdev = stdev
    self.N = N
############################## 

############################## 
class norm_stats:
 
  ############################## 
  def __init__(self, xml_fname):
    xmldoc = minidom.parse(xml_fname)

    self.norm_stats = {}
    if(xmldoc.firstChild == None):
      print "Unable to parse the file: " + xml_in
  
    head_node = xmldoc.childNodes[0]
    if(not head_node.tagName == "datasets"):
      print "Expected to find datasets as the head node tag; found " + \
        head_node.tagName

    self.norm_stats = \
      self.process_datasets(xmldoc.getElementsByTagName("dataset"))

    xmldoc.unlink()
  ############################## 

  ############################## 
  def process_datasets(self, datasets):
    dset_vals = {}
    for dataset in datasets:
      dataset_name = dataset.attributes.getNamedItem("name").firstChild.data
      dset_vals[dataset_name] = \
        self.process_queries(dataset.getElementsByTagName("query"))
    return dset_vals
  ############################## 
  
  ############################## 
  def process_queries(self, queries):
    q_vals = {}
    for query in queries:
      query_name = query.attributes.getNamedItem("name").firstChild.data
      q_vals[query_name] = self.process_SFs(query.getElementsByTagName("SF"))
    return q_vals
  ############################## 
  
  ############################## 
  def process_SFs(self, SFs):
    SFs_stats = {}
    for SF in SFs:
      tmp_info = score_info()

      for node in SF.childNodes:
        if(not node.nodeType == node.ELEMENT_NODE):
          continue
  
        if(node.hasChildNodes()):
          if(node.tagName == "mean_score"):
            tmp_info.mean = float(node.firstChild.data)
          elif(node.tagName == "stdev_score"):
            tmp_info.stdev = float(node.firstChild.data)
          elif(node.tagName == "number_of_scores"):
            tmp_info.N = float(node.firstChild.data)
      
      SF_name = SF.attributes.getNamedItem("name").firstChild.data
      SFs_stats[SF_name] = tmp_info
    return SFs_stats
  ############################## 

############################## 
class scale_factors:
 
  ############################## 
  def __init__(self, data_dir, tag=None):
    self.scale_factors = {}

    from os import listdir
    for csv_fname in listdir(data_dir):
      if(not csv_fname.endswith("_s.csv")):
        continue

      if(not tag == None):
        q_id = csv_fname[:csv_fname.find(tag) - 1]
      else:
        q_id = csv_fname[:-6]
  
      csv_file = file(data_dir + "/" + csv_fname, "r")
      grab = False
      for line in csv_file:
        if(grab):
          self.scale_factors[q_id] = [ float(s) for s in line.split() ]
          break;
  
        if(line.startswith("<SimSite3D features scales>")):
          grab = True
      csv_file.close()


#my_stats = norm_stats("norm_stats.xml")
#for dataset in my_stats.norm_stats:
#  print dataset
#  for query_name in my_stats.norm_stats[dataset]:
#    print query_name
#    query = my_stats.norm_stats[dataset][query_name]
#    for SF in query:
#      print SF, query[SF].mean, query[SF].stdev
