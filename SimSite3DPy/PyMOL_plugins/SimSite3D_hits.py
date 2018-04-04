import os
import Tkinter
from Tkinter import *
import Pmw
import traceback
from SimSite3D import *

def __init__(self):
  self.menuBar.addmenuitem('Plugin', 'command', 'SimSite3D_hits',
                           label='SimSite3D hits',
                           command = lambda s=self : SimSite3D_hits(s))

class LabeledRealCounter(Pmw.Counter):
  
  def __init__(self, parent=None, minval=-10.0, maxval=10.0, **kw):

    # Define the validate options
    validate = {'validator':'real', 'min':minval, 'max':maxval}
    kw["datatype"] = "real"
    kw["entryfield_validate"] = validate
    kw["entryfield_value"] = -1.5
    kw["labelpos"] = "w"
    kw["increment"] = 0.01

    apply(Pmw.Counter.__init__, (self, parent), kw)

class LabeledIntCounter(Pmw.Counter):
  
  def __init__(self, parent=None, minval=0.0, maxval=1000, **kw):

    # Define the validate options
    validate = {'validator':'integer', 'min':minval, 'max':maxval}
    kw["orient"] = "horizontal"
    kw["datatype"] = "integer"
    kw["entryfield_validate"] = validate
    kw["entryfield_value"] = 20
    kw["labelpos"] = "w"
    kw["increment"] = 1

    apply(Pmw.Counter.__init__, (self, parent), kw)

class ScoreThreshold:
 
  def __init__(self, parent, minval=-10.0, maxval=10.0, **kw):
    self.dialog = Pmw.LabeledRealCounter(parent,
      label_text = "Select the maximum score you would like to see",
      minval=minval, maxval=maxval)
    self.dialog.withdraw()



class SimSite3D_hits:

  def set_hits_fname(self, fname):
    print "SimSite3D: set hits file to", fname
    self.hits_fname.setvalue(fname)

  def set_query_lig_fname(self, fname):
    self.query_lig_fname.setvalue(fname)

  def __init__(self, app):
    self.app = app
    parent = app.root
    self.parent = parent

    # Create the dialog
    self.dialog = Pmw.Dialog(parent,
                             buttons = ('Display', 'Exit'),
                             title = 'SimSite3D hits',
                             command = self.execute)
    self.dialog.withdraw()
    Pmw.setbusycursorattributes(self.dialog.component('hull'))

    w = Tkinter.Label(self.dialog.interior(), 
                      text = 'PyMOL visualization of SimSite3D hits\n' + \
                      'Jeffrey Van Voorst & Leslie Kuhn, 2008-2010\n',
                      background = 'black', foreground = 'white',)
    w.pack(expand = 1, fill = 'both', padx = 4, pady = 4)
    self.notebook = Pmw.NoteBook(self.dialog.interior())
    self.notebook.pack(fill='both',expand=1,padx=10,pady=10)

    def quickFileValidation(s):
      if s == '': return Pmw.PARTIAL
      elif os.path.isfile(s): return Pmw.OK
      elif os.path.exists(s): return Pmw.PARTIAL 
      else: return Pmw.PARTIAL
    
    # Set up the Main page
#    page = self.notebook.add('Main')
#    group = Pmw.Group(page,tag_text='Main options')
    page = self.notebook.add('Hits')
    group = Pmw.Group(page,tag_text='Load hits file')
    group.pack(fill = 'both', expand = 1, padx = 10, pady = 5)

    FileButton = PmwExtras.FileDialogButtonClassFactory 
    self.hits_fname = \
      Pmw.EntryField(group.interior(), labelpos = 'w',
                     label_pyclass = FileButton.get(self.set_hits_fname), 
                     validate = {'validator':quickFileValidation,}, 
                     label_text = 'SimSite3D hits file')
    self.query_lig_fname = \
      Pmw.EntryField(group.interior(), labelpos = 'w', 
                     label_pyclass = FileButton.get(self.set_query_lig_fname), 
                     validate = {'validator':quickFileValidation,}, 
                     label_text = 'Query ligand file') 

# enable these for now to have a working plugin
    self.query_C_color = \
      Pmw.EntryField(group.interior(), labelpos = 'w',
                     label_text = 'Query carbon color', 
                     value = 'green')
    self.hbond_stick_rad = \
      Pmw.EntryField(group.interior(), labelpos = 'w',
                     label_text = 'Hbond vector radius',
                     value=0.015, validate = {'validator':'real', 'min':0.0,},)


    self._score_counter = \
      LabeledRealCounter(group.interior(), label_text="Score threshold")
    self._max_num_hits_counter = \
      LabeledIntCounter(group.interior(), label_text="Max # hits to load")

    fields = (self.hits_fname, self.query_lig_fname, self.query_C_color,
              self.hbond_stick_rad, self._score_counter)
    fields = (self.hits_fname, self.query_lig_fname, self._score_counter, 
              self._max_num_hits_counter)
    Pmw.alignlabels(fields)
    for field in fields:
      field.pack(fill='x', expand=1, padx=20, pady=5)

    self.max_num_hits_per_pair = 1
    #self.max_num_hits_per_pair = 10
    load_hits = Tkinter.Button(page, text="Load", command=self.load_hits)
    load_hits.pack(expand=1)
    clear_hits = Tkinter.Button(page, text="Clear")
    clear_hits.pack(expand=1)

    # Show the initial dialog
    self.notebook.setnaturalsize()
    self.dialog.show()

  def load_hits(self):
    """
    Load the hits file from path stored in self.hits_fname
    
    1) If ....
    """
    #from SimSite3D import hits
    hits_fname = self.hits_fname.getvalue()
    try:
      self.hits_meta = hits(hits_fname, self.max_num_hits_per_pair)
    except IOError, (errno, strerror):
      msg = "Unable to open the file: %s\nerror(%s): %s" % \
        (hits_fname, errno, strerror), 
      errmsg = Pmw.MessageDialog(None, message_text=msg,
        title="Unable to open hits file", buttons=("Dismiss",), 
        defaultbutton="Dismiss")
      errmsg.command=errmsg.destroy

  def execute(self, result):
    print "result of execute is:", result

    if result == "Display":
      # Error here if self.hits_meta does not exist

      hits_meta = self.hits_meta
      hits = hits_meta.hits

      score_threshold = float(self._score_counter.getvalue())
      max_num_hits = int(self._max_num_hits_counter.getvalue())
      #max_num_hits = 200
      #score_threshold = -1.5
      if(len(hits) < max_num_hits): max_num_hits = len(hits)
      print "max num hits to display is:", max_num_hits
      score_threshold = min(hits[max_num_hits - 1][1], score_threshold)
      print "score threshold is:", score_threshold

      aligns = assess_alignments(hits, hits_meta.query_dir, hits_meta.query, 
                                 self.query_lig_fname.getvalue(),
                                 hits_meta.dbase_dir, hits_meta.ligs_dir,
                                 self.query_C_color.getvalue(), 
                                 hbond_stick_rad = \
                                   float(self.hbond_stick_rad.getvalue()),
                                 score_tol = score_threshold)


      x = ResultsTable(self.app, hits, hits_meta.query, aligns.q_site, 
                       score_threshold)
    else:
      self.dialog.withdraw()



class ResultsTable:

  def __init__(self, app, hits, query_id, q_site, score_tol = -1.5):
    parent = app.root
    self.parent = parent

    # Create the dialog
    self.dialog = Pmw.Dialog(parent,
                             buttons = ('Hide',),
                             title = 'SimSite3D Query Match Prints',
                             command = self.execute)
    self.dialog.withdraw()
    Pmw.setbusycursorattributes(self.dialog.component('hull'))

    fixedFont = Pmw.logicalfont('Fixed')
    my_label = "%s SimSite3D hits" % (query_id)
    self.st = Pmw.ScrolledText(self.dialog.interior(), # borderframe = 1,
                               labelpos = 'n', 
                               label_text=my_label,
                               columnheader = 1, rowheader = 1, 
                               rowcolumnheader = 1, usehullsize = 1, 
                               hull_width = 400, hull_height = 300, 
                               text_wrap='none', text_font = fixedFont, 
                               Header_font = fixedFont, 
                               Header_foreground = 'blue', rowheader_width = 5,
                               rowcolumnheader_width = 5, 
                               text_padx = 4, text_pady = 4, 
                               Header_padx = 4, rowheader_pady = 4,)
    self.st.pack(padx = 5, pady = 5, fill = 'both', expand = 1)

    # Create the header for the row headers
    self.st.component('rowcolumnheader').insert('end', 'state')

    # Create the column headers
    headerLine = 'site score  '
    # set id len based on max length of hit id
    id_len = 0
    for hit in hits:
      if(hit[0].endswith("mol2")): n = len(hit[0][:-13])
      else: n = len(hit[0])
      if(n > id_len): id_len = n
    headerLine += "ID%s   " % (" " * (id_len - 2))
    # set matchprint length based on hit matchprint leng
    N = len(hits[0][4])
    if(N > 10): headerLine += "Matchprint%s" % (" " * (N - 10))
    else: headerLine += "Matchprint"
    self.st.component('columnheader').insert('0.0', headerLine)

    self.st.tag_configure('A', background = 'red')
    self.st.tag_configure('D', background = 'blue')
    self.st.tag_configure('N', background = 'white')
    self.st.tag_configure('H', background = 'green')
    # add the query row
    print "id_len", id_len
    print "N", N
    self.st.insert('end', 
                   "            query%s%s\n" % (" " * (id_len - 2), "1" * N))
    self.st.component('rowheader').insert('end', "    1\n")

   # Add tags to color query matchprint
    pt_lists = { 'A':[], 'D':[], 'N':[], 'H':[] }
    state = 1
    #  Matchprints start at score (9) + 3 + id + 3 = id + 15
    q_pt_types = []
    pos = id_len + 15
    for pt in q_site.points: 
      tag1 = '%d.%d' % (state, pos)
      pos += 1
      tag2 = '%d.%d' % (state, pos)

      if(pt.tempFactor == 0.0): 
        pt_lists["A"].extend([tag1, tag2])
        q_pt_types.append("A")
      elif(pt.tempFactor == 50.0): 
        pt_lists["D"].extend([tag1, tag2])
        q_pt_types.append("D")
      elif(pt.tempFactor == 25.0): 
        pt_lists["N"].extend([tag1, tag2])
        q_pt_types.append("N")
      elif(pt.tempFactor == 100.0): 
        pt_lists["H"].extend([tag1, tag2])
        q_pt_types.append("H")
      else: print "bad tempFactor %f" % (pt.tempFactor)

    # Add hits rows
    state = 2
    for hit in hits:
      # Assume hits are sorted by score
      if(hit[1] > score_tol): break

      # The state holding the database object
      header = "%5d" % (state) 
      # The score of this match
      dataLine = "   %6.2f   " % (hit[1])
      # The database protein id
# dumb little hack

      if(hit[0].endswith(".mol2")): id = hit[0][:-13]
      else: id = hit[0]
      x = len(id)
      dataLine += "%s%s   " % (id, " " * (id_len - x))
      # The matchprint
      dataLine += hit[4]

      if(not hit == hits[-1]):
        dataLine += '\n' 
        header += '\n' 

      # Add tags to color query matchprint
      pos = id_len + 15
      for m, pt_type in zip(hit[4], q_pt_types):
        tag1 = '%d.%d' % (state, pos)
        pos += 1
        tag2 = '%d.%d' % (state, pos)
        if(m == "1"): pt_lists[pt_type].extend([tag1, tag2])
        
      # Add row 
      self.st.insert('end', dataLine)
      self.st.component('rowheader').insert('end', header)

      state += 1

    # Apply the tags
    for tag in pt_lists:
      # Need to check for the length since we may not have points in all of the 
      # categories
      if(len(pt_lists[tag])):
        apply(self.st.tag_add, (tag,) + tuple(pt_lists[tag]))

    # Prevent users' modifying text and headers
    self.st.configure(text_state = 'disabled', Header_state = 'disabled',)

    self.dialog.show()

  def execute(self, result):
    if result == "Hide": self.dialog.withdraw()
    else: self.dialog.destroy()


# Create demo in root window for testing.
if __name__ == '__main__':
    class App:
        def my_show(self,*args,**kwargs):
            pass
    app = App()
    app.root = Tkinter.Tk()
    Pmw.initialise(app.root)
    app.root.title('Some Title')

    widget = SimSite3D_hits(app)
    exitButton = Tkinter.Button(app.root, text = 'Exit2', command = app.root.destroy)
    exitButton.pack()
    app.root.mainloop()
