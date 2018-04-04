/******************************************************************************
 * Copyright (c) 2006,2007, Michigan State University (MSU) Board of Trustees.
 *   All rights reserved.
 *
 * This file is part of the ASCbase Software project.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *
 * Authors: Jeffrey Van Voorst, vanvoor4@msu.edu
 *          Leslie Kuhn, Ph.D., KuhnL@msu.edu 
 *****************************************************************************/

/*
 * $Source: /psa/share/repository/pfizer_proj/src/basics/simple_graph.H,v $
 * $Revision: 1.2 $
 * $Author: vanvoor4 $
 * $Date: 2007-08-21 15:55:52 $
 * 
 * $Log: not supported by cvs2svn $
 * Revision 1.1  2007/06/07 20:50:12  vanvoor4
 * initial checkin
 *
 * 
 * 
 */

#ifndef _SIMPLE_GRAPH_HEADER_FILE
#define _SIMPLE_GRAPH_HEADER_FILE

#include <types.H>
#include <deque>
#include <vector>

namespace ASCbase{
  // Taken from "A new way to enumerate cycles in graph" in 
  // IEEE Proc of the Advanced International Conference on Telecom and ... 
  // by Liu, Hongbo and Wang, Jiaxin
  template<typename _Tp> class simple_graph{

  public:
    typedef typename std::vector<_Tp> _Tp_vec;
    typedef typename _Tp_vec::const_iterator _Tp_const_iter;
    typedef typename std::vector<_Tp_const_iter> adj_nodes;
    typedef typename adj_nodes::const_iterator adj_node_iter;
    typedef typename std::vector<adj_nodes> adj_list;
    typedef typename std::deque<_Tp_const_iter> path;
    typedef typename path::iterator path_iter;
    typedef typename path::const_iterator path_citer;
    typedef typename std::deque<path> path_queue;
    typedef typename path_queue::iterator path_queue_iter;
  
    simple_graph(const std::vector<_Tp>& V, const adj_list& E, 
                 uint min_len_in = 5, uint max_len_in = 6)
    {
      max_cycle_len = max_len_in;
      min_cycle_len = min_len_in;
      const _Tp_const_iter V_beg = V.begin(); 
      const _Tp_const_iter V_end = V.end(); 

      // Initialize the Q with paths of length zero 
      path_queue Q(V.size());
      path_queue_iter q = Q.begin();
      for(_Tp_const_iter v = V_beg; v < V_end; ++v, ++q) q->push_back(v);

      while(Q.size()){
        path& P = Q.front();    
	const _Tp_const_iter head = P.front();
	const _Tp_const_iter tail = P.back();

        // Check if there is an edge linking tail to head
        const adj_nodes& edges = E[head - V_beg];
	adj_node_iter hnbr;
	for(hnbr = edges.begin(); hnbr < edges.end(); ++hnbr){
          if(tail == *hnbr){
	    if(P.size() >= min_cycle_len && P.size() <= max_cycle_len) 
	      cycles.push_back(P);
            break;
	  }
	}

	//
	// add all edges in adj_list for tail which have order > order(head)
	// to Q (1 path for each) such that each added vertex is not in P
	//
        if(P.size() < max_cycle_len){
	  adj_node_iter tnbr;
	  const adj_nodes& tail_edges = E[tail - V_beg];
	  // Find nodes adjacent to tail with order > order(head)
	  for(tnbr = tail_edges.begin(); tnbr < tail_edges.end(); ++tnbr)
            if(head < *tnbr){
	      // Check if this node is already included in the current path
              path_iter p;
              for(p = P.begin(); p < P.end() && *tnbr != *p; ++p);
              if(p != P.end()) continue;
  
	      // Add path P to Q and append the node *tnbr to the path P -- do 
	      // not push back onto P itself; doing so will mess all iters to P
              Q.push_back(P);
	      Q.back().push_back(*tnbr);
	    }
        }
        Q.pop_front();
      }

      remove_duplicate_cycles();
    }

    //! How many cycles with length between min and max cycle length are in 
    //! the graph?
    inline unsigned int num_cycles() const { return cycles.size(); }

    //! List the numbers associated with the cycles of the graph that contain 
    //! the given node 
    inline std::vector<uint> cycles_that_contain(const _Tp_const_iter node)
    {
      std::vector<uint> cycle_nums;
      for(path_queue_iter c = cycles.begin(); c < cycles.end(); ++c)
        for(path_iter p = c->begin(); p < c->end(); ++p)
          if(*p == node){
            cycle_nums.push_back(c - cycles.begin());
            break;
          }
      return cycle_nums;
    }

    inline void get_cycle(uint cycle_num, adj_nodes* vertices)
    {
      path_queue_iter c = cycles.begin() + cycle_num;
      vertices->resize(c->size());
      std::copy(c->begin(), c->end(), vertices->begin());
    }

  private:
    uint max_cycle_len;  //!< Maximum length of enumerated cycles 
    uint min_cycle_len;  //!< Minimum length of enumerated cycles 
    path_queue cycles;   //!< Enumerated cycles (nodes)

    //! Both a "forward" and "backward" path is enumerated for each cycle --
    //! we only want one path for each cycle.
    void remove_duplicate_cycles()
    {
      std::vector<bool> dup(cycles.size(), false);
      for(uint i = 0; i < cycles.size(); ++i){
        if(dup[i]) continue;

        for(uint j = i+1; j < cycles.size(); ++j){
          if(dup[j]) continue;

          uint a;
          for(a = 0; a < cycles[i].size(); ++a){
            uint b;
            for(b = 0; b < cycles[j].size(); ++b){
              if(cycles[i][a] == cycles[j][b]) break;
            }
            if(b >= cycles[j].size()) break;
          }
          if(a >= cycles[i].size()){
            dup[j] = true;
            cycles[j].clear();
          }
        }
      }
      
      for(path_queue_iter c = cycles.begin(); c < cycles.end(); ++c){
        if(c->size() == 0){
          cycles.erase(c);
          c = cycles.begin();
        }
      }
    }

  };
}
#endif
