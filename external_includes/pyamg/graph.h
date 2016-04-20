#ifndef GRAPH_H
#define GRAPH_H

#include <algorithm>
#include <stack>
#include <cassert>
#include <limits>

/*
 *  Compute a maximal independent set for a graph stored in CSR format
 *  using a greedy serial algorithm
 *
 *  Parameters
 *      num_rows   - number of rows in A (number of vertices)
 *      Ap[]       - CSR row pointer
 *      Aj[]       - CSR index array
 *      active     - value used for active vertices        (input)
 *       C         - value used to mark non-MIS vertices   (output)
 *       F         - value used to mark MIS vertices       (output)
 *      x[]        - state of each vertex
 *
 *  
 *  Returns:
 *      The number of nodes in the MIS.
 *
 *  Notes:
 *      Only the vertices with values with x[i] == active are considered 
 *      when determining the MIS.  Upon return, all active vertices will
 *      be assigned the value C or F depending on whether they are in the 
 *      MIS or not.
 *
 */
template<class IndexType, class IntegerType, class IndexContainerType, class IntegerContainerType>
IndexType maximal_independent_set_serial(const IndexType num_rows,
                                 const IndexContainerType& Ap, 
                                 const IndexContainerType& Aj, 
                                 const IntegerType  active,
                                 const IntegerType  C,
                                 const IntegerType  F,
                                        IntegerContainerType&  x)
{
    IndexType N = 0;
    
    for(IndexType i = 0; i < num_rows; i++){
        if(x[i] != active) continue;

        x[i] = C;
        N++;

        for(IndexType jj = Ap[i]; jj < Ap[i+1]; jj++){
            const IndexType j = Aj[jj];
            if(x[j] == active) {
                x[j] = F;
            }
        }

    }

    return N;
}

/*
 *  Compute a maximal independent set for a graph stored in CSR format
 *  using a variant of Luby's parallel MIS algorithm
 *
 *  Parameters
 *      num_rows   - number of rows in A (number of vertices)
 *      Ap[]       - CSR row pointer
 *      Aj[]       - CSR index array
 *      active     - value used for active vertices        (input)
 *       C         - value used to mark non-MIS vertices   (output)
 *       F         - value used to mark MIS vertices       (output)
 *      x[]        - state of each vertex
 *      y[]        - random values for each vertex
 *      max_iters  - maximum number of iterations 
 *                   by default max_iters=-1 and no limit 
 *                   is imposed
 *  
 *  Returns:
 *      The number of nodes in the MIS.
 *
 *  Notes:
 *      Only the vertices with values with x[i] == active are considered 
 *      when determining the MIS.  Upon return, all active vertices will
 *      be assigned the value C or F depending on whether they are in the 
 *      MIS or not.
 *  
 */
template<class IndexType, class IndexContainerType, class IntegerType, 
            class IntegerContainerType, class ValueType = double, class ValueContainerType>
IndexType maximal_independent_set_parallel(const IndexType num_rows,
                                   const IndexContainerType& Ap, 
                                   const IndexContainerType& Aj,
                                   const IntegerType active,
                                   const IntegerType  C,
                                   const IntegerType  F,
                                         IntegerContainerType& x,
                                   const ValueContainerType&  y,
                                   const IntegerType max_iters=-1)
{
    IndexType N = 0;
    IndexType num_iters = 0;

    bool active_nodes = true;

    while(active_nodes && (max_iters == -1 || num_iters < max_iters)){
        active_nodes = false;

        num_iters++;
        
        for(IndexType i = 0; i < num_rows; i++){
            const ValueType yi = y[i];

            if(x[i] != active) continue;
            
            const IndexType row_start = Ap[i];
            const IndexType row_end   = Ap[i+1];
    
            IndexType jj;

            for(jj = row_start; jj < row_end; jj++){
                const IndexType j  = Aj[jj];
                const IntegerType xj = x[j];

                if(xj == C) {
                    x[i] = F;                      //neighbor is MIS
                    break;  
                }
                
                if(xj == active){
                    const ValueType yj = y[j];
                    if(yj > yi)
                        break;                     //neighbor is larger 
                    else if (yj == yi && j > i)
                        break;                     //tie breaker goes to neighbor
                }
            }
   
            if(jj == row_end){
                for(jj = row_start; jj < row_end; jj++){
                    const IndexType j  = Aj[jj];
                    if(x[j] == active)
                        x[j] = F;
                }
                N++;
                x[i] = C;
            } else {
                active_nodes = true;
            }
        }
    } // end while
        
    //std::cout << std::endl << "Luby's finished in " << num_iters << " iterations " << std::endl;

    return N;
}

/*
 *  Compute a vertex coloring for a graph stored in CSR format.
 *
 *  The coloring is computed by removing maximal independent sets
 *  of vertices from the graph.  Specifically, at iteration i an
 *  independent set of the remaining subgraph is constructed and
 *  assigned color i.
 *
 *  Returns the K, the number of colors used in the coloring.
 *  On return x[i] \in [0,1, ..., K - 1] will contain the color
 *  of the i-th vertex.
 *
 */
template<class IndexType, class IndexContainerType, class IntegerContainerType>
IndexType vertex_coloring_mis(const IndexType num_rows,
                      const IndexContainerType& Ap, 
                      const IndexContainerType& Aj, 
                            IntegerContainerType&  x)
{
//    std::fill( x, x + num_rows, -1);
    std::fill( x.begin(), x.end(), -1);

    IndexType N = 0;
    IndexType K = 0;

    while(N < num_rows){
        N += maximal_independent_set_serial(num_rows,Ap,Aj,-1-K,K,-2-K,x);
        K++;
    }

    return K;
}


/*
 *  Applies the first fit heuristic to a graph coloring.
 *
 *  For each vertex with color K the vertex is assigned the *first* 
 *  available color such that no neighbor of the vertex has that
 *  color.  This heuristic is used to reduce the number of color used
 *  in the vertex coloring.
 *
 */
template<class IndexType, class IntegerType, class IndexContainerType, class IntegerContainerType>
void vertex_coloring_first_fit(const IndexType num_rows,
                               const IndexContainerType& Ap, 
                               const IndexContainerType& Aj, 
                                     IntegerContainerType& x,
                               const IntegerType K)
{
    for(IndexType i = 0; i < num_rows; i++){
        if(x[i] != K) continue;
        std::vector<bool> mask(K,false);
        for(IndexType jj = Ap[i]; jj < Ap[i+1]; jj++){
            const IndexType j = Aj[jj];
            if(  i == j  ) continue; //ignore diagonal
            if( x[j] < 0 ) continue; //ignore uncolored vertices
            mask[x[j]] = true;
        }
        x[i] = std::find(mask.begin(), mask.end(), false) - mask.begin();            
    }
}



/*
 * Compute a vertex coloring of a graph using the Jones-Plassmann algorithm
 *
 *  Parameters   
 *      num_rows   - number of rows in A (number of vertices)
 *      Ap[]       - CSR row pointer
 *      Aj[]       - CSR index array
 *      x[]        - color of each vertex
 *      y[]        - initial random values for each vertex
 *
 *  Notes:
 *      Arrays x and y will be overwritten
 *
 *  References:
 *      Mark IntegerType. Jones and Paul E. Plassmann
 *      A Parallel Graph Coloring Heuristic
 *      SIAM Journal on Scientific Computing 14:3 (1993) 654--669
 *      http://citeseer.ist.psu.edu/jones92parallel.html
 *      
 */
template<class IndexType, class IndexContainerType, class IntegerContainerType, class ValueContainerType, class IntegerType = int>
IndexType vertex_coloring_jones_plassmann(const IndexType num_rows,
                                  const IndexContainerType& Ap, 
                                  const IndexContainerType& Aj, 
                                         IntegerContainerType& x,
                                         ValueContainerType& y)
{
//    std::fill( x, x + num_rows, -1);
    std::fill( x.begin(), x.end(), -1);

    for(IndexType i = 0; i < num_rows; i++){
        y[i] += Ap[i+1] - Ap[i];
    }

    IndexType N = 0;
    IndexType K = 0; //iteration number

    while(N < num_rows){
        N += maximal_independent_set_parallel(num_rows,Ap,Aj,-1,static_cast<IntegerType>(K),-2,x,y,1);
        for(IndexType i = 0; i < num_rows; i++){
            if(x[i] == -2)
                x[i] = -1;
        }
        vertex_coloring_first_fit(num_rows,Ap,Aj,x,K);
        K++;
    }

//    return *std::max_element(x, x + num_rows);
    return *std::max_element(x.begin(), x.end());
}


/*
 * Compute a vertex coloring of a graph using the parallel 
 * Largest-Degree-First (LDF) algorithm
 *
 *  Parameters   
 *      num_rows   - number of rows in A (number of vertices)
 *      Ap[]       - CSR row pointer
 *      Aj[]       - CSR index array
 *      x[]        - color of each vertex
 *      y[]        - initial random values for each vertex
 *
 *   References:
 *     J. ValueType. Allwright and ValueType. Bordawekar and P. D. Coddington and K. Dincer and C. L. Martin
 *     A Comparison of Parallel Graph Coloring Algorithms
 *     DRAFT SCCS-666
 *     http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.45.4650
 *
 */
template<class IndexType, class ValueType = double, class IntegerType = int,
                class IndexContainerType, class IntegerContainerType, class ValueContainerType>
IndexType vertex_coloring_LDF(const IndexType num_rows,
                      const IndexContainerType& Ap, 
                      const IndexContainerType& Aj, 
                             IntegerContainerType&  x,
                      const ValueContainerType& y)
{
//    std::fill( x, x + num_rows, -1);
    std::fill( x.begin(), x.end(), -1);

    std::vector<ValueType> weights(num_rows);

    IndexType N = 0;
    IndexType K = 0; //iteration number

    while(N < num_rows){
        // weight is # edges in induced subgraph + random value
        for(IndexType i = 0; i < num_rows; i++){
            if(x[i] != -1) continue;
            IndexType num_neighbors = 0;
            for(IndexType jj = Ap[i]; jj < Ap[i+1]; jj++){
                IndexType j = Aj[jj];
                if(x[j] == -1 && i != j)
                    num_neighbors++;
            }
            weights[i] = y[i] + num_neighbors;
        }

//        N += maximal_independent_set_parallel(num_rows,Ap,Aj,-1,static_cast<IntegerType>(K),-2,x,&weights[0],1);
        N += maximal_independent_set_parallel(num_rows,Ap,Aj,-1,static_cast<IntegerType>(K),-2,x,weights,1);
        for(IndexType i = 0; i < num_rows; i++){
            if(x[i] == -2)
                x[i] = -1;
        }
        vertex_coloring_first_fit(num_rows,Ap,Aj,x,K);
        K++;
    }

//    return *std::max_element(x, x + num_rows);
    return *std::max_element(x.begin(), x.end());
}


/* 
 * Apply one iteration of Bellman-Ford iteration on a distance
 * graph stored in CSR format.
 *
 *  Parameters
 *      num_rows   - number of rows in A (number of vertices)
 *      Ap[]       - CSR row pointer
 *      Aj[]       - CSR index array
 *      Ax[]       - CSR data array (edge lengths)
 *      x[]        - (current) distance to nearest center
 *      y[]        - (current) index of nearest center
 *
 *  References:
 *      http://en.wikipedia.org/wiki/Bellman-Ford_algorithm
 */
template<class IndexType, class ValueType = double, class IndexContainerType, class ValueContainerType>
void bellman_ford(const IndexType num_rows,
                  const IndexContainerType& Ap, 
                  const IndexContainerType& Aj, 
                  const ValueContainerType& Ax,
                        ValueContainerType&  x,
                        IndexContainerType&  y)
{
    for(IndexType i = 0; i < num_rows; i++){
        ValueType xi = x[i];
        IndexType yi = y[i];
        for(IndexType jj = Ap[i]; jj < Ap[i+1]; jj++){
            const IndexType j = Aj[jj];
            const ValueType d = Ax[jj] + x[j];
            if(d < xi){
                xi = d;
                yi = y[j];
            }
        }
        x[i] = xi;
        y[i] = yi;
    }
}
 
    
/*
 * Perform Lloyd clustering on a distance graph
 *
 *  Parameters   
 *      num_rows       - number of rows in A (number of vertices)
 *      Ap[]           - CSR row pointer
 *      Aj[]           - CSR index array
 *      Ax[]           - CSR data array (edge lengths)
 *      x[num_rows]    - distance to nearest seed
 *      y[num_rows]    - cluster membership
 *      z[num_centers] - cluster centers
 *
 *  References
 *      Nathan Bell
 *      Algebraic Multigrid for Discrete Differential Forms
 *      PhD thesis (UIUC), August 2008
 *
 */
template<class IndexType, class ValueType = double, class IndexContainerType, class ValueContainerType>
void lloyd_cluster(const IndexType num_rows,
                   const IndexContainerType& Ap, 
                   const IndexContainerType& Aj, 
                   const ValueContainerType& Ax,
                   const IndexType num_seeds,
                         ValueContainerType&  x,
                         IndexContainerType&  y,
                         IndexContainerType&  z)
{
    for(IndexType i = 0; i < num_rows; i++){
        x[i] = std::numeric_limits<ValueType>::max();
        y[i] = -1;
    }
    for(IndexType i = 0; i < num_seeds; i++){
        IndexType seed = z[i];
        assert(seed >= 0 && seed < num_rows);
        x[seed] = 0;
        y[seed] = i;
    }

    std::vector<ValueType> old_distances(num_rows);

    // propagate distances outward
    do{
        std::copy(x, x+num_rows, old_distances.begin());
        bellman_ford(num_rows, Ap, Aj, Ax, x, y);
    } while ( !std::equal( x, x+num_rows, old_distances.begin() ) );

    //find boundaries
    for(IndexType i = 0; i < num_rows; i++){
        x[i] = std::numeric_limits<ValueType>::max();
    }
    for(IndexType i = 0; i < num_rows; i++){
        for(IndexType jj = Ap[i]; jj < Ap[i+1]; jj++){
            IndexType j = Aj[jj];
            if( y[i] != y[j] ){
                x[i] = 0;
                break;
            }
        }
    }

    // propagate distances inward
    do{
        std::copy(x, x+num_rows, old_distances.begin());
        bellman_ford(num_rows, Ap, Aj, Ax, x, y);
    } while ( !std::equal( x, x+num_rows, old_distances.begin() ) );


    // compute new seeds
    for(IndexType i = 0; i < num_rows; i++){
        const IndexType seed = y[i];

        if (seed == -1) //node belongs to no cluster
            continue;
        
        assert(seed >= 0 && seed < num_seeds);

        if( x[z[seed]] < x[i] )
            z[seed] = i;
    }
}




/*
 * Propagate (key,value) pairs across a graph in CSR format.
 *
 * Each vertex in the graph looks at all neighboring vertices
 * and selects the (key,value) pair such that the value is 
 * greater or equal to every other neighboring value.  If
 * two (key,value) pairs have the same value, the one with 
 * the higher index is chosen
 *
 * This method is used within a parallel MIS-k method to 
 * propagate the local maximia's information to neighboring
 * vertices at distance K > 1 away.
 *
 */
template<typename IndexType, typename ValueType = double, class IndexContainerType, class ValueContainerType>
void csr_propagate_max(const IndexType  num_rows,
                       const IndexContainerType&  Ap, 
                       const IndexContainerType&  Aj,
                       const IndexContainerType&  i_keys,
                             IndexContainerType&  o_keys,
                       const ValueContainerType&  i_vals,
                             ValueContainerType&  o_vals)
{
    for(IndexType i = 0; i < num_rows; i++){

        IndexType k_max = i_keys[i];
        ValueType v_max = i_vals[i];

        for(IndexType jj = Ap[i]; jj < Ap[i+1]; jj++){
            const IndexType j   = Aj[jj];
            const IndexType k_j = i_keys[j];
            const ValueType v_j = i_vals[j];

            if( k_j == k_max ) continue;
            if( v_j < v_max ) continue;
            if( v_j > v_max || k_j > k_max ){
                k_max = k_j;
                v_max = v_j;
            }
        }

        o_keys[i] = k_max;
        o_vals[i] = v_max;
    }
}

/*
 *  Compute a distance-k maximal independent set for a graph stored
 *  in CSR format using a parallel algorithm.  An MIS-k is a set of 
 *  vertices such that all vertices in the MIS-k are separated by a
 *  path of at least K+1 edges and no additional vertex can be added
 *  to the set without destroying this property.  A standard MIS
 *  is therefore a MIS-1.
 *
 *  Parameters
 *      num_rows   - number of rows in A (number of vertices)
 *      Ap[]       - CSR row pointer
 *      Aj[]       - CSR index array
 *      k          - minimum separation between MIS vertices
 *      x[]        - state of each vertex (1 if in the MIS, 0 otherwise)
 *      y[]        - random values used during parallel MIS algorithm 
 *      max_iters  - maximum number of iterations to use (default, no limit)
 *  
 */
template<class IndexType, class IntegerType, class ValueType = double, 
            class IndexContainerType, class ValueContainerType, class IntegerContainerType>
void maximal_independent_set_k_parallel(const IndexType num_rows,
                                        const IndexContainerType& Ap, 
                                        const IndexContainerType& Aj,
                                        const IndexType  k,
                                              IntegerContainerType& x,
                                        const ValueContainerType&  y,
                                        const IntegerType max_iters=-1)
{
    std::vector<bool> active(num_rows,true);

    std::vector<IndexType> i_keys(num_rows);
    std::vector<IndexType> o_keys(num_rows);
    std::vector<ValueType> i_vals(num_rows); 
    std::vector<ValueType> o_vals(num_rows); 

    for(IndexType i = 0; i < num_rows; i++){
        i_keys[i] = i;
        i_vals[i] = y[i];
        x[i] = 0;
    }

    for(IndexType iter = 0; max_iters == -1 || iter < max_iters; iter++){
        for(IndexType i = 0; i < k; i++){
            csr_propagate_max(num_rows, Ap, Aj, &(i_keys[0]), &(o_keys[0]), &(i_vals[0]), &(o_vals[0]));
            std::swap(i_keys, o_keys);
            std::swap(i_vals, o_vals);
        }

        for(IndexType i = 0; i < num_rows; i++){
            if( i_keys[i] == i && active[i]){
                x[i] = 1; // i is a MIS-k node
            } 
            
            i_keys[i] = i;
            i_vals[i] = x[i];
        }
       
        IndexType rank = 0;
        //while(rank < k && 2*(k - rank) > k){
        //    csr_propagate_max(num_rows, Ap, Aj, &(i_keys[0]), &(o_keys[0]), &(i_vals[0]), &(o_vals[0]));
        //    std::swap(i_keys, o_keys);
        //    std::swap(i_vals, o_vals);
        //    rank++;
        //}
        
        while(rank < k){
            csr_propagate_max(num_rows, Ap, Aj, &(i_keys[0]), &(o_keys[0]), &(i_vals[0]), &(o_vals[0]));
            std::swap(i_keys, o_keys);
            std::swap(i_vals, o_vals);
            rank++;
        }

        bool work_left = false;

        for(IndexType i = 0; i < num_rows; i++){
            if(i_vals[i] == 1){
                active[i] =  false;
                i_vals[i] = -1;
            } else {
                i_vals[i] = y[i];
                work_left = true;
            } 
            i_keys[i] = i;
        }

        if( !work_left )
            return;
    }

}

/*
 *  Compute a breadth first search of a graph in CSR format
 *  beginning at a given seed vertex.
 *
 *  Parameters
 *      num_rows         - number of rows in A (number of vertices)
 *      Ap[]             - CSR row pointer
 *      Aj[]             - CSR index array
 *      order[num_rows]  - records the order in which vertices were searched
 *      level[num_rows]  - records the level set of the searched vertices (i.e. the minimum distance to the seed)
 *
 *  Notes:
 *      The values of the level must be initialized to -1
 *
 */
template <class IndexType, class IndexContainerType, class IntegerType = int, class IntegerContainerType>
void breadth_first_search(const IndexContainerType& Ap, 
                          const IndexContainerType& Aj,
                          const IndexType seed,
                                IndexContainerType& order,
                                IntegerContainerType& level)
{
    // initialize seed
    order[0]    = seed;
    level[seed] = 0;
   
    IntegerType N = 1;
    IntegerType level_begin = 0;
    IntegerType level_end   = N;

    IntegerType current_level = 1;

    while(level_begin < level_end){
        // for each node of the last level
        for(IndexType ii = level_begin; ii < level_end; ii++){
            const IndexType i = order[ii];

            // add all unmarked neighbors to the queue
            for(IndexType jj = Ap[i]; jj < Ap[i+1]; jj++){
                const IndexType j = Aj[jj];
                if(level[j] == -1){
                    order[N] = j;
                    level[j] = current_level;
                    N++;
                }
            }
        }

        level_begin = level_end;
        level_end   = N;
        current_level++;
    }

}


/*
 *  Compute the connected components of a graph stored in CSR format.
 *
 *  Vertices belonging to each component are marked with a unique integer
 *  in the range [0,K), where K is the number of components.
 *
 *  Parameters
 *      num_rows             - number of rows in A (number of vertices)
 *      Ap[]                 - CSR row pointer
 *      Aj[]                 - CSR index array
 *      components[num_rows] - component labels
 *
 */
template <class IndexType, class IndexContainerType, class IntegerContainerType>
IndexType connected_components(const IndexType num_nodes,
                       const IndexContainerType& Ap, 
                       const IndexContainerType& Aj,
                             IntegerContainerType& components)
{
//    std::fill(components, components + num_nodes, -1);
    std::fill(components.begin(), components.end(), -1);
    std::stack<IndexType> DFS;
    IndexType component = 0;

    for(IndexType i = 0; i < num_nodes; i++)
    {
        if(components[i] == -1)
        {
            DFS.push(i);
            components[i] = component;

            while (!DFS.empty())
            {
                IndexType top = DFS.top();
                DFS.pop();
    
                for(IndexType jj = Ap[top]; jj < Ap[top + 1]; jj++){
                    const IndexType j = Aj[jj];
                    if(components[j] == -1){
                        DFS.push(j);
                        components[j] = component;
                    }
                }
            }

            component++;
        }
    }

    return component;
}

#endif
