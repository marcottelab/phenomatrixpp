/* 
 * File:   adjacency_list.cpp
 * Author: jwoods
 * 
 * Created on July 15, 2009, 7:53 PM
 */

#include "adjacency_list.h"

// Helper function which gets the number of rows in some column of a cell_test_set.
// Note that cell test sets are ordered by column.
set<size_t> rows_in_column_of_cell_test_set(size_t j, const multimap<size_t,size_t>& cells) {
    set<size_t> ret;

    pair<multimap<size_t,size_t>::const_iterator,multimap<size_t,size_t>::const_iterator>
                range = cells.equal_range(j);

    for (multimap<size_t,size_t>::const_iterator it = range.first; it != range.second; ++it)
        ret.insert(it->second);

    return ret;
}

// Basically just union of two sets, but without all the copying if it's not
// necessary.
set<size_t> combine_rows_from_cell_test_set(set<size_t> jset, set<size_t> kset) {
    if (jset.size() == 0) return kset;
    if (kset.size() == 0) return jset;

    set<size_t> ret;
    set_union(jset.begin(), jset.end(), kset.begin(), kset.end(),
            insert_iterator<set<size_t> >(ret, ret.begin()));

    return ret;
}


// Third argument is a dummy variable, necessary so this can work as a drop-in
// replacement for UBLAS matrices.
// First argument is also a dummy variable.
adjacency_list::adjacency_list(size_t rows, size_t columns, bool weighted = false) :
    list_(columns, set<size_t>()),
    weights_(NULL)
{
    // This will be true when any oracle is used except mindist.
    if (weighted)
        weights_ = new matrix<double>(rows, columns);
}


adjacency_list::adjacency_list(const adjacency_list& orig) {
    cerr << "Error: Copy constructor not written yet." << endl;
    throw;
}


adjacency_list::~adjacency_list() {
    delete weights_;
}

const double adjacency_list::operator()(size_t i, size_t j) const {
    if (weights_ == NULL)
        // Convert bool to double
        return (double)(list_[j].find(i) != list_[j].end());
    else {
        // Weighted matrix; return the weights
        if (list_[j].find(i) != list_[j].end())    return (*weights_)(i,j);
        else                                       return 0.0;
    }
}
