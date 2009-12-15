/* 
 * File:   adjacency_list.h
 * Author: jwoods
 *
 * Created on July 15, 2009, 7:53 PM
 */

#ifndef _ADJACENCY_LIST_H
#define	_ADJACENCY_LIST_H

#include "constants.h"
using std::set_intersection;
using std::set_union;
using std::insert_iterator;



// Helper function which gets the number of rows in some column of a cell_test_set.
// Note that cell test sets are ordered by column.
set<size_t> rows_in_column_of_cell_test_set(size_t j, const multimap<size_t,size_t>& cells);

// Basically just union of two sets, but without all the copying if it's not
// necessary.
set<size_t> combine_rows_from_cell_test_set(set<size_t> jset, set<size_t> kset);


// Goal is to duplicate the interface of boost UBLAS matrices as much as
// possible, so this can be a drop-in replacement.
// This particular adjacency list is bipartite--and asymmetric--thus the BGL
// (Boost) adjacency list class will not work
class adjacency_list {
public:
    adjacency_list(size_t rows, size_t columns, bool weighted);
    adjacency_list(const adjacency_list& orig);
    virtual ~adjacency_list();

    // Value of the association between column j and row i?
    const double operator()(size_t i, size_t j) const;

    bool has_edge(size_t i, size_t j) {
#ifdef DEBUG
        if (j < list_.size())
#endif
        if (list_[j].find(i) != list_[j].end())
            return true;
        else
            return false;
#ifdef DEBUG
        else {
            cerr << "Warning: has_edge called, j out of range." << endl;
            return false;
        }
#endif
    }

    double edge_weight(size_t i, size_t j) {
        if (!has_edge(i,j))
            return 0.0;
        
        if (weights_ != NULL)
            return (*weights_)(i,j);
        else
            // No weight matrix allocated; just return 0 or 1.
            return (double)(has_edge(i,j));
    }

    // If the edge does not exist, create it (weight 1). If the edge does exist,
    // add 'amount' to its weight.
    // Only add the edge if weights is not allocated.
    double add_to_edge(size_t i, size_t j, double amount = 1) {
        if (!has_edge(i,j)) {
            add_edge(i,j);
            if (weights_ != NULL)
                (*weights_)(i,j) = amount;
        } else if (weights_ != NULL) {
            (*weights_)(i,j) += amount;
            return (*weights_)(i,j);
        }

        // If weights is not allocated or there was no edge before this function
        // was called, just return a bool.
        return 1.0;
    }

    double divide_edge_by(size_t i, size_t j, double div) {
        if (has_edge(i,j))
            return ((*weights_)(i,j) /= div);
        else
            return 0.0;
    }

    // Set edge value and truth. Setting to 0 means no edge any longer.
    void set_edge(size_t i, size_t j, double amount) {
        if (amount > 0) {
            if (!has_edge(i,j))
                add_edge(i,j);
        } else {
            remove_edge(i,j);
        }

        // Set weight if allocated.
        if (weights_ != NULL)
            (*weights_)(i,j) = amount;
    }

    // Remove an association.
    void remove_edge(size_t i, size_t j) {
        list_[j].erase(i);
        if (weights_ != NULL)
            (*weights_)(i,j) = 0;
    }

    // Returns the intersection of columns j and k's sets
    set<size_t> common_edges(size_t j, size_t k) const {

        // Find the intersection.
        set<size_t> return_set;
        set_intersection(list_[j].begin(), list_[j].end(),
                list_[k].begin(), list_[k].end(),
                insert_iterator<set<size_t> >(return_set, return_set.begin()));

        return return_set;
    }

    // Calls the other common_edges function, removes from the result those rows
    // found in row_test_set.
    set<size_t> common_edges(size_t j, size_t k, const set<size_t>& row_test_set) const {
        set<size_t> return_set;
        set<size_t> intersection = common_edges(j,k);

        // Simple case requires one less set copy operation.
        if (row_test_set.size() == 0 || intersection.size() == 0) return intersection;

        // row-based cross-validation
        set_difference(intersection.begin(), intersection.end(),
                       row_test_set.begin(), row_test_set.end(),
                       insert_iterator<set<size_t> >(return_set, return_set.begin()));
        return return_set;
    }

    // Returns the cardinality of the intersection of columns j and k's sets.
    // row_test_set is an optional argument which excludes certain rows from the
    // calculation.
    size_t num_common_edges(size_t j, size_t k, const set<size_t>& row_test_set, const multimap<size_t,size_t>& cell_test_set) const {
        if (j == k) return list_[j].size(); // simple case
        if (row_test_set.size() == 0) {
            if (cell_test_set.size() == 0) {
                // no test set supplied.
                return common_edges(j,k).size(); // no test set supplied
            } else {
                // cell-based test_set supplied
                set<size_t> rows_in_j_or_k_test_set = combine_rows_from_cell_test_set(
                        rows_in_column_of_cell_test_set(j, cell_test_set),
                        rows_in_column_of_cell_test_set(k, cell_test_set)
                        );
                // pass the result as if it's a row_test_set
                return common_edges(j, k, rows_in_j_or_k_test_set).size();
            }
        } else {
            // row_test_set supplied
            return common_edges(j, k, row_test_set).size();
        }
    }
    size_t num_common_edges(size_t j, size_t k) const {
        if (j == k) return list_[j].size(); // simple case
        return common_edges(j,k).size(); // no test set supplied
    }

    // Get the number of rows associated with a given column.
    size_t num_rows(size_t j) const { return list_[j].size(); }
    size_t size1(size_t j) const { return num_rows(j); }

    size_t num_columns() const { return list_.size(); }
    size_t size2() const { return num_columns(); }
protected:


    // Deprecated; use add_to_edge. This is still used internally.
    // Add an association between column j and row i
    // Does not affect weight
    void add_edge(size_t i, size_t j) {
#ifdef DEBUG
        if (list_.size() <= j) {
            cerr << "Error: Adjacency list columns argument was too small." << endl;
            set<size_t> i_set; i_set.insert(i);
            list_.push_back( i_set );
        } else
#endif
            list_[j].insert(i);
    }
    
    // The vector is columns. Each set within the vector represents e.g., a
    // set of genes that are associated with a phenotype.
    vector< set<size_t> > list_;

    // Rows are the set items, columns are the vector items.
    // Dynamically allocated when requested.
    matrix<double>* weights_;
};

#endif	/* _ADJACENCY_LIST_H */

