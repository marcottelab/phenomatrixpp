/* 
 * File:   genephene.h
 * Author: jwoods
 *
 * Matrix class with some species' genes as rows and one or more species'
 * phenotypes as columns.
 *
 * This is a templated class. T is the data type stored in it.
 *
 * Created on June 19, 2009, 4:04 PM
 */

#ifndef _GENEPHENE_H
#define	_GENEPHENE_H

#include "utilities.h"
#include "adjacency_list.h"




class genephene {
public:

    // NEW CONSTRUCTOR:
    // Doesn't read a file; file must be read in advance. 'identifier' controls how the distance and common matrices
    // are saved in the filesystem.
    genephene(double (*distance_fn)(size_t,size_t,size_t,size_t),
              const set<gene_id_t>& genes,
              size_t num_phenes,
              const multimap<gene_id_t,phene_id_t>& assn,
              const string& identifier,
              const set<phene_id_t>& local_phenes);
    // NOTE: phene_filename is just the phenotypes for the source species. phenes is the total
    // matrix size, which will be source species + human phenotypes.
    genephene(double (*distance_fn)(size_t,size_t,size_t,size_t), size_t phenes, string gene_filename, string phene_filename, string assn_filename);
    genephene(double (*distance_fn)(size_t,size_t,size_t,size_t), const set<gene_id_t>& genes, const set<phene_id_t>& phenes, bool weighted_predictions);
    genephene(const genephene& orig);
    virtual ~genephene();
    
    // Load or save the contents of the adjacency list.
    void load(const string& filename);
    void save(const string& filename);

    void load_genes(const string& filename);

    bool save_common_items(string suffix = "") const;
    bool save_distances(string suffix, string distance_fun) const;

    // Attempt to load matrices -- return true if not possible.
    bool load_common_items(string suffix = "") const;
    bool load_distances(string suffix, string distance_fun) const;

    const dist_t operator()(const gene_id_t& g, const phene_id_t& p) const;

    // Accessors for contents of phene_distance and _common_items, by phenotype ID.
    const dist_t distance(uint p1, uint p2) const;
    const size_t common_items(uint p1, uint p2) const;

    // Sometimes we want to know if this matrix even has a phenotype before we,
    // for example, call nearest.
    bool has_phene(const phene_id_t& p) const {
        return (phene_to_column.find(p) != phene_to_column.end());
    }
    // Same thing for genes.
    bool has_gene(const gene_id_t& g) const {
        return (gene_to_row.find(g) != gene_to_row.end());
    }

    void add_association(const gene_id_t& g, const phene_id_t& p) {
        gene_by_phene.add_to_edge(gene_to_row[g], phene_to_column[p]);
    }

    bool has_association(const gene_id_t& g, const phene_id_t& p) {
        return gene_by_phene.has_edge(gene_to_row[g], phene_to_column[p]);
    }

    dist_t add_to_association(const gene_id_t& g, const phene_id_t& p, dist_t amt = 1) {
        return gene_by_phene.add_to_edge(gene_to_row[g], phene_to_column[p], amt);
    }

    dist_t divide_association_by(const gene_id_t& g, const phene_id_t& p, const dist_t& div) {
        return gene_by_phene.divide_edge_by(gene_to_row[g], phene_to_column[p], div);
    }

    void calculate_distances(
            const set<gene_id_t>&                   r_test_set,
            const multimap<gene_id_t, phene_id_t>&  c_test_set,
            string                                  suffix,
            string                                  distance_fun);

    // Size functions.
    size_t num_phenes() const { return _num_phenes; }
    size_t num_genes() const  { return _num_genes; }
    size_t num_genes_with_phene(uint p) const;

    // Get the genes in this matrix in the form of a vector
    const vector<gene_id_t>& gene_vector() const { return row_to_gene; }
    set<gene_id_t> gene_set() const;
    const vector<phene_id_t>& phene_vector() const { return column_to_phene; }
    set<phene_id_t> phene_set() const;

    // Get species information from the file.
    string species1() const { return species1_; }
    void species1(const string& s1) { species1_ = s1; }
    string species2() const { return species2_; }
    void species2(const string& s2) { species2_ = s2; }

    // Get the nearest phenotype to some other phenotype
    pair<phene_id_t,dist_t> nearest(phene_id_t p) const;

    // Return all of the phenotypes with < upperbound distance from p.
    multiset<dist_phene_t>* sorting(phene_id_t p, dist_t upperbound = 1) const;
    void add_to_sorting(dist_phene_multiset_t& to, phene_id_t p, size_t k, dist_t upperbound = 1) const;

    void check_species() const;

    set<phene_id_t> local_phenes() const;

    phene_id_t col_to_phene(const size_t& col) const { return column_to_phene[col]; }

    // Gives an error if translation fails.
    size_t safely_translate_phenotype_to_column(const phene_id_t& p) const;

    void test_distance_function() const;

protected:
    size_t _num_phenes;
    size_t _num_genes;
    
    // Map DB accessions to rows/columns
    unordered_map<gene_id_t, size_t>  gene_to_row;
    unordered_map<phene_id_t,size_t> phene_to_column;

    // Map rows/columns to DB accessions
    vector<gene_id_t>          row_to_gene;
    vector<phene_id_t>         column_to_phene;

    // Map genes to phenotypes (genes are rows, phenotypes are columns),
    // most likely for just one species pair.
    adjacency_list gene_by_phene;

    // Allow allocation of a matrix which stores distances between phenotypes.
    // This stores them as 1 - the distance.
    symmetric_matrix<dist_t>* phene_distance;
    
    // Also allow allocation for certain distance functions of a common_items
    // matrix, keeping track of the cardinality of set intersections of different
    // columns.
    compressed_matrix<size_t>* _common_items;

    // Phenotypes which only belong to this species.
    set<size_t> local_columns;

    string original_filename;
    string species1_;
    string species2_;

    // Distance function functor to use for calculations.
    double (*distance_function)(size_t, size_t, size_t, size_t);

    set<size_t> test_set_genes_to_rows(const set<gene_id_t>& test_set) const;
    multimap<size_t,size_t> test_set_genes_and_phenes_to_cells(const multimap<gene_id_t, phene_id_t>& test_set) const;

    // Called by calculate_distances() if using Hypergeometric.
    void calculate_common_items(const set<size_t>& row_test_set, const multimap<size_t,size_t>& cell_test_set, string suffix = "");

    // Filenames for output files.
    string save_common_items_filename(string suffix = "") const {
        if (suffix == original_filename) return suffix + ".common";
        else return original_filename + "." + suffix + ".common";
    }
    string save_pcommon_items_filename(string suffix = "") const {
        if (suffix == original_filename) return suffix + ".pcommon";
        else return original_filename + "." + suffix + ".pcommon";
    }
    string save_pdistances_filename(string suffix = "", string distance_fun="hypergeometric") const {
        if (suffix == original_filename) return suffix + "." + distance_fun + ".pdistances";
        else return original_filename + "." + suffix + "." + distance_fun + ".pdistances";
    }
    string save_distances_filename(string suffix = "", string distance_fun="hypergeometric") const {
        if (suffix == original_filename) return suffix + "." + distance_fun + ".distances";
        else return original_filename + "." + suffix + "." + distance_fun + ".distances";
    }
    string save_col_to_phene_filename(string suffix = "") const {
        if (suffix == original_filename) return suffix + ".col_to_phene";
        else return original_filename + "." + suffix + ".col_to_phene";
    }
};

typedef shared_ptr<genephene> genephene_ptr;
typedef genephene phenomatrix;
typedef unordered_map<species_id_t, genephene_ptr> sp_genephene_map_t;




template <typename M, typename T>
void save_upper_matrix(const string& filename, const M& m, T default_value, const genephene& gp, bool phenotypes=false) {
    ofstream out(filename.c_str());
    out.precision(dbprec1::digits10);
    out.setf( std::ios::scientific,std::ios::floatfield );
    
    for (size_t i = 0; i < m.size1(); ++i) {
        for (size_t j = i; j < m.size2(); ++j) {
            if (m(i,j) != default_value) {
                if (phenotypes)  out << gp.col_to_phene(i) << '\t' << gp.col_to_phene(j) << '\t' << m(i,j) << endl;
                else             out << i << '\t' << j << '\t' << m(i,j) << endl;
            }
                
        }
    }

    out.close();
}



template <typename M, typename T>
void save_matrix(const string& filename, const M& m, T default_value) {
    ofstream out(filename.c_str());
    out.precision(dbprec1::digits10);
    out.setf(std::ios::scientific,std::ios::floatfield);

    for (size_t i = 0; i < m.size1(); ++i) {
        for (size_t j = 0; j < m.size2(); ++j) {
            if (m(i,j) != default_value)
                out << i << '\t' << j << '\t' << m(i,j) << endl;
        }
    }

    out.close();
}

template <typename M, typename T>
void save_upper_prob_matrix(const string& filename, const M& m, T default_value, const genephene& gp, bool phenotypes = false, T threshold = 0) {
    ofstream out(filename.c_str());
    out.precision(dbprec1::digits10);
    out.setf(std::ios::scientific,std::ios::floatfield);
    
    // We don't want to read certain values, e.g., when distance is so close to
    // 1.0 that it won't matter.
    T thresh_value;
    if (are_equal<T>(default_value, 1)) {
        thresh_value = default_value - threshold;
    } else if (are_equal<T>(default_value, 0)) {
        thresh_value = default_value + threshold;
    } else {
        cerr << "Error: threshold needs to be 0 or 1." << endl;
        throw;
    }

    out.precision(dbprec1::digits10);
    out.setf(std::ios::scientific,std::ios::floatfield);


    for (size_t i = 0; i < m.size1(); ++i) {
        for (size_t j = i; j < m.size2(); ++j) {

            T val = m(i,j);

            if ((are_equal<T>(default_value, 0) && val > thresh_value) || (are_equal<T>(default_value, 1) && val < thresh_value)) {
                if (phenotypes) out << gp.col_to_phene(i) << '\t' << gp.col_to_phene(j) << '\t' << val << endl;
                else            out << i << '\t' << j << '\t' << val << endl;
            }
            // Otherwise, don't write.

        }
    }

    out.close();
}


// Note that this function is not responsible for initializing positions that
// are not in the file.
template <typename M, typename T>
void load_matrix(const string& filename, M& m) {
    ifstream in(filename.c_str());

    size_t i = 0, j = 0;
    T val = 0;

    in >> i;
    in >> j;
    in >> val;


    while (in) {
        // Insert in matrix
        m(i,j) = val;

        in >> i;
        in >> j;
        in >> val;
    }

    in.close();
}


template <typename M, typename T>
void init_load_upper_matrix(const string& filename, M& m, T initial_value) {
    // Matrix must be square, at least.
    assert(m.size1() == m.size2());
    
    for (size_t i = 0; i < m.size1(); ++i) {
        for (size_t j = i; j < m.size2(); ++j) {
            m(i,j) = initial_value;
        }
    }

    load_matrix<M,T>(filename, m);
}


template <typename M, typename T>
void init_load_matrix(const string& filename, M& m, T initial_value) {
    for (size_t i = 0; i < m.size1(); ++i) {
        for (size_t j = 0; j < m.size2(); ++j) {
            m(i,j) = initial_value;
        }
    }

    load_matrix<M,T>(filename, m);
}


// Note that this function is not responsible for initializing positions that
// are not in the file.
template <typename M, typename T>
void load_prob_matrix(const string& filename, M& m) {
    ifstream in(filename.c_str());

    size_t i = 0, j = 0;
    T val = 1;

    in >> i;
    in >> j;
    in >> val;


    while (in) {
        // Insert in matrix
        m(i,j) = val;

        // Next read.
        in >> i;
        in >> j;
        in >> val;
    }

    in.close();
}



// This function is slower than necessary. Skip it if you're using a compressed
// matrix and your initial_value is 0. If your initial_value is 1, skip it but
// store values inverted (1 - x instead of x).
// If you skip it, call load_prob_matrix instead.
template <typename M, typename T>
void init_load_upper_prob_matrix(const string& filename, M& m, T initial_value) {
    for (size_t i = 0; i < m.size1(); ++i) {
        for (size_t j = i; j < m.size2(); ++j) {
            m(i,j) = initial_value;
        }
    }
    load_prob_matrix<M,T>(filename, m);
}


#endif	/* _GENEPHENE_H */

