/* 
 * File:   genephene.cpp
 * Author: jwoods
 * 
 * Created on June 19, 2009, 4:04 PM
 */

#include "genephene.h"


genephene::genephene(
        double  (*distance_fn)(size_t,size_t,size_t,size_t),
        const   set<gene_id_t>& genes,
        size_t  num_phenes,
        const   multimap<gene_id_t,phene_id_t>& assn,
        const   string& identifier,
        const   set<phene_id_t>& local_phenes
)
        
        : _num_phenes( num_phenes ),
        _num_genes( genes.size() ),
        gene_to_row( genes.size() ),
        phene_to_column( num_phenes ),
        gene_by_phene( genes.size(), num_phenes, false ),
        phene_distance(NULL),
        _common_items(NULL),
        original_filename(identifier),
        distance_function(distance_fn)
{

    // Get all of the phenes that were passed in.
    set<phene_id_t> phenes;
    for (multimap<gene_id_t,phene_id_t>::const_iterator it = assn.begin(); it != assn.end(); ++it)
        phenes.insert(it->second);

    // Build the conversion tables for genes and phenes.
    for (set<gene_id_t>::iterator it = genes.begin(); it != genes.end(); ++it) {
        row_to_gene.push_back(*it);
        gene_to_row[*it] = row_to_gene.size() - 1;
    }
    for (set<phene_id_t>::iterator it = phenes.begin(); it != phenes.end(); ++it) {
        column_to_phene.push_back(*it);
        phene_to_column[*it] = column_to_phene.size() - 1;
    }

    size_t throwcount = 0;

    // Add to local_columns.
    for (set<phene_id_t>::const_iterator it = local_phenes.begin(); it != local_phenes.end(); ++it) {
        
        // Check that each local_phene is in phene_to_column
        unordered_map<phene_id_t,size_t>::iterator pc = phene_to_column.find(*it);
        if (pc == phene_to_column.end()) {
            if (throwcount < 10)
                cerr << "Error: Phenotype " << *it << " in local_phenes not found in phenomatrix." << endl;
            throwcount++;
        } else local_columns.insert(pc->second);
        
    }

    if (throwcount > 0) throw;

    // Finally, add associations to the matrix.
    for (multimap<gene_id_t,phene_id_t>::const_iterator it = assn.begin(); it != assn.end(); ++it) {
        add_association(it->first, it->second);
    }

    // ERROR CHECK:
    // Now check that every column contains at least two genes.

    for (size_t i = 0; i < column_to_phene.size(); ++i) {
        // Don't bother with phenotypes which are too small for the destination species.
        if ((local_columns.size() == 0 || local_columns.find(i) != local_columns.end()) && gene_by_phene.num_rows(i) < 2) {
            cerr << "Error: column " << i << " (phene " << column_to_phene[i] << ") has " << gene_by_phene.num_rows(i) << " genes." << endl;
            throwcount++;
        }
    }

    if (throwcount > 0) {
        cerr << "Total errors: " << throwcount << "; raising exception..." << endl;
        cerr << "local_columns" << endl;
        for (set<size_t>::iterator i = local_columns.begin(); i != local_columns.end(); ++i)
            cerr << *i << "\t";
        cerr << endl;
        throw;
    }
    
}


genephene::genephene(
        double  (*distance_fn)(size_t,size_t,size_t,size_t),
        const set<gene_id_t>& genes,
        const set<phene_id_t>& phenes,
        bool weighted_predictions = false)
: _num_phenes(phenes.size()),
  _num_genes(genes.size()),
  gene_to_row(genes.size()),
  phene_to_column(phenes.size()),
  gene_by_phene(genes.size(), phenes.size(), weighted_predictions),
  phene_distance(NULL),
  _common_items(NULL),
  distance_function(distance_fn)
{

    // Create the two mappings from gene to row and phenotype to column (and vice-versa)
    for (set<gene_id_t>::const_iterator i = genes.begin(); i != genes.end(); ++i) {
        row_to_gene.push_back(*i);
        gene_to_row[*i] = row_to_gene.size() - 1;
    }

    for (set<phene_id_t>::const_iterator i = phenes.begin(); i != phenes.end(); ++i) {
        column_to_phene.push_back(*i);
        phene_to_column[*i] = column_to_phene.size() - 1;
    }


}


void genephene::save(const string& filename) {
    cerr << "Need to finish writing save()" << endl;
    throw;
}


void genephene::load_genes(const string& filename) {
    ifstream fin(filename.c_str());

    // Ignore header 'gene'
    string header;
    getline(fin, header);

    gene_id_t gene = GENE_NONE;
    fin >> gene;
    while (fin) {
        // Insert gene in the look-up tables.
        row_to_gene.push_back(gene);
        gene_to_row[gene] = row_to_gene.size() - 1;

        fin >> gene;
    }
}


void genephene::load(const string& filename) {
    ifstream fin(filename.c_str());

    string header_line1, header_line2;
    cerr << "header_line1 uncoupled" << endl; throw;
    // Ignore header lines
    getline(fin, header_line1);
    cout << "Ignoring line:" << header_line1 << endl;
    getline(fin, header_line2);
    cout << "Ignoring line:" << header_line2 << endl;

    //getline(fin, line);

    gene_id_t gene = GENE_NONE;
    uint phene = 0;
    // Parse the line -- priming read
    read_gene_phene_line(fin, gene, phene);

    // Read until end of the file, put directly into the matrix.
    while (fin) {
        if (has_gene(gene)) {
            // If this phenotype has not already been read, add it to the hash and vector.
            // Check the hash for it since phenotypes might be inserted non-sequentially.
            hash_iterator column_iter = phene_to_column.find(phene);
            if (column_iter == phene_to_column.end()) {
                column_to_phene.push_back(phene);
                phene_to_column[phene] = column_to_phene.size() - 1;
#ifdef DEBUG
                cout << "Added column " << column_to_phene.size()-1 << " phene " << phene << endl;
                assert(column_to_phene[phene_to_column[phene]] == phene);
#endif
            }

            // Now insert the gene/phenotype linkage in the matrix.
            add_association(gene, phene);
        }

        // Next priming read
        read_gene_phene_line(fin, gene, phene);
    }

#ifdef DEBUG
    cout << "Debugging enabled." << endl;
    // Quick check that row_to_gene and gene_to_row are correct
    for (size_t i = 0; i < row_to_gene.size(); ++i) {
        assert(i == gene_to_row[row_to_gene[i]]);
    }
    for (size_t i = 0; i < column_to_phene.size(); ++i) {
        assert(i == phene_to_column[column_to_phene[i]]);
    }
    cout << gene_to_row.size() << endl;
    cout << row_to_gene.size() << endl;
    cout << phene_to_column.size() << endl;
    cout << column_to_phene.size() << endl;
#endif

    // Close the file
    fin.close();
}



genephene::genephene(const genephene& orig)
    : _num_phenes(orig._num_phenes),
      _num_genes(orig._num_genes),
      gene_to_row(orig.gene_to_row),
      phene_to_column(orig.phene_to_column),
      row_to_gene(orig.row_to_gene),
      column_to_phene(orig.column_to_phene),
      gene_by_phene(orig.gene_by_phene),
      phene_distance(NULL),
      _common_items(NULL),
      original_filename(orig.original_filename), // remember to change this if we save a modified matrix.
        species1_(orig.species1_),
        species2_(orig.species2_),
        distance_function(orig.distance_function)
{
    // Make a new copy of the distance matrix if necessary.
    if (orig.phene_distance != NULL)
        phene_distance = new symmetric_matrix<dist_t>(*(orig.phene_distance));

    // Make a new copy of the common_items matrix if necessary
    if (orig._common_items != NULL)
        _common_items = new compressed_matrix<size_t>(*(orig._common_items));
}


genephene::~genephene() {
    if (phene_distance != NULL)        delete phene_distance;
    if (_common_items != NULL)         delete _common_items;
}


bool genephene::save_common_items(string suffix) const {
    
    if (_common_items != NULL) {
        string filename = save_common_items_filename(suffix);
        string pfilename = save_pcommon_items_filename(suffix);

        if (!exists(filename)) {
            cout << "Saving common items matrix." << endl;

            // Save only the upper half (since it's symmetric) and don't bother
            // storing 0s.
            save_upper_matrix<compressed_matrix<size_t>, size_t>(filename, *_common_items, 0, *this, false);
            save_upper_matrix<compressed_matrix<size_t>, size_t>(pfilename, *_common_items, 0, *this, true);
            
            cout << "Done saving common items." << endl;

            return true;

        } else {
            cout << "Not saving common items matrix. File already exists." << endl;
        }
    } else {
        cerr << "Warning: Save requested for _common_items, but matrix is not allocated." << endl;
    }
    return false;
}


bool genephene::save_distances(string suffix, string distance_fun) const {
    if (phene_distance != NULL) {
        string filename              = save_distances_filename(suffix, distance_fun);
        string pfilename             = save_pdistances_filename(suffix, distance_fun);
        //string col_to_phene_filename = save_col_to_phene_filename(suffix);

        //cout << "Saving (if absent) column-to-phene vector. File: " << col_to_phene_filename << endl;
        // Can save a little time by touching this filename:
        // save_vector_if_absent<phene_id_t>(col_to_phene_filename, column_to_phene, true);
        // The above line is only needed when we want to manually decode the
        // distance matrices, or if we're using the phenologdb ruby-on-rails
        // front end I wrote. -- JW 9/25/09
        //cout << "Done saving column-to-phene vector." << endl;

        
        if (!exists(filename)) {
            cout << "Saving distance matrix." << endl;

            save_upper_prob_matrix<symmetric_matrix<dist_t>, dist_t>(filename,  *phene_distance, 1.0, *this, false, 1e-7);
            // The second one saves by phenotype instead of by column.
            save_upper_prob_matrix<symmetric_matrix<dist_t>, dist_t>(pfilename, *phene_distance, 1.0, *this, true,  1e-7);
            
            cout << "Done saving distance matrix." << endl;

            return true;
            
        } else {
            cout << "Not saving distances matrix. File already exists." << endl;
        }

        
    } else {
        cerr << "Warning: Save requested for phene_distance, but matrix is not allocated." << endl;
    }
    
    return false;
}


bool genephene::load_common_items(string suffix) const {
    if (_common_items != NULL) {
        string filename = save_common_items_filename(suffix);

        if (!exists(filename)) {
            cout << "No common items matrix file detected." << endl;
        } else {
            cout << "Loading common items matrix file." << endl;
            // No need to call init_load, because this is a compressed matrix --
            // zeros will not be stored.
            load_matrix<compressed_matrix<size_t>, size_t>(filename, *_common_items);
            cout << "Done loading common items matrix." << endl;
            return true;
        }
        
    } else {
        cerr << "Warning: Save requested for _common_items, but matrix is not allocated." << endl;
    }
    return false;
}


bool genephene::load_distances(string suffix, string distance_fun) const {
    if (phene_distance != NULL) {
        string filename = save_distances_filename(suffix, distance_fun);

        if (!exists(filename)) {
            cout << "No distance matrix file detected." << endl;
        } else {
            cout << "Loading distance matrix file." << endl;
            // Load an inverted probability matrix, unloaded values are set to 1.
            init_load_upper_prob_matrix<symmetric_matrix<dist_t>, dist_t>(filename, *phene_distance, 1);
            cout << "Done loading distance matrix." << endl;
            return true;
        }

    } else {
        cerr << "Warning: Save requested for phene_distance, but matrix is not allocated." << endl;
    }
    return false;
}


const dist_t genephene::operator()(const gene_id_t& g, const phene_id_t& p) const {
    if (gene_to_row.find(g) == gene_to_row.end())
        return (dist_t)(false);
    else if (phene_to_column.find(p) == phene_to_column.end())
        return (dist_t)(false);
    else
        return gene_by_phene( (gene_to_row.find(g))->second,
                              (phene_to_column.find(p))->second );
}


const dist_t genephene::distance(phene_id_t p1, phene_id_t p2) const {
    if (phene_distance == NULL) {
        cerr << "Error: Phenotype distance matrix is not allocated." << endl;
        throw;
    }

    return (*phene_distance)((phene_to_column.find(p1))->second, (phene_to_column.find(p2))->second);
}


const size_t genephene::common_items(phene_id_t p1, phene_id_t p2) const {
    if (_common_items == NULL) {
        cerr << "Warning: Common items matrix is not allocated. Must calculate common_items on the fly." << endl;
        cerr << "Also be aware that the test set will not be taken into account, if there is one!" << endl;
        return gene_by_phene.num_common_edges((phene_to_column.find(p1))->second, (phene_to_column.find(p2))->second);
    }

    return (*_common_items)((phene_to_column.find(p1))->second, (phene_to_column.find(p2))->second);
}


// For a row-based cross-validation test set, convert genes to row indeces.
set<size_t> genephene::test_set_genes_to_rows(const set<gene_id_t>& test_set) const {
    set<size_t> row_test_set;
    for (set<gene_id_t>::const_iterator i = test_set.begin(); i != test_set.end(); ++i) {
        unordered_map<gene_id_t, size_t>::const_iterator row_it = gene_to_row.find(*i);
        if (row_it != gene_to_row.end())
            row_test_set.insert(row_it->second);
    }
    return row_test_set;
}


// For a cell-based cross-validation test set, convert genes to row indeces and
// phenes to column indeces.
// Phenes and genes are switched for the returned multimap. This makes computing
// intersections faster.
multimap<size_t,size_t> genephene::test_set_genes_and_phenes_to_cells(const multimap<gene_id_t, phene_id_t>& test_set) const {
    multimap<size_t,size_t> cell_test_set;

    for (multimap<gene_id_t, phene_id_t>::const_iterator i = test_set.begin(); i != test_set.end(); ++i) {
        unordered_map<gene_id_t,  size_t>::const_iterator row_it  =      gene_to_row.find(i->first);
        unordered_map<phene_id_t, size_t>::const_iterator col_it  =  phene_to_column.find(i->second);
        if (row_it != gene_to_row.end() && col_it != phene_to_column.end())
            cell_test_set.insert( pair<size_t,size_t>(col_it->second, row_it->second) );
    }
    return cell_test_set;
}


void genephene::calculate_distances(const set<gene_id_t>& r_test_set, const multimap<gene_id_t, phene_id_t>& c_test_set, string suffix, string distance_fun) {

    cout << "calculate_distances: before conversion -- rows in test set = " << r_test_set.size() << endl;
    // If there's a test set, change it to matrix rows/columns
    set<size_t>             row_test_set  = test_set_genes_to_rows(r_test_set);
    multimap<size_t,size_t> cell_test_set = test_set_genes_and_phenes_to_cells(c_test_set);
    cout << "calculate_distances: after conversion -- rows in test set = " << row_test_set.size() << endl;
    cout << "calculate_distances: cells in test set = " << cell_test_set.size() << endl;

    // Need common_items in order to calculate distances with hypergeometric fn.
    // This also does calculate_total_nonzeros, but without the extra iterations.
    calculate_common_items(row_test_set, cell_test_set, suffix);

    // Ensure that the distance function is working properly.
    test_distance_function();

    if (phene_distance == NULL) {
        phene_distance = new symmetric_matrix<dist_t>(_num_phenes);

        // Load from file rather than calculate if possible.
        if (load_distances(suffix, distance_fun)) {
            cout << "Distances matrix successfully loaded." << endl;
            return;
        }

        cout << "Calculating distances." << endl;
        
        // Calculate the distances.
        for (size_t j = 0; j < _num_phenes; ++j) {
            // Diagonal will always be zero.
            (*phene_distance)(j,j) = 0;
            
            for (size_t k = j+1; k < _num_phenes; ++k) {
                size_t common_items_j_k = (*_common_items)(j,k);
                if (common_items_j_k == 0)
                    (*phene_distance)(j,k) = 1; // Don't want any wishy-washy 9.9999e-1s.
                else
                    (*phene_distance)(j,k) = (*distance_function)((*_common_items)(j,j),  // defective
                                                               (*_common_items)(k,k),  // drawn
                                                               common_items_j_k,       // drawn defective
                                                               _num_genes - row_test_set.size()); // total
            }

            if (j % 200 == 0)
                cout << "j = " << j << " of " << _num_phenes << endl;
        }

    } else {
        cerr << "Warning: Requested genephene::calculate_distances, but the distance matrix is already allocated." << endl;
    }
}


void genephene::calculate_common_items(const set<size_t>& row_test_set, const multimap<size_t,size_t>& cell_test_set, string suffix) {
    cout << "calculate_common_items:" << endl;

    if (_common_items == NULL) {
        _common_items = new compressed_matrix<size_t>(_num_phenes, _num_phenes, 0);

        // Attempt to load if possible, rather than calculate.
        // Matrix is initalized to zeros by default, since it's compressed.
        if (load_common_items(suffix)) {
            cout << "Successfully loaded existing common items matrix." << endl;
            return;
        }
    } else {
        cerr << "Warning: Requested genephene::calculate_common_items, but the _common_items matrix is already allocated." << endl;
        // throw; // Eventually we can fix this, but for now, let's give an error.
    }

    cout << "Calculating common items matrix." << endl;
    // Calculate the number of 1s in common in each column.
    // TODO: Speed this up.
    for (size_t j = 0; j < _num_phenes; ++j) {
        
        for (size_t k = j; k < _num_phenes; ++k) {

            // Calculate using set intersection on the adjacency list object.
            (*_common_items)(j,k) = gene_by_phene.num_common_edges(j,k, row_test_set, cell_test_set);
        }
    }
    cout << "Done calculating common items matrix." << endl;
}


size_t genephene::num_genes_with_phene(phene_id_t p) const {
    if (_common_items == NULL) {
        cerr << "Error: Common items matrix is not allocated." << endl;
        throw;
    }

    return gene_by_phene.num_rows((phene_to_column.find(p))->second);
}


set<gene_id_t> genephene::gene_set() const {
    // Convert the contents of row_to_gene vector into a set.
    set<gene_id_t> ret;

    // Initial insert, so we have an iterator clue for the next one.
    vector<gene_id_t>::const_iterator i = row_to_gene.begin();
    set<gene_id_t>::iterator set_iter   = ret.insert(*i).first;

    // Start loop at next position
    ++i;
    for (; i != row_to_gene.end(); ++i) {
        // Use each insertion as a clue for the next, since list is sorted.
        set_iter = ret.insert(set_iter, *i);
    }

    return ret;
}


set<phene_id_t> genephene::phene_set() const {
    // Convert the contents of column_to_phene vector into a set.
    set<phene_id_t> ret;

    // Initial insert, so we have an iterator clue for the next one.
    vector<phene_id_t>::const_iterator i = column_to_phene.begin();
    set<phene_id_t>::iterator set_iter   = ret.insert(*i).first;

    // Start loop at next position
    ++i;
    for (; i != column_to_phene.end(); ++i) {
        // Use each insertion as a clue for the next, since list is sorted.
        set_iter = ret.insert(set_iter, *i);
    }

    return ret;
}


// Find the closest phenotype to p in the matrix.
// Only returns one -- the one with the most genes in it.
// Also returns the distance; if distance is 1, we might not want to bother with
// this particular phene.
// Requires that p be for the prediction species (e.g., exists in each matrix).
pair<phene_id_t,dist_t> genephene::nearest(phene_id_t p) const {

    // Translate phenotype to column
    size_t col = safely_translate_phenotype_to_column(p);

    dist_t min_value = 1.0;
    size_t min_index = 0;
    size_t min_genes = 0;

#ifdef DEBUG
    check_species();
#endif

    if (species1_ != species2_) {
        // Make sure this phenotype is not a local phenotype.
        if (local_columns.find(col) != local_columns.end()) {
            cerr << "Error: Species' phenes file does not seem to match the associations file." << endl;
            cerr << "species1 = " << species1_ << endl;
            cerr << "species2 = " << species2_ << endl;
            throw;
        }
        
        // Only check local columns; no need to check the entire matrix.
        for (set<size_t>::iterator jt = local_columns.begin(); jt != local_columns.end(); ++jt) {
            size_t j = *jt; // get the actual column
            if (col == j) continue; // skip diagonal

            dist_t cur_value = (*phene_distance)(j,col);
            if (cur_value < min_value) {
                min_value = cur_value;
                min_index = j;
            } else if (are_equal(cur_value, min_value) && min_value < 1) {
                // If this one has the same value but more genes, keep it.
                size_t cur_genes = (*_common_items)(j,j);
                if (cur_genes > min_genes) min_index = j;
            }
        }
    } else {
        // Exactly the same thing, but all phenotypes.
        for (size_t j = 0; j < _num_phenes; ++j) {
            if (col == j) continue;
            
            dist_t cur_value = (*phene_distance)(j,col);
            if (cur_value < min_value) {
                min_value = cur_value;
                min_index = j;
            } else if (are_equal(cur_value, min_value) && min_value < 1) {
                // If this one has the same value but more genes, keep it.
                size_t cur_genes = (*_common_items)(j,j);
                if (cur_genes > min_genes) min_index = j;
            }
        }
    }

    // Return phene id and distance.
    return pair<phene_id_t,dist_t>(column_to_phene[min_index], min_value);
}


// Like nearest but returns them all. Sorting is automatic, as we're using a multiset.
// Multiset instead of set because two distances MAY be the same.
// Will not return anything >= upperbound.
multiset<dist_phene_t>* genephene::sorting(phene_id_t p, dist_t upperbound) const {

    // Values to be returned.
    multiset<dist_phene_t>* ret = new multiset<dist_phene_t>();

    // Translate phenotype to column
    size_t col = safely_translate_phenotype_to_column(p);

    check_species();

    if (species1_ != species2_) {
        // Make sure this phenotype is not a local phenotype.
        if (local_columns.find(col) != local_columns.end()) {
            cerr << "Error: Species' phenes file does not seem to match the associations file." << endl;
            cerr << "species1 = " << species1_ << endl;
            cerr << "species2 = " << species2_ << endl;
            cerr << "In other words, this column is not in local_columns." << endl;
            throw;
        }

        // Only check local columns; no need to check the entire matrix.
        for (set<size_t>::iterator jt = local_columns.begin(); jt != local_columns.end(); ++jt) {
            size_t j = *jt; // get the actual column
            if (col == j) continue; // skip diagonal

            dist_t cur_value = (*phene_distance)(j,col);
            if (cur_value < upperbound) // Don't want to return things that are 'infinite' (1.0) distance away.
                ret->insert(dist_phene_t(cur_value, column_to_phene[j]));
        }
    } else {
        // Exactly the same thing, but all phenotypes.
        for (size_t j = 0; j < _num_phenes; ++j) {
            if (col == j) continue;

            dist_t cur_value = (*phene_distance)(j,col);
            if (cur_value < upperbound) // Don't want to return things that are 'infinite' (1.0) distance away.
                ret->insert(dist_phene_t(cur_value, column_to_phene[j]));
        }
    }

    return ret;
}

size_t genephene::safely_translate_phenotype_to_column(const phene_id_t& p) const {
    unordered_map<phene_id_t,size_t>::const_iterator find_p = phene_to_column.find(p);
    if (find_p == phene_to_column.end()) {
        cerr << "Error: Attempted to get nearest for a phenotype not in this matrix." << endl;
        throw;
    }
    return find_p->second;
}


void genephene::add_to_sorting(dist_phene_multiset_t& to, phene_id_t p, size_t k, dist_t upperbound) const {
    // Translate phenotype to column
    size_t col = safely_translate_phenotype_to_column(p);

    check_species();

    if (species1_ != species2_) {
        // Make sure this phenotype is not a local phenotype.
        if (local_columns.find(col) != local_columns.end()) {
            cerr << "Error: Species' phenes file does not seem to match the associations file." << endl;
            cerr << "species1 = " << species1_ << endl;
            cerr << "species2 = " << species2_ << endl;
            cerr << "In other words, this column is not in local_columns." << endl;
            throw;
        }

        // Only check local columns; no need to check the entire matrix.
        for (set<size_t>::iterator jt = local_columns.begin(); jt != local_columns.end(); ++jt) {
            size_t j = *jt; // get the actual column
            if (col == j) continue; // skip diagonal

            // Don't want to return things that are 'infinite' (1.0) distance away.
            dist_t cur_value = (*phene_distance)(j,col);
            if (cur_value < upperbound)
                to.insert(dist_phene_t(cur_value, column_to_phene[j]));

        }
    } else {
        // Exactly the same thing, but all phenotypes.
        for (size_t j = 0; j < column_to_phene.size(); ++j) {
            if (col == j) continue;

            // Don't want to return things that are 'infinite' (1.0) distance away.
            dist_t cur_value = (*phene_distance)(j,col);
            if (cur_value < upperbound)
                to.insert(dist_phene_t(cur_value, column_to_phene[j]));
        }
    }

    if (to.size() > k && k > 0) {
        // Delete all but the first k items (or k+m where the m items after k
        // are equal to item k's value).
        size_t kpos                             = 0;
        double prev_item                        = 0;
        dist_phene_multiset_t::iterator kt      = to.begin();
        dist_phene_multiset_t::iterator prev_kt = to.begin();
        
        while (kpos < k || are_equal(kt->first(), prev_item)) {
            prev_item = kt->first();
            prev_kt = kt;
            
            ++kpos;
            ++kt;
        }
        // Found the last position. Perform the actual deletion.
        for (size_t i = 0; i < to.size() - kpos && kt != to.end(); ++i) {
            to.erase(kt);
            kt = prev_kt;
            ++kt;
        }
    }
    // cout << "to.size() = " << to.size() << endl;
/*    if (to.size() > 200) {
        for (dist_phene_multiset_t::iterator it = to.begin(); it != to.end(); ++it) {
            cout << it->first() << ' ' << it->second() << ", ";
        }
        cout << endl;
        throw;
    }
*/
}

void genephene::check_species() const {
    if (species1_.size() < 2 || species2_.size() < 2) {
        cerr << "Error: species pair not set on phenomatrix." << endl;
        cerr << "species1_ = '" << species1_ << "'; species2_ = '" << species2_ << "'" << endl;
        throw;
    }
}

set<phene_id_t> genephene::local_phenes() const {
    set<phene_id_t> ret;

    // If source and destination species is the same, all columns are local.
    if (species1_ == species2_) {
        for (vector<phene_id_t>::const_iterator i = column_to_phene.begin(); i != column_to_phene.end(); ++i)
            ret.insert(*i);
    } else {
        for (set<size_t>::const_iterator i = local_columns.begin(); i != local_columns.end(); ++i)
            ret.insert(column_to_phene[*i]);
    }
    return ret;
}


void genephene::test_distance_function() const {
    using std::cout;
    using std::flush;
    using std::endl;

    cout << "Testing distance function:" << endl;
    cout << " 1. Testing for zero result..." << flush;
    if (are_equal<double>( (*distance_function)(10, 10, 10, 100), 0.0 )) cout << "OK" << endl;
    else {
        cout << "FAIL: " << (*distance_function)(10, 10, 0, 100) << endl;
        cerr << "Distance function must yield 0 for columns that have no items in common." << endl;
        throw;
    }

    cout << " 1. Testing for non-zero result..." << flush;
    if (!are_equal<double>( (*distance_function)(10, 10, 0, 100), 0.0 )) cout << "OK" << endl;
    else {
        cout << "FAIL: " << (*distance_function)(10, 10, 0, 100) << endl;
        cerr << "Warning: Distance function for cross-validation and prediction must yield exactly 1.0 for columns that are the same, and the current function does not. This is okay if you're calculating a distribution." << endl;
        // Don't throw, since we're probably calculating a distribution if this test fails.
    }
}