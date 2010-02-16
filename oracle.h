/* 
 * File:   oracle.h
 * Author: jwoods
 *
 * Created on July 16, 2009, 1:57 PM
 */

#ifndef _ORACLE_H
#define	_ORACLE_H

#include <utility>
#include <list>
using std::list;
using std::make_pair;

// #include "genephene.h"
#include "marshall.h"
#include "utilities.h"

using boost::filesystem::create_directory;

typedef vector<phene_id_t> pvector;
typedef vector<gene_id_t> gvector;
typedef set<gene_id_t> gset;
typedef multimap<gene_id_t, phene_id_t> gptable;
typedef multimap<phene_id_t, gene_id_t> pgtable;
typedef unordered_map<gene_id_t,dist_t> gdistmap;

void create_directory_if_necessary(const string&);
progress_display* save_progress_indicator(const pvector&);

// Pure virtual base class for predictions. Derived directly from this are
// polyoracle and mindist.
// This class contains much of the interface for the prediction stuff. The data
// structures are specializations.
// It does contain the data structures representing the different genephene
// matrices.
class oracle {
public:
    // Provide a list of filenames containing gene-phenotype associations. These
    // will be passed to genephene.
    oracle(list<marshall>& mar, const marshall& predmar, const gset& row_test_set_, const gptable& cell_test_set_, const string& identifier, const string& distance_function, bool weighted_predictions = false);
    oracle(const oracle& orig);
    virtual ~oracle();

    virtual string identify_column() const = 0;
    virtual string identify_method() const = 0;

    virtual void save_predictions(const string&) const;
    virtual void save_row_test_predictions (const string&) const;
    virtual void save_cell_test_predictions(const string&) const;

    // pure virtual function depends entirely on the
    // child class.
    virtual void predict() = 0;


    size_t get_num_species_with_gene(const gene_id_t& g) const;

    // Return all of the phenotypes with < upperbound distance from p.
    multiset<dist_phene_t>* sorting(phene_id_t p, dist_t upperbound = 1) const;
    void add_to_sorting(dist_phene_multiset_t& to, phene_id_t p, size_t k, dist_t upperbound);

    pgtable reverse_of_cell_test_set() const;

protected:
    void write_predictions_by_phenotype(const string&, phene_id_t, const gdistmap&) const;
    gdistmap find_cell_test_predictions_by_phenotype(const pgtable&, phene_id_t) const;
    gdistmap find_row_test_predictions_by_phenotype(phene_id_t, size_t) const;
    gdistmap find_predictions_by_phenotype(phene_id_t, const gvector&) const;
    
    // Hash by species identifier to genephene matrix (pointer).
    sp_genephene_map_t genephenes;

    // Hash by phenotype identifier to species identifier.
    phene_sp_map_t phene_sp;

    // dest stands for destination: phenotypes for which we want to
    // make predictions.
    shared_ptr<genephene> dest;

    // For cross-validation. If empty, test everything.
    gset    row_test_set; // row test set
    gptable cell_test_set;
};

#endif	/* _ORACLE_H */

