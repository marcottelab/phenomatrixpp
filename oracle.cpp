/* 
 * File:   oracle.cpp
 * Author: jwoods
 * 
 * Created on July 16, 2009, 1:57 PM
 */

#include "oracle.h"

void create_directory_if_necessary(const string& dirname) {
    if (!exists(dirname)) {
        cout << "Creating directory '" << dirname << "'." << endl;
        create_directory(dirname);
    } else if (!is_directory(dirname)) {
        cerr << "Error: Cannot create directory; name " << dirname << " is taken!";
        throw;
    }
}


progress_display* save_progress_indicator(const vector<phene_id_t>& vec) {
    // For each phenotype, get the matrix column we want to print, and also
    // what is known.
    cout << "Total items to save: " << vec.size() << endl << endl;
    return new progress_display(vec.size());
}


oracle::oracle(
        list<marshall>& mar,
        const marshall& predmar,
        const set<gene_id_t>& row_test_set_,
        const multimap<gene_id_t, phene_id_t>& cell_test_set_,
        const string& identifier,
        const string& distance_function,
        bool weighted_predictions)
        :
        dest( predmar() ),
        row_test_set( row_test_set_ ),
        cell_test_set( cell_test_set_ )
{

    // Get the constructed matrices.
    for (list<marshall>::iterator i = mar.begin(); i != mar.end(); ++i) {
        genephenes.insert(make_pair(i->species(), (*i)()));

        // Get the phenotypes that map to this species
        set<phene_id_t> phenes = genephenes[i->species()]->local_phenes();
        for (set<phene_id_t>::const_iterator j = phenes.begin(); j != phenes.end(); ++j)
            phene_sp[*j] = i->species();

        genephenes[i->species()]->calculate_distances(row_test_set, cell_test_set, identifier, distance_function);
        genephenes[i->species()]->save_distances(identifier, distance_function);
        genephenes[i->species()]->save_common_items(identifier);
    }
    mar.clear();
}


oracle::oracle(const oracle& orig)
: genephenes(orig.genephenes), phene_sp(orig.phene_sp)
{
    cerr << "Copy constructor not written." << endl;
    throw;
}

oracle::~oracle() { }


gdistmap oracle::find_predictions_by_phenotype(phene_id_t p, const gvector& genes) const {

    // Variable to return.
    gdistmap predicted(genes.size());

    // Iterate down the rows (we're within a column iteration)
    for (gvector::const_iterator i = genes.begin(); i != genes.end(); ++i)
        predicted[*i] = (*dest)(*i,p);

    return predicted;
}

gdistmap oracle::find_row_test_predictions_by_phenotype(phene_id_t p, size_t genes_size) const {

    // Variable to return.
    gdistmap predicted(genes_size);

    // Iterate down the rows (we're within a column iteration)
    for (gset::const_iterator i = row_test_set.begin(); i != row_test_set.end(); ++i)
        predicted[*i] = (*dest)(*i,p);

    return predicted;
}

// As arguments, takes the cell test set (reversed, so it's phenes - genes instead of genes - phenes)
// and the specific phenotype ID.
gdistmap oracle::find_cell_test_predictions_by_phenotype(const pgtable& rev_cell_test_set, phene_id_t p) const {

    // Variable to return.
    gdistmap predicted(rev_cell_test_set.size());

    // Iterate through the cells in the test set.
    for (pgtable::const_iterator i = rev_cell_test_set.lower_bound(p); i != rev_cell_test_set.upper_bound(p); ++i)
        predicted[i->second] = (*dest)(i->second, p);

    return predicted;
}


// Reverse the cell test set for writing predictions.
// TODO: Speedup.
pgtable oracle::reverse_of_cell_test_set() const {
    pgtable res;
    for (gptable::const_iterator i = cell_test_set.begin(); i != cell_test_set.end(); ++i)
        res.insert(make_pair<phene_id_t,gene_id_t>(i->second, i->first));
    return res;
}


void oracle::write_predictions_by_phenotype(const string& dirname, phene_id_t p, const gdistmap& predicted) const {
    if (!create_prediction_file(dirname, p, identify_method(), identify_column(), predicted))
        add_column_to_prediction_file(dirname, p, identify_method(), identify_column(), predicted);
}


void oracle::save_predictions(const string& dirname) const {

    cout << "Saving all predictions." << endl;

    // Get phenotypes and genes to predict for
    const gvector& genes  = dest->gene_vector(); assert(genes.size() > 0);
    const pvector& phenes = dest->phene_vector();

    create_directory_if_necessary(dirname);

    // For each phenotype, get the matrix column we want to print, and also
    // what is known.
    cout << "Total items to save: " << phenes.size() << endl << endl;
    progress_display show_progress(phenes.size());

    for (pvector::const_iterator j = phenes.begin(); j != phenes.end(); ++j, ++show_progress)
        // Map genes to the prediction value and write them
        write_predictions_by_phenotype(dirname, *j, find_predictions_by_phenotype(*j, genes));
}

void oracle::save_row_test_predictions(const string& dirname) const {

    cout << "Saving row test set predictions." << endl;

    // Get phenotypes and genes to predict for
    const gvector& genes  = dest->gene_vector(); assert(genes.size() > 0);
    const pvector& phenes = dest->phene_vector();

    create_directory_if_necessary(dirname);

    // For each phenotype, get the matrix column we want to print, and also
    // what is known.
    cout << "Total items to save: " << phenes.size() << endl << endl;
    progress_display show_progress(phenes.size());

    for (pvector::const_iterator j = phenes.begin(); j != phenes.end(); ++j, ++show_progress)
        // Map genes to the prediction value and write them
        write_predictions_by_phenotype(dirname, *j, find_row_test_predictions_by_phenotype(*j, genes.size()));

}

void oracle::save_cell_test_predictions(const string& dirname) const {

    cout << "Saving cell test set predictions." << endl;

    // Get phenotypes and genes to predict for
    const gvector& genes  = dest->gene_vector(); assert(genes.size() > 0);
    const pvector& phenes = dest->phene_vector();

    create_directory_if_necessary(dirname);

    // For each phenotype, get the matrix column we want to print, and also
    // what is known.
    cout << "Total items to save: " << phenes.size() << endl << endl;
    progress_display show_progress(phenes.size());

    pgtable rev = reverse_of_cell_test_set();

    for (pvector::const_iterator j = phenes.begin(); j != phenes.end(); ++j, ++show_progress)
        // Map genes to the prediction value and write them
        write_predictions_by_phenotype(dirname, *j, find_cell_test_predictions_by_phenotype(rev, *j) );
}



size_t oracle::get_num_species_with_gene(const gene_id_t& g) const {
    size_t num_species_with_gene = 0;
    // Check each species.
    for (sp_genephene_map_t::const_iterator sp = genephenes.begin(); sp != genephenes.end(); ++sp) {
        if (sp->second->has_gene(g))
            num_species_with_gene++;
    }
    return num_species_with_gene;
}

// Return all of the phenotypes with < upperbound distance from p, from each of
// the phenomatrices.
multiset<dist_phene_t>* oracle::sorting(phene_id_t p, dist_t upperbound) const {
    // Combine all species' sortings into one single one and return it.
    multiset<dist_phene_t>* combined = new multiset<dist_phene_t>();

    // Get the set from each species' phenomatrix.
    for (phene_sp_map_t::const_iterator it = phene_sp.begin(); it != phene_sp.end(); ++it) {
        sp_genephene_map_t::const_iterator spit = genephenes.find(it->second);
        if (!spit->second->has_phene(p))
            continue; // Can't get a sorting if this phene isn't in the matrix.
        multiset<dist_phene_t>* this_species = spit->second->sorting(p, upperbound);

        combined->insert(this_species->begin(), this_species->end());

        // Clean up the pointer
        delete this_species;
    }

    return combined;
}

void oracle::add_to_sorting(dist_phene_multiset_t& to, phene_id_t p, size_t k, dist_t upperbound) {
    cout << "Called add_to_sorting on oracle with p = " << p << endl;
    // Get the set from each species' phenomatrix.

    for (sp_genephene_map_t::const_iterator gpit = genephenes.begin(); gpit != genephenes.end(); ++gpit) {

        // cout << "add_to_sorting: Looking at species " << it->second << endl;
        if (gpit->second->has_phene(p))
            gpit->second->add_to_sorting(to, p, k, upperbound);
    //    else
    //        cout << "Species does not have that phenotype." << endl;
    }
}

