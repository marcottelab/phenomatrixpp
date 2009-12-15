/* 
 * File:   mindist.cpp
 * Author: jwoods
 * 
 * Created on July 16, 2009, 2:32 PM
 */

#include "mindist.h"
mindist::mindist(list<marshall>& mar, const marshall& predmar, const set<gene_id_t>& row_test_set, const multimap<gene_id_t,phene_id_t>& cell_test_set, const string& identifier, const string& distance_function)
        : oracle(mar, predmar, row_test_set, cell_test_set, identifier, distance_function, false)
{

    typedef std::pair<phene_id_t,double> pd;
    // Initialize oracle in constructor list.
    // Now that that's done, we can calculate the minimums and stick them in the
    // vector.

    set<phene_id_t> to_predict = dest->local_phenes();

    for (set<phene_id_t>::const_iterator i = to_predict.begin(); i != to_predict.end(); ++i) {
        // Find the minimum for each of these phenes
        pd min(0, 1);

        for (sp_genephene_map_t::iterator j = genephenes.begin(); j != genephenes.end(); ++j) {
            if (j->second->has_phene(*i)) {
                // Which is smallest?
                pd cur(j->second->nearest(*i));
                if (cur.second < min.second)
                    min = cur;
            }
        }

        if (min.first > 0)
            // Save it.
            nearest[*i] = min.first;
    }
}

mindist::mindist(const mindist& orig)
: oracle(orig), nearest(orig.nearest)
{
    cerr << "need to test copy constructors for oracle derivatives." << endl;
    throw;
}

mindist::~mindist() {
}


// Predict for all phenotypes in the dest matrix.
void mindist::predict() {
    // Get the two items over which we will iterate
    const vector<phene_id_t>& phenes = dest->phene_vector();
    const vector<gene_id_t>&  genes  = dest->gene_vector();

    // Predict for each phenotype
    for (vector<phene_id_t>::const_iterator j = phenes.begin(); j != phenes.end(); ++j) {

        // Only bother if there's a nearest.
        map<phene_id_t,phene_id_t>::const_iterator closest = nearest.find(*j);
        
        if (closest != nearest.end()) {
            // Get the phenotype identifier
            phene_id_t from_phenotype = closest->second;
            unordered_map<phene_id_t,species_id_t>::const_iterator find_p = phene_sp.find(from_phenotype);
            if (find_p == phene_sp.end()) {
                cerr << "Error: Could not find phenotype " << from_phenotype << " from phene_sp size = " << phene_sp.size() << endl;
                cerr << "Error: That ID was produced with nearest.find(*j), with *j = " << *j << " and nearest size = " << nearest.size() << endl;
                throw;
            }
            // Using the phenotype, figure out which species.
            species_id_t from_species = find_p->second;

            for (vector<gene_id_t>::const_iterator i = genes.begin(); i != genes.end(); ++i) {

                // Assign the new value for each gene
                // No need to check has_gene because this particular entry would
                // only be marked as nearest if has_gene is true.
                if (genephenes[from_species]->has_association( *i, from_phenotype ))
                    dest->add_association(*i, *j);
                // (*dest)( *i, *j ) = (*(genephenes[from_species]))( *i, from_phenotype );
            }
        }
    }
}
