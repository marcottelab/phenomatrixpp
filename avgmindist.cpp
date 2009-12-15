/* 
 * File:   avgmindist.cpp
 * Author: jwoods
 * 
 * Created on July 22, 2009, 4:13 PM
 */

#include "avgmindist.h"

avgmindist::avgmindist(list<marshall>& mar, const marshall& predmar, const set<gene_id_t>& row_test_set, const multimap<gene_id_t,phene_id_t>& cell_test_set, const string& identifier, const string& distance_function)

        // true specifies that we should use weighted predictions. Takes more memory,
        // useless when only doing nearest, but necessary for all other methods.
: oracle(mar, predmar, row_test_set, cell_test_set, identifier, distance_function, true)
{
    typedef std::pair<phene_id_t,double> pd;

    set<phene_id_t> to_predict = dest->local_phenes();

    // Find each species' closest and stick them all in the vector.
    for (set<phene_id_t>::const_iterator i = to_predict.begin(); i != to_predict.end(); ++i) {
        // Initialize the map so it behaves like a vector and we don't have to
        // use find all the time.
        nearest[*i] = list<phene_id_t>();

        for (sp_genephene_map_t::iterator j = genephenes.begin(); j != genephenes.end(); ++j) {
            if (j->second->has_phene(*i)) {

                // Get the nearest and add it to the map.
                pd cur(j->second->nearest(*i));
                nearest[*i].push_back(cur.first);
            }
        }
    }
}
/*
avgmindist::avgmindist(const avgmindist& orig) : oracle({
    cerr << "need to test copy constructors for oracle derivatives." << endl;
    throw;
} */

avgmindist::~avgmindist() {
}


// Predict for all phenotypes in the dest matrix.
void avgmindist::predict() {

    // Get the two items over which we will iterate
    const vector<phene_id_t>& phenes = dest->phene_vector();
    const vector<gene_id_t>&  genes  = dest->gene_vector();

    cout << "Getting total number of species with each gene." << endl;
    // Get the total number of species with each gene.
    unordered_map<gene_id_t,size_t> num_species_with_gene(genes.size());
    for (vector<gene_id_t>::const_iterator n = genes.begin(); n != genes.end(); ++n)
        num_species_with_gene[*n] = get_num_species_with_gene(*n);

    cout << "Predicting for each phenotype." << endl;
    // Predict for each phenotype
    for (vector<phene_id_t>::const_iterator j = phenes.begin(); j != phenes.end(); ++j) {

        // Only bother if there's a nearest.
        nearest_t::const_iterator closest = nearest.find(*j);

        if (closest != nearest.end()) {

            // For each of the n nearest phenotypes, get the information from
            // the genes in that phenotype's species matrix.
            for (list<phene_id_t>::const_iterator c = closest->second.begin(); c != closest->second.end(); ++c) {
                unordered_map<phene_id_t,species_id_t>::const_iterator find_p = phene_sp.find(*c);
                if (find_p == phene_sp.end()) {
                    cerr << "Error: Could not find phenotype " << *c << " from phene_sp size = " << phene_sp.size() << endl;
                    cerr << "Error: That ID was produced with nearest.find(*j), with *j = " << *j << " and nearest size = " << nearest.size() << endl;
                    throw;
                }
                // Using the phenotype, figure out which species.
                species_id_t from_species = find_p->second;

                for (vector<gene_id_t>::const_iterator i = genes.begin(); i != genes.end(); ++i) {

                    
                    // Assign the new value for each gene
                    // No need to check has_gene because this particular entry would
                    // only be marked as nearest if has_gene is true. // removed: dest->species1() != from_species &&
                    if ((genephenes)[from_species]->has_association(*i,*c))
                        dest->add_to_association(*i, *j);
                    // Divisor causes us to take the average.
                }
            }

            // Now take the average.
            for (vector<gene_id_t>::const_iterator i = genes.begin(); i != genes.end(); ++i) {
                if (num_species_with_gene[*i] > 0)
                    dest->divide_association_by(*i, *j, num_species_with_gene[*i]);
            }
        }
    }
}
