/* 
 * File:   marshall.cpp
 * Author: jwoods
 * 
 * Created on July 28, 2009, 4:24 PM
 */

#include "marshall.h"

// Actual constructor for Hs x nonHs and Hs x Hs matrices (source, not prediction).
// Takes two association files. The first argument is a file containing a list of
// orthologs between the two species.
// The first of each pair is the species identifier (Hs, Mm, etc).
void marshall::construct() {

    // Check correctness.
    size_t throwcount = 0;
    if (!exists(orthologs_filename.c_str())) {
        cerr << "Error: File with list of orthologs for source/dest combination " << source_species_info.first << "/" << dest_species_info.first << " not found." << endl;
        cerr << "Filename should be " << orthologs_filename << endl;
        throwcount++;
    }
    if (!exists(source_species_info.second)) {
        cerr << "Error: Source species associations file '" << source_species_info.second << "' not found." << endl;
        throwcount++;
    }
    if (dest_species_info.second != "" && !exists(dest_species_info.second)) {
        cerr << "Error: Destination species association file '" << dest_species_info.second << "' not found." << endl;
        throwcount++;
    }
    if (throwcount > 0) throw;

    source_species_ = source_species_info.first;

    cout << "Opening orthologs file: " << orthologs_filename << endl;

    // Read the set of orthologs
    ifstream oin(orthologs_filename.c_str());
    ignore_header(oin);

    set<gene_id_t> Orthologs = read_set<gene_id_t>(oin);

    oin.close();
    

    cout << "Opening source gene-phenotype file: " << source_species_info.second << endl;

    // First read the organism which is the primary source of information.
    // This sets how we will read the destination organism.
    ifstream sin(source_species_info.second.c_str());
    ignore_header(sin);

    set<phene_id_t> local_phenes;
    // Get all of them
    multimap<gene_id_t,phene_id_t> Smap(read_associations<gene_id_t,phene_id_t>(sin));

    if (dest_species_info.first != source_species_info.first) {
        // Keep track of phenes which are the source species.
        for (multimap<gene_id_t,phene_id_t>::iterator i = Smap.begin(); i != Smap.end(); ++i)
            local_phenes.insert(i->second);
    }

    if (Smap.size() == 0) {
        cerr << "Error: Source species had 0 associations." << endl;
        throw;
    } else cout << "Read " << Smap.size() << " associations." << endl;

    map<phene_id_t, size_t> Dphene_genes; // keep track of number of genes per phene.
    multimap<phene_id_t,gene_id_t> rDmap; // reversed from Smap so we can count entries.

    sin.close();

    if (dest_species_info.second != "") {
        cout << "Opening destination gene-phenotype file: " << dest_species_info.second << endl;

        ifstream din(dest_species_info.second.c_str());
        ignore_header(din);

        gene_id_t dgene;
        phene_id_t dphene;
        din >> dgene;
        din >> dphene;
        while (din) {
            rDmap.insert(make_pair(dphene, dgene));

            // Read the next line
            din >> dgene;
            din >> dphene;
        }


        // Take elements from rDmap and insert them in Smap.
        for (multimap<phene_id_t,gene_id_t>::const_iterator i = rDmap.begin(); i != rDmap.end(); ++i) {
            if (rDmap.count(i->first) >= min_genes && Smap.find(dgene) != Smap.end()) {
                // Only allow this entry to be included in the matrix if:
                // * there are at least min_genes in the dest phenotype
                // * the gene exists in the source species already.
                Smap.insert(make_pair(i->second, i->first));
            }
        }
    }

    // Get the unique values in the maps.
    set<phene_id_t> phenes;

    for (multimap<gene_id_t,phene_id_t>::iterator i = Smap.begin(); i != Smap.end(); ++i)
        phenes.insert(i->second);

    // Figure out which distance function to use.
    double (*distfn)(size_t,size_t,size_t,size_t) = switch_distance_function(distance_measure);

    cout << "Creating matrix for species pair: " << dest_species_info.first << ", " << source_species_info.first << endl;
    // Create matrix.
    // Set species information.
    if (dest_species_info.first == source_species_info.first)
        Matrix = shared_ptr<genephene>(new genephene(distfn, Orthologs, phenes.size(), Smap, matrix_identifier, set<phene_id_t>()));
    else
        Matrix = shared_ptr<genephene>(new genephene(distfn, Orthologs, phenes.size(), Smap, matrix_identifier, local_phenes));
    Matrix->species1(dest_species_info.first);
    Matrix->species2(source_species_info.first);
    source_species_ = source_species_info.first;
}


// Constructor for prediction matrix
void marshall::construct(bool weighted = false) {

    // Check correctness.
    size_t throwcount = 0;
    if (!exists(orthologs_filename)) {
        cerr << "Error: File with list of orthologs for " << dest_species_info.first << ", '" << orthologs_filename << "' << not found." << endl;
        throwcount++;
    }
    if (!exists(phenotypes_filename)) {
        cerr << "Error: File with list of phenotypes for " << dest_species_info.second << ", '" << phenotypes_filename << "' << not found." << endl;
        throwcount++;
    }

    if (throwcount > 0) throw;

    // Read the set of orthologs
    ifstream oin(orthologs_filename.c_str());
    ignore_header(oin);

    set<gene_id_t> Orthologs = read_set<gene_id_t>(oin);

    oin.close();

    ifstream pin(phenotypes_filename.c_str());
    ignore_header(pin);

    set<phene_id_t> Phenotypes = read_set<phene_id_t>(pin);

    pin.close();

    // Figure out which distance function to use.
    double (*distfn)(size_t,size_t,size_t,size_t) = switch_distance_function(distance_measure);

    // Create matrix.
    // Set species information.
    Matrix = shared_ptr<genephene>(
       new genephene(distfn, Orthologs, Phenotypes, weighted)
            );
    Matrix->species1(dest_species_info.first);
    Matrix->species2(dest_species_info.first);
    source_species_ = dest_species_info.first;
}

/*
marshall::marshall(const marshall& orig) {
    throw;
} */

marshall::~marshall() {
}

void marshall::ignore_header(istream& in) const {
    // Figure out what line we're on.
    char c = in.get();
    char d = in.get();
    while ((c >= 'A' && c <= 'Z' && d >= 'a' && d <= 'z') || (c >= 'a' && d <= 'z')) {
        // Put them back so we can say what we ignored.
        in.putback(d);
        in.putback(c);

        // Ignore the header
        string header;
        getline(in, header);

        cerr << "Ignored header line '" << header << "'" << endl;

        if (has_numeric_component(header)) {
            cerr << "Error: Header line contains numeric component. This is like a gene or phenotype label, and ignore_header needs to be modified." << endl;
            throw;
        }

        // Get the next two characters.
        c = in.get();
        d = in.get();
    }
    // When the loop fails put them back.
    in.putback(d);
    in.putback(c);

    // Now we're done.
}
