/* 
 * File:   marshall.h
 * Author: jwoods
 *
 * Created on July 28, 2009, 4:24 PM
 */

#ifndef _MARSHALL_H
#define	_MARSHALL_H

#include "genephene.h"
#include "utilities.h"
#include "distance.h"

class marshall {
    typedef pair<string,string> str_pair_t;
public:
    marshall() { }
    // Constructor for prediction matrix.
    marshall(const string& distance_measure, const string& species, const string& genefn, const string& phenefn, bool weighted=false);

    virtual ~marshall();
    string species() const { return source_species_; }

    // Get the created matrix.
    shared_ptr<genephene> operator()() const {
        return Matrix;
    }

    // Make these public to avoid obnoxious constructors.
    string     distance_measure;
    string     matrix_identifier;
    string     orthologs_filename;
    string     phenotypes_filename;
    str_pair_t source_species_info;
    str_pair_t dest_species_info;

    // Build it and give the user the pointer.
    void construct();       // build a source matrix
    void construct(bool);   // build a prediction matrix

protected:
    void ignore_header(istream& in) const;

    // Create the matrix and store it here until user asks for the pointer.
    shared_ptr<genephene> Matrix;
    string source_species_;
};

#endif	/* _MARSHALL_H */
