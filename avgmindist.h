/* 
 * File:   avgmindist.h
 * Author: jwoods
 *
 * Created on July 22, 2009, 4:13 PM
 */

#ifndef _AVGMINDIST_H
#define	_AVGMINDIST_H

#include "oracle.h"

class avgmindist : public oracle {
    typedef map<phene_id_t,list<phene_id_t> > nearest_t;
public:
    avgmindist(list<marshall>& mar, const marshall& predmar, const set<gene_id_t>& row_test_set, const multimap<gene_id_t,phene_id_t>& cell_test_set, const string& identifier, const string& distance_function);
    //avgmindist(const avgmindist& orig);
    virtual ~avgmindist();

    virtual string identify_method() const {
        return "nearest of each species (averaged)";
    }
    virtual string identify_column() const {
        return "avgmindist";
    }

    virtual void predict();
protected:
    // Keep track of the nearest items. Key is a dest-org phenotype, probably
    // human (the one to be predicted).
    nearest_t nearest;
};

#endif	/* _AVGMINDIST_H */

