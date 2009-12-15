/* 
 * File:   mindist.h
 * Author: jwoods
 *
 * Created on July 16, 2009, 2:32 PM
 */

#ifndef _MINDIST_H
#define	_MINDIST_H

#include "oracle.h"

class mindist : public oracle {
public:
    mindist(list<marshall>& mar, const marshall& predmar, const set<gene_id_t>& row_test_set, const multimap<gene_id_t,phene_id_t>& cell_test_set, const string& identifier, const string& distance_function);
    mindist(const mindist& orig);
    virtual ~mindist();

    virtual string identify_method() const {
        return "nearest";
    }
    virtual string identify_column() const {
        return "mindist";
    }

    virtual void predict();
protected:
    // Keep track of the nearest. Key is a dest-org phenotype,
    // probably human (the one to be predicted).
    map<phene_id_t,phene_id_t> nearest;
};

#endif	/* _MINDIST_H */

