/* 
 * File:   partialbayes.h
 * Author: jwoods
 *
 * Created on August 1, 2009, 7:35 PM
 */

#ifndef _PARTIALBAYES_H
#define	_PARTIALBAYES_H

#include "knearest.h"

class partialbayes : public knearest {
public:
    partialbayes(list<marshall>& mar, marshall& predmar, size_t kval, const set<gene_id_t>& test_set, const multimap<gene_id_t, phene_id_t>& cell_test_set, const string& identifier, const string& distance_function);
    virtual ~partialbayes();

    virtual void predict();
private:

};

#endif	/* _PARTIALBAYES_H */

