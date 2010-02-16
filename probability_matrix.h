/* 
 * File:   probability_matrix.h
 * Author: jwoods
 *
 * Created on July 22, 2009, 1:28 PM
 */

#ifndef _PROBABILITY_MATRIX_H
#define	_PROBABILITY_MATRIX_H

class probability_matrix {
public:
    probability_matrix();
    probability_matrix(probability_matrix& orig);
    virtual ~probability_matrix();
protected:
    bool create_mode;
};

#endif	/* _PROBABILITY_MATRIX_H */

