#include "model.h"

/// Apply distortion model to every column of xy (must have two rows)
void Model::applyMat(libNumerics::matrix<flnum>& xy) const {
    assert(xy.nrow()==2);
    for(int i=0; i<xy.ncol(); i++)
        apply(xy(0,i), xy(1,i));
}
