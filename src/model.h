#ifndef MODEL_H
#define MODEL_H

#include "libNumerics/numerics.h"

/// Base class for distortion model
class Model {
public:
    typedef libNumerics::flnum flnum;
public:
    virtual ~Model() {}
    virtual void apply(flnum& x, flnum& y) const =0;
    virtual void applyMat(libNumerics::matrix<flnum>& xy) const;
    virtual void evaluate(libNumerics::matrix<flnum>& /*ptsU*/,
                          libNumerics::matrix<flnum>& /*ptsD*/) {}
    virtual void change_order(int /*order*/) {}
};

#endif
