#ifndef MODELPOLY_H
#define MODELPOLY_H

#include "model.h"

/// Bivariate polynomial model
class ModelPoly : public Model {
    int order_;
    std::vector<flnum> a_, b_;
    void powers(flnum x, flnum* xi) const;
public:
    explicit ModelPoly(int order);
    void change_order(int order);

    void apply(flnum& x, flnum& y) const;
    void evaluate(libNumerics::matrix<flnum>& ptsU,
                  libNumerics::matrix<flnum>& ptsD);
};

#endif
