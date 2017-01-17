#ifndef MODELDIVISION_H
#define MODELDIVISION_H

#include "model.h"

/// Classic radial model (center variable or fixed at origin)
class ModelDivision : public Model {
    std::vector<flnum> k_;
    flnum center_[2];
    bool varCenter_;
    bool even_;
public:
    ModelDivision(int order, bool varCenter=false, bool onlyEvenOrder=false);
    void change_order(int order);

    void apply(flnum& x, flnum& y) const;
    void evaluate(libNumerics::matrix<flnum>& ptsU,
                  libNumerics::matrix<flnum>& ptsD);
};

#endif
