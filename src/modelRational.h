#ifndef MODELRATIONAL_H
#define MODELRATIONAL_H

#include "model.h"

/// Bivariate rational model
class ModelRational : public Model {
    int order_;
    std::vector<flnum> a_, b_, c_;
    mutable std::vector<flnum> xi_, yi_; ///< Temporary storage used frequently
public:
    explicit ModelRational(int order);
    void change_order(int order);

    void apply(flnum& x, flnum& y) const;
    void evaluate(libNumerics::matrix<flnum>& ptsU,
                  libNumerics::matrix<flnum>& ptsD);
};

#endif
