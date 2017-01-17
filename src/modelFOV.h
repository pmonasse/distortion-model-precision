#ifndef MODELFOV_H
#define MODELFOV_H

#include "model.h"

/// Classic radial model (fixed center at origin)
/// TODO: variable center
class ModelFOV : public Model {
    flnum tanOmega_;
    std::vector<flnum> k_;
public:
    ModelFOV(int order);
    void change_order(int order);

    void apply(flnum& x, flnum& y) const;
    void evaluate(libNumerics::matrix<flnum>& ptsU,
                  libNumerics::matrix<flnum>& ptsD);
};

#endif
