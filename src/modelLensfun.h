#ifndef MODELLENSFUN_H
#define MODELLENSFUN_H

#include "model.h"

class ModelLensfun : public Model {
public:
    flnum k0_, k1_, k2_, k3_, k4_;
    /// Identity constructor
    ModelLensfun(): k0_(1), k1_(0), k2_(0), k3_(0), k4_(0) {}
    /// poly3 constructor
    ModelLensfun(flnum k1): k0_(1-k1), k1_(0), k2_(k1), k3_(0), k4_(0) {}
    /// poly5 constructor
    ModelLensfun(flnum k1,flnum k2): k0_(1), k1_(0), k2_(k1), k3_(0), k4_(k2) {}
    /// ptlens constructor
    ModelLensfun(flnum a, flnum b, flnum c)
    : k0_(1-a-b-c), k1_(c), k2_(b), k3_(a), k4_(0) {}

    void apply(flnum& x, flnum& y) const;
};

#endif
