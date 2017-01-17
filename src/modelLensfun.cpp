#include "modelLensfun.h"

/// Apply distortion to point (x,y)
void ModelLensfun::apply(flnum& x, flnum& y) const {
    flnum ru = sqrt(x*x+y*y);
    flnum rd = k0_ + k1_*ru + k2_*ru*ru + k3_*ru*ru*ru + k4_*ru*ru*ru*ru;
    x *= rd;
    y *= rd;
}
