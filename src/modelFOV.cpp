#include <cmath>
#include "modelFOV.h"
#include <cassert>
using namespace libNumerics;

/// Constructor. No term in r*r^2, so order cannot be 3.
ModelFOV::ModelFOV(int order)
: tanOmega_( tan(M_PI/6) ), k_(order, 0) {
    assert(order>=0 && order!=3);
    k_[0] = 1;
}

/// Change order, normally used for expansion
void ModelFOV::change_order(int order) {
    assert(order>=0 && order!=3);
    k_.resize(order, 0);
}

void ModelFOV::apply(flnum& x, flnum& y) const {
    const flnum r1 = sqrt(x*x+y*y);
    flnum r2_r1=0; // r2/r1, polynomial of degree order in r1
    std::vector<flnum>::const_reverse_iterator it=k_.rbegin();
    for(int i=k_.size()-1; it!=k_.rend(); ++it, --i) {
        if(i==1)
            r2_r1 *= r1;
        r2_r1 = *it + r1*r2_r1;
    }
    r2_r1 += (r1>1e-5? tan(r1*tanOmega_)/(r1*tanOmega_): 1);
    x = x*r2_r1;
    y = y*r2_r1;
}

/// For Levenberg-Marquardt minimization of the FOV model
class MinFOV : public MinLM {
    vector<flnum> ptsU_;
    vector<flnum> ptsD_;
public:
    MinFOV(const matrix<flnum>& ptsU, const matrix<flnum>& ptsD);

    flnum eval(flnum& tanOmega, std::vector<flnum>& k);

    virtual void modelData(const vector<flnum>& P,
                           vector<flnum>& ymodel) const;
    virtual void modelJacobian(const vector<flnum>& P,
                               matrix<flnum>& J) const;
};

/// Direct evaluation if fixed center, LM otherwise.
void ModelFOV::evaluate(matrix<flnum>& ptsU, matrix<flnum>& ptsD) {
    MinFOV(ptsU, ptsD).eval(tanOmega_, k_);
}

/// Constructor
MinFOV::MinFOV(const matrix<flnum>& ptsU, const matrix<flnum>& ptsD)
: ptsU_(2*ptsU.ncol()), ptsD_(2*ptsD.ncol()) {
    assert(ptsU.nrow()==2 && ptsD.nrow()==2);
    assert(ptsU.ncol()==ptsD.ncol());
    const int n=ptsU.ncol();
    for(int j=0, k=0; j<n; j++, k+=2) {
        ptsU_(k+0) = ptsU(0,j);
        ptsU_(k+1) = ptsU(1,j);
        ptsD_(k+0) = ptsD(0,j);
        ptsD_(k+1) = ptsD(1,j);
    }
}

/// Evaluation
flnum MinFOV::eval(flnum& tanOmega, std::vector<flnum>& k) {
    int nCoeff = 1 + static_cast<int>(k.size());
    vector<flnum> P(nCoeff);
    P(0) = tanOmega;
    std::copy(k.begin(), k.end(), P.begin()+1);
    float err = minimize(P, ptsD_, 1.0e-17, 1000);
    tanOmega = P(0);
    std::copy(P.begin()+1, P.end(), k.begin());
    return err/M_SQRT2;
}

/// P contains tan(omega) followed by the coefficients of the polynomial
void MinFOV::modelData(const vector<flnum>& P, vector<flnum>& ymodel) const {
    const int order = P.nrow();
    for(int i=0; i<ptsU_.nrow(); i+=2) {
        const flnum x=ptsU_(i+0), y=ptsU_(i+1);
        const flnum r = sqrt(x*x+y*y);
        flnum r2_r1=0; // r2/r1
        for(int j=order-2; j>=0; j--) {
            if(j==1)
                r2_r1 *= r;
            r2_r1 = P(j+1) + r*r2_r1;
        }
        if(r*P(0) >= 1e-5)
            r2_r1 += tan(r*P(0))/(r*P(0));
        ymodel(i+0) = x*r2_r1;
        ymodel(i+1) = y*r2_r1;
    }
}

/// P contains tan(omega) followed by the coefficients of the polynomial
void MinFOV::modelJacobian(const vector<flnum>& P, matrix<flnum>& J) const {
    const int order = P.nrow()-1;
    for(int i=0; i<ptsU_.nrow(); i+=2) {
        const flnum x=ptsU_(i+0), y=ptsU_(i+1);
        const flnum r = sqrt(x*x+y*y);
        if(r*P(0)<1e-10)
            J(i+0,0) = J(i+1,0) = 0;
        else {
            flnum d = tan(r*P(0));
            d = (r*P(0)*(1+d*d)-d)/(r*P(0)*P(0));
            J(i+0,0) = x*d;
            J(i+1,0) = y*d;
        }

        flnum rj=1; // rj: j-th power of r
        for(int j=0; j<order; j++, rj*=r) {
            J(i+0,j+1) = x*rj;
            J(i+1,j+1) = y*rj;
            if(j==1) // Skip r2
                rj *= r;
        }
    }
}
