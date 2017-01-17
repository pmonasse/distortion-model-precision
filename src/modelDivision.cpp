#include "modelDivision.h"
#include <cassert>
using namespace libNumerics;

/// Constructor
ModelDivision::ModelDivision(int order, bool varCenter, bool onlyEvenOrder)
: k_(order+1, 0), varCenter_(varCenter), even_(onlyEvenOrder) {
    assert(order>=0);
    k_[0] = 1;
    center_[0]=center_[1]=0;
}

/// Change order, normally used for expansion
void ModelDivision::change_order(int order) {
    assert(order>=0);
    k_.resize(order+1, 0);
}

void ModelDivision::apply(flnum& x, flnum& y) const {
    x -= center_[0];
    y -= center_[1];
    const flnum r1 = sqrt(x*x+y*y);
    flnum r2_r1=0; // r2*r1, polynomial of degree order in r1
    for(std::vector<flnum>::const_reverse_iterator it=k_.rbegin();
        it!=k_.rend(); ++it)
        r2_r1 = *it + r1*r2_r1;
    x = center_[0] + x/r2_r1;
    y = center_[1] + y/r2_r1;
}

/// For Levenberg-Marquardt minimization of the division model
class MinDivision : public MinLM {
    vector<flnum> ptsU_;
    vector<flnum> ptsD_;
    bool center_;
    bool even_;
public:
    MinDivision(const matrix<flnum>& ptsU, const matrix<flnum>& ptsD);

    flnum eval(std::vector<flnum>& k, bool varCenter, flnum center[2],
               bool even);

    virtual void modelData(const vector<flnum>& P,
                           vector<flnum>& ymodel) const;
    virtual void modelJacobian(const vector<flnum>& P,
                               matrix<flnum>& J) const;
};

// Copy even index elements of v1 into v2
static void copy_step2_in(const std::vector<flnum>& v1, vector<flnum>& v2) {
    std::vector<flnum>::const_iterator in=v1.begin();
    vector<flnum>::iterator out=v2.begin();
    for(bool ok=true; in!=v1.end(); ++in, ok=!ok)
        if(ok) {
            assert(out != v2.end());
            *out++ = *in;
        }
}

// Copy v1 into even index elements of v2
static void copy_step2_out(vector<flnum>::const_iterator in,
                           vector<flnum>::const_iterator end,
                           std::vector<flnum>& v2) {
    std::vector<flnum>::iterator out=v2.begin();
    for(; in!=end; ++in, ++out) {
        assert(out != v2.end());
        *out++ = *in;
        if(out == v2.end()) // Do not go beyond v2.end()
            break;
    }
}

/// Direct evaluation if fixed center, LM otherwise.
void ModelDivision::evaluate(matrix<flnum>& ptsU, matrix<flnum>& ptsD) {
    MinDivision(ptsU, ptsD).eval(k_, varCenter_, center_, even_);
}

/// Constructor
MinDivision::MinDivision(const matrix<flnum>& ptsU, const matrix<flnum>& ptsD)
: ptsU_(2*ptsU.ncol()), ptsD_(2*ptsD.ncol()), center_(false), even_(false) {
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
flnum MinDivision::eval(std::vector<flnum>& k, bool varCenter, flnum center[2],
                        bool even) {
    center_ = varCenter;
    even_ = even;
    int nCoeff = static_cast<int>(even? (k.size()+1)/2: k.size());
    vector<flnum> P(nCoeff + (center_? 2: 0)); // +2: center of distortion
    if(even) copy_step2_in(k, P);
    else std::copy(k.begin(), k.end(), P.begin());
    if(center_) std::copy(center, center+2, P.begin()+nCoeff);
    float err = minimize(P, ptsD_, 1.0e-17, 1000);
    if(even) copy_step2_out(P.begin(), P.begin()+nCoeff, k);
    else std::copy(P.begin(), P.begin()+nCoeff, k.begin());
    if(center_) std::copy(P.begin()+nCoeff, P.end(), center);
    return err/M_SQRT2;
}

/// P contains the coefficients of the polynomial, followed by the center coords
/// if it is variable.
void MinDivision::modelData(const vector<flnum>& P,
                            vector<flnum>& ymodel) const {
    const int order = P.nrow()-(center_? 2: 0);
    const float cx=(center_? P(order): 0), cy=(center_? P(order+1): 0);
    for(int i=0; i<ptsU_.nrow(); i+=2) {
        const flnum x=ptsU_(i+0)-cx, y=ptsU_(i+1)-cy;
        const flnum r = even_? x*x+y*y: sqrt(x*x+y*y);
        flnum r2_r1=0; // r2*r1, polynomial of degree order-1 in r1
        for(int j=order-1; j>=0; j--)
            r2_r1 = P(j) + r*r2_r1;
        ymodel(i+0) = cx + x/r2_r1;
        ymodel(i+1) = cy + y/r2_r1;
    }
}

/// P contains the coefficients of the polynomial, followed by the center coords
/// if it is variable.
void MinDivision::modelJacobian(const vector<flnum>& P,
                              matrix<flnum>& J) const {
    const int order = P.nrow()-(center_? 2: 0);
    const float cx=(center_? P(order): 0), cy=(center_? P(order+1): 0);
    for(int i=0; i<ptsU_.nrow(); i+=2) {
        const flnum x=ptsU_(i+0)-cx, y=ptsU_(i+1)-cy;
        const flnum r = even_? x*x+y*y: sqrt(x*x+y*y);
        flnum d=0, rj=1; // rj: j-th power of r or of r2
        for(int j=0; j<order; j++, rj*=r)
            d += P(j)*rj;
        d = 1/d;
        rj=1;
        for(int j=0; j<order; j++, rj*=r) {
            J(i+0,j) = -x*rj*d*d;
            J(i+1,j) = -y*rj*d*d;
        }
        if(center_) {
            flnum n = 0;
            if(x*x+y*y>1e-10) {
                rj = 1/(x*x+y*y);
                for(int j=0; j<order; j++, rj*=r) {
                    int k=even_? 2*j: j;
                    n += P(j)*rj*k;
                }
            }
            n *= d*d;
            J(i+0,order+0)= 1-d+n*x*x;
            J(i+1,order+0)=     n*x*y;
            J(i+0,order+1)=     n*x*y;
            J(i+1,order+1)= 1-d+n*y*y;
        }
    }
}
