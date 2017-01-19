#include "modelLensfun.h"
#include "modelPoly.h"
#include "modelRadial.h"
#include "modelRational.h"
#include "modelDivision.h"
#include "modelFOV.h"
#include "cmdLine.h"
#include <numeric>
#include <functional>

/// Number of points per line and per column
const int nbPoints=20;
using libNumerics::flnum;

/// Coordinate with function create_model.
const char* modeles[] = {
    "radial",
    "radial_center",
    "radial_odd",
    "radial_center_odd",

    "division",
    "division_center",
    "division_even",
    "division_center_even",

    "FOV",

    "poly",

    "rational"
};

/// Create model corresponding to its name
Model* create_model(std::string name, int order) {
    if(name == "radial")
        return new ModelRadial(order, false, false);
    if(name == "radial_center")
        return new ModelRadial(order, true, false);
    if(name == "radial_odd")
        return new ModelRadial(order, false, true);
    if(name == "radial_center_odd")
        return new ModelRadial(order, true, true);
    if(name == "poly")
        return new ModelPoly(order);
    if(name == "division")
        return new ModelDivision(order, false, false);
    if(name == "division_center")
        return new ModelDivision(order, true, false);
    if(name == "division_even")
        return new ModelDivision(order, false, true);
    if(name == "division_center_even")
        return new ModelDivision(order, true, true);
    if(name == "FOV")
        return new ModelFOV(order);
    if(name == "rational")
        return new ModelRational(order);
    return 0;
}

/// Display list of available models for user
void list_models() {
    std::cerr << " list of models:";
    int nModels = sizeof(modeles)/sizeof(modeles[0]);
    for(int i=0; i<nModels; i++)
        std::cerr << "  " << modeles[i];
    std::cerr << std::endl;
}

/// List of integers separated by commas
std::vector<int> decode_orders(std::string list) {
    std::stringstream s; s<<list; 
    std::vector<int> orders;
    while(true) {
        int order=0; s>>order;
        if(!order) break;
        orders.push_back(order);
        char c; s>>c; // Eat the separator (normally a comma)
    }
    return orders;
}

/// Read floating point parameter following field in name
bool decode_param(std::string name, std::string field, flnum& p) {
    size_t pos = name.find(field);
    if(pos == std::string::npos)
        return false;
    std::stringstream s; s<<name.substr(pos+field.size());
    return (s>>p);
}

/// Read a line of lensfun xml database to create appropriate model
Model* decode_model(std::string name) {
    std::string field="model=";
    size_t pos = name.find(field);
    if(pos == std::string::npos)
        return 0;
    std::stringstream s(name.substr(pos+field.size()));
    std::string str; s>>str;
    if(str == "ptlens") {
        flnum a=0, b=0, c=0;
        if(! decode_param(name, "a=", a) ||
           ! decode_param(name, "b=", b) ||
           ! decode_param(name, "c=", c))
            std::cerr << "An expected field a, b, or c could not be decoded in "
                      << "string:\n" << name << std::endl;
        return new ModelLensfun(a, b, c);
    }
    if(str == "poly3" || str == "poly5") {
        flnum k1=0;
        if(! decode_param(name, "k1=", k1))
            return 0;
        if(str == "poly3")
            return new ModelLensfun(k1);
        flnum k2=0;
        if(! decode_param(name, "k2=", k2))
            return 0;
        return new ModelLensfun(k1,k2);
    }
    std::cerr << "Unknown model: " << str << std::endl;
    return 0;
}

/// Synthesize points on a grid of nbPoints x nbPoints pixels with coordinates
/// in [-1,1].
void synthesize_points(libNumerics::matrix<flnum>& ptsU,
                       libNumerics::matrix<flnum>& ptsD,
                       flnum dx, flnum dy,
                       const Model& m=ModelLensfun()) {
    assert(ptsU.nrow()==2 && ptsU.ncol()==nbPoints*nbPoints);
    assert(ptsD.nrow()==2 && ptsD.ncol()==nbPoints*nbPoints);
    const flnum delta=2.0/nbPoints;
    for(int i=0, n=0; i<nbPoints; i++)
        for(int j=0; j<nbPoints; j++) {
            ptsU(0,n) = -1+dx+delta*i;
            ptsU(1,n) = -1+dy+delta*j;
            ++n;
        }
    ptsD = ptsU;
    m.applyMat(ptsD);
}

void compare(const libNumerics::matrix<flnum>& A, libNumerics::matrix<flnum> B){
    B -= A;
    // Take square of each element
    std::transform(B.begin(), B.end(), B.begin(), B.begin(),
                   std::multiplies<flnum>());
    //    std::cout << B << std::endl;
    std::cout << "RMSE/max: "
              << sqrt(std::accumulate(B.begin(), B.end(), 0.0)
                      /(B.nrow()*B.ncol())) << ' ' 
              << sqrt(*std::max_element(B.begin(), B.end())) << std::endl;
}

int main(int argc, char* argv[]) {
    CmdLine cmd;
    bool reverse=false;
    std::string distort_model;
    cmd.add( make_option('r',reverse,"reverse").doc("Find inverse model") );
    cmd.add( make_option('d',distort_model,"distort")
             .doc("XML line of lensfun db") );
    try { cmd.process(argc, argv); }
    catch(std::string str) {
        std::cerr << "Error: " << str << std::endl;
        return 1;
    }
    if(argc!=3) {
        std::cerr <<"Usage: " <<argv[0] <<" "
                  << "[options] model order,order,..." << std::endl;
        cmd.prefixDoc = "\t";
        std::cerr << "Options:\n" << cmd;
        list_models();
        return 1;
    }

    std::vector<int> orders = decode_orders(argv[2]);
    if(orders.empty()) {
        std::cerr << "Invalid orders, must be a comma separated list of "
                  << "positive integers" << std::endl;
        return 1;
    }

    Model* model = create_model(argv[1], orders.front());
    if(! model) {
        list_models();
        return 1;
    }

    Model* m=0;
    if(distort_model.empty()) {
        flnum a=0.013182, b=-0.036019, c=0;
        m = new ModelLensfun(a, b, c);
    } else
        m = decode_model(distort_model);
    if(!m) {
        std::cerr << "Could not decode distortion model" << std::endl;
        return 1;
    }

    const int n=nbPoints*nbPoints;

    // Synthesize (undistorted,distorted) point pairs
    libNumerics::matrix<flnum> ptsU(2,n), ptsD(2,n);

    for(std::vector<int>::const_iterator it=orders.begin();
        it!=orders.end(); ++it) {
        model->change_order(*it);
        synthesize_points(ptsU, ptsD, 1.0/nbPoints, 1.0/nbPoints, *m);
        if(reverse)
            swap(ptsU, ptsD);
        model->evaluate(ptsU, ptsD);
        model->applyMat(ptsU);

        // Evaluate on independent data
        synthesize_points(ptsU, ptsD, 0, 0, *m);
        if(reverse)
            swap(ptsU, ptsD);
        model->applyMat(ptsU);
        std::cout << *it << ' '; compare(ptsU, ptsD);
    }

    delete m;
    delete model;
    return 0;
}
