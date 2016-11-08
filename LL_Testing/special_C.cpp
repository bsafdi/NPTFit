#include "/group/hepheno/heptools/boost_1_61_0/boost/python.hpp"
#include <iostream>
#include <ctime>
#include <cmath>
#include "/group/hepheno/heptools/boost_1_61_0/boost/math/special_functions/gamma.hpp"
#include "/group/hepheno/heptools/boost_1_61_0/boost/multiprecision/float128.hpp"


using namespace std;
using namespace boost::multiprecision;
using namespace boost::math;


//int sp_max = 5000;

std::vector<double> gammainc_up_fct_ary(signed long m_max, double z, double a) {
    std::vector<double> vec(m_max + 1 );
    int index = 0 ;
    while (index < m_max +1 ) {
        double val = gamma_q(z + index, a);
        std::cout << val << " \n";
        vec.push_back(val); //
        index += 1 ;
    }
    return vec;
}

std::vector<double> gammainc_lo_fct_ary(signed long m_max, double z, double a) {
    std::vector<double> vec(m_max + 1 );
    int index = 0 ;
    while (index < m_max + 1) {
        double val = gamma_p(z + index, a);
        std::cout << val << " \n";
        vec.push_back(val); //
        index += 1 ;
    }
    return vec;
}

double log_factorial(int k){
    return lgamma(k+1); // std::
}

int main(){
    std::cout << "C++ method loaded " << std::endl;
};


BOOST_PYTHON_MODULE(special){
    using namespace boost::python;
    def("gammainc_up_fct_ary", gammainc_up_fct_ary);
    def("log_factorial",log_factorial);
    def("gammainc_lo_fct_ary",gammainc_lo_fct_ary);
}
