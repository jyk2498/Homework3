#include "GF_field.hh"
#include <cmath>
#include<iostream>
#include<iomanip>
#include<cassert>

void GF_FIELD::print_alphaTobin() 
{
    for(const int& elem : alphaTobinary) std::cout << std::setw(3) << elem << " ";
    std::cout << std::endl;
}

void GF_FIELD::print_binToalpha() 
{
    for(const int& elem : binaryToalpha) std::cout << std::setw(3) << elem << " ";
    std::cout << std::endl;
}

GF_FIELD::GF_FIELD(int primitivePoly, int symbolSize) 
    : primitivePolynomial(primitivePoly), 
      symbolSize(symbolSize)  
{
    fieldSize = (int)pow(2, symbolSize);
    alphaTobinary.resize((int)pow(2, symbolSize));
    binaryToalpha.resize((int)pow(2, symbolSize));

    int checkingMask = (0x01 << symbolSize);
    int tempVal = 0x01;
    for(int i = 0; i < ((int)pow(2, symbolSize) - 1); i++) {
        alphaTobinary[i] = tempVal;
        binaryToalpha[tempVal] = i;

        tempVal = (tempVal << 1);
        if(tempVal & checkingMask) { // very first bit is 1 -> modulo operation
            tempVal ^= primitivePoly;
        }
    }

    alphaTobinary[(int)pow(2, symbolSize) - 1] = 0;
    binaryToalpha[0] = -1;
}

int GF_FIELD::alpha_to_binary(int exp) const
{
    if(exp == -1) return 0;
    else return alphaTobinary[exp];
}

int GF_FIELD::binary_to_alpha(int bin) const
{
    return binaryToalpha[bin];
}

GF_ELEMENT::GF_ELEMENT(int exp, int binary, GF_FIELD* basefield, int fieldSize)
    : exp(exp), binary(binary), basefield(basefield), fieldSize(fieldSize) 
{}

void GF_ELEMENT::set_value(int exp_, int binary_)
{
    exp = exp_;
    binary = binary_;
}

int GF_ELEMENT::pow(int operandexp)
{
    return ((exp * operandexp) % (fieldSize - 1));
}

GF_ELEMENT GF_ELEMENT::operator+(const GF_ELEMENT& b) const
{
    int temp_result = binary ^ b.get_binary();
    return GF_ELEMENT(basefield->binary_to_alpha(temp_result), temp_result, basefield, basefield->get_fieldsize()); 
}

GF_ELEMENT GF_ELEMENT::operator*(const GF_ELEMENT& b) const
{
    int temp_exp = ((exp + b.get_exp()) % (fieldSize - 1));
    return GF_ELEMENT(temp_exp, basefield->alpha_to_binary(temp_exp), basefield, basefield->get_fieldsize());
}

std::ostream& operator<<(std::ostream& os, const GF_ELEMENT& b)
{
    os << b.get_exp();
    // os << b.get_binary();
    return os;
}

void GF_ELEMENT::clear()
{
    binary = 0;
    exp = -1;
    fieldSize = -1;
    if(basefield != nullptr) delete basefield;
    basefield = nullptr;
}