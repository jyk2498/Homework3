#include <iostream>
#include"ELP.hh"
#include<cassert>
#include<iomanip>

ELP::ELP()
    : degree{-1}, maxDegree{-1}
{

}

ELP::ELP(int maxDegree_, GF_FIELD* ELPBaseField) : degree{-1}, maxDegree{maxDegree_}, ELPBaseField{ELPBaseField} 
{
    coefficients.resize(maxDegree_ + 1);
    for(auto& elem : coefficients) {
        elem.set_value(-1, 0);
        elem.set_field(ELPBaseField);
    }
}

ELP::~ELP()
{
    // delete ELPBaseField; this field is from QPC_Codeword. cannot delete
}

void ELP::print_poly()
{
    std::cout << "ELP: ";
    for(int i = 0; i < degree; i++) {
        int currentExp = coefficients[i].get_exp();

        if(currentExp == -1) std::cout << " ";
        else std::cout << currentExp << " ";
    }
    std::cout << std::endl;
}

void ELP::change_coeff(int nthTerm, GF_ELEMENT* coefficient)
{
    assert(nthTerm <= maxDegree);

    coefficients[nthTerm].set_value(coefficient->get_exp(), coefficient->get_binary());
    for(int i = maxDegree; i >= 0; i--) {
        if(coefficients[i].get_exp() != -1) {
            degree = i;
            return;
        }
    }
}

void ELP::set_maxdegree(int maxDegree_)
{
    coefficients.resize(maxDegree_ + 1);
    maxDegree = maxDegree_;
}


void ELP::clear()
{
    degree = -1;
    for(auto& elem : coefficients) elem.set_value(-1, 0);
}

ELP ELP::operator+(const ELP& otherelp) const
{
    if(degree == -1) {
        ELP Result(maxDegree, ELPBaseField);
        GF_ELEMENT tempGFelem(-1, 0);

        Result.set_degree(otherelp.get_degree());
        for(int i = 0; i <= otherelp.get_degree(); i++) {
            tempGFelem.set_value(otherelp.get_coeff_exp(i), otherelp.get_coeff_binary(i));
            Result.change_coeff(i, &tempGFelem);
        }
        return Result;
    } else if(otherelp.get_degree() == -1) {
        ELP Result(maxDegree, ELPBaseField);
        GF_ELEMENT tempGFelem(-1, 0);

        Result.set_degree(degree);
        for(int i = 0; i <= degree; i++) {
            tempGFelem.set_value(this->get_coeff_exp(i), this->get_coeff_binary(i));
            Result.change_coeff(i, &tempGFelem);
        }
        return Result;
    }

    ELP tempResult(maxDegree, ELPBaseField);
    GF_ELEMENT tempGFelem(-1, 0);

    if(degree == (otherelp.get_degree())){ // same degree
        int lastNonzeroDegree = -1;

        for(int i = 0; i <= degree; i++) {
            int tempBin = (coefficients[i].get_binary() ^ otherelp.get_coeff_binary(i));
            if(tempBin != 0) lastNonzeroDegree = i;

            tempGFelem.set_value(ELPBaseField->binary_to_alpha(tempBin), tempBin);
            tempResult.change_coeff(i, &tempGFelem);
        }
        tempResult.set_degree(lastNonzeroDegree);
    } else { // degree is not same
        int biggerDegree;
        int smallerDegree;
        const ELP* biggerelp;

        if(degree > otherelp.get_degree()) {
            biggerDegree = degree;
            smallerDegree = otherelp.get_degree();
            biggerelp = this;
        } else {
            biggerDegree = otherelp.get_degree();
            smallerDegree = degree;
            biggerelp = &otherelp;
        }

        for(int i = 0; i <= smallerDegree; i++) {
            int tempBin = (coefficients[i].get_binary() ^ otherelp.get_coeff_binary(i));
            tempGFelem.set_value(ELPBaseField->binary_to_alpha(tempBin), tempBin);
            tempResult.change_coeff(i, &tempGFelem);
        }

        for(int i = smallerDegree + 1; i <= biggerDegree; i++) {
            tempGFelem.set_value((biggerelp->get_coeff_exp(i)), (biggerelp->get_coeff_binary(i)));
            tempResult.change_coeff(i, &tempGFelem);
        }
    }

    return tempResult;
}

ELP ELP::operator*(const ELP& otherelp) const
{
    if((degree == -1) || (otherelp.get_degree() == -1)) {
        ELP tempResult(maxDegree, ELPBaseField);
        return tempResult;
    }

    ELP tempResult(maxDegree, ELPBaseField);
    GF_ELEMENT tempGFelem(-1, 0);
    ELP tempElp(maxDegree, ELPBaseField);

    for(int i = 0; i <= degree; i++) {
        tempElp.clear();
        for(int j = 0; j <= otherelp.get_degree(); j++) {
            int myexp = get_coeff_exp(i);
            int otherexp = otherelp.get_coeff_exp(j);

            int tempExpResult;
            if((myexp == -1) || (otherexp == -1)) tempExpResult = -1;
            else tempExpResult = ((myexp + otherexp) % (ELPBaseField->get_fieldsize() - 1));

            tempGFelem.set_value(tempExpResult, ELPBaseField->alpha_to_binary(tempExpResult));
            tempElp.change_coeff(i + j, &tempGFelem);
        }
        tempResult = (tempResult + tempElp);
    }

    for(int i = maxDegree; i > 0; i--) {
        if(tempResult.get_coeff_exp(i) != -1) {
            tempResult.set_degree(i);
            break;
        }
    }
    return tempResult;
}

ELP ELP::derivative()
{
    ELP tempResult(maxDegree, ELPBaseField);
    if((degree == -1) || (degree == 0)) return tempResult;

    GF_ELEMENT tempGFelem(-1, 0);
    for(int i = 1; i <= degree; i += 2) {
        tempGFelem.set_value(coefficients[i].get_exp(), coefficients[i].get_binary());
        tempResult.change_coeff(i - 1, &tempGFelem);
    }

    for(int i = maxDegree; i > 0; i--) {
        if(tempResult.get_coeff_exp(i) != -1) {
            tempResult.set_degree(i);
            break;
        }
    }

    return tempResult;
}

int ELP::assign(GF_ELEMENT* term)
{
    if(this->get_degree() == -1) return -1;
    else if(this ->get_degree() == 0) return this->get_coeff_exp(0);

    int tempBin = this->get_coeff_binary(0);
    for(int i = 1; i <= degree; i++) {
        int coefficientExp = coefficients[i].get_exp();
        if(coefficientExp == -1) continue;

        tempBin ^= ELPBaseField->alpha_to_binary((coefficientExp + (term->pow(i))) % (ELPBaseField->get_fieldsize() - 1));
    }

    return ELPBaseField->binary_to_alpha(tempBin);
}

std::ostream& operator<<(std::ostream& os, ELP& otherelp)
{
    for(int i = 0; i <= otherelp.get_degree(); i++) {
        os << std::setw(3) << otherelp.get_coeff_exp(i) << " ";
    }
    return os;
}
