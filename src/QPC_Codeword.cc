#include"QPC_Codeword.hh"
#include"QPC_config.hh"
#include"DRAM_Config.hh"
#include"ECC.hh"
#include<iostream>
#include<vector>

GF_FIELD* QPC_Codeword::QPCBaseField = new GF_FIELD;
std::array<std::array<GF_ELEMENT, QPC_N>, QPC_N - QPC_K>* QPC_Codeword::QPC_hmatrix = new std::array<std::array<GF_ELEMENT, QPC_N>, QPC_N - QPC_K>;

QPC_Codeword::QPC_Codeword() 
    : decodedErrNum{0}, errNum{0}
{
    ErrevaLPoly.set_maxdegree(QPC_MAX_POLY_DEGREE);
    ErrLocationPoly.set_maxdegree(QPC_MAX_POLY_DEGREE);
}

void QPC_Codeword::calculate_syndrome()
{
    for(int i = 0; i < (QPC_N - QPC_K); i++) {
        int tempBin = 0;
        for(int j = 0; j < QPC_N; j++) {
            int operand1 = codeword[j].get_exp();        
            int tempExp;

            if(operand1 == -1) continue;
            else {
                int operand2 = (*QPC_hmatrix)[i][j].get_exp();
                tempExp = ((operand1 + operand2) % (QPCBaseField->get_fieldsize() - 1));
            }

            tempBin ^= QPCBaseField->alpha_to_binary(tempExp);
        }
        syndrome[i].set_value(QPCBaseField->binary_to_alpha(tempBin), tempBin);
    }
}

void QPC_Codeword::set_err_value(int symbolPos, int exp)
{

    if(codeword[symbolPos].get_exp() == -1) {
        errNum++;
    }

    errInfo[symbolPos] = exp;
    codeword[symbolPos].set_value(exp, QPCBaseField->alpha_to_binary(exp));
}

void QPC_Codeword::set_err_value_bin(int symbolPos, int binary) 
{
    if(codeword[symbolPos].get_exp() == -1) errNum++;

    errInfo[symbolPos] = QPCBaseField->binary_to_alpha(binary);
    codeword[symbolPos].set_value(QPCBaseField->binary_to_alpha(binary), binary);
}

void QPC_Codeword::add_err_value(int symbolPos, int binary)
{
    if(binary == 0) return;
    else {
        int originalVal = codeword[symbolPos].get_binary();
        int temp = binary ^ originalVal;
        if(temp == 0) return;

        codeword[symbolPos].set_value(QPCBaseField->binary_to_alpha(temp), temp);
        if(originalVal == 0) errNum++;
        errInfo[symbolPos] = QPCBaseField->binary_to_alpha(temp);
    }
}

Decodingresult QPC_Codeword::decode()
{
    bool NEflag = true;
    for(int i = 0; i < QPC_N - QPC_K; i++) {
        if(syndrome[i].get_exp() != -1) {
            NEflag = false;
        }
    }
    if(NEflag) return NE;

    // location find
    apply_euclidean();
    find_location_chien();

    if((decodedErrNum <= 0) || (decodedErrNum) >= QPC_T) return DUE;
    if(ErrLocationPoly.get_degree() != decodedErrNum) return DUE;

    int chipPos[CHIP_NUM] = {0, };
    for(auto it = decodedPos.begin(); it != decodedPos.end(); it++) {
        int tempPos = (*it);

        if(tempPos >= QPC_N) return DUE; 
        chipPos[tempPos / PIN_CONFIG] += 1;
    }

    int tempSum = 0;
    for(int i = 0; i < CHIP_NUM; i++) {
        if(chipPos[i] != 0) tempSum++;
    }

    /*
    if(tempSum > 2) return DUE;
    else if(tempSum == 2) {
        for(int i = 0; i < CHIP_NUM; i++) {
            if(chipPos[i] >= 3) return DUE;
        }
    }
    */


    // value find
    find_err_value_forney();
    if(decodedErrInfo != errInfo) return SDC;
    return CE;
}

int* QPC_Codeword::get_decoded_positions()
{
    int* resultarr = new int[CHIP_NUM];
    for(const auto& elem : decodedErrInfo) {
        resultarr[elem.first / PIN_CONFIG] = 1;
    }
    return resultarr;
}

void QPC_Codeword::clear()
{
    for(auto& element : codeword) {
        element.set_value(-1, 0);
    }

    for(auto& element : syndrome) {
        element.set_value(-1, 0);
    }

    errNum = 0;
    errInfo.clear();

    decodedErrNum = 0;
    decodedErrInfo.clear();
    decodedPos.clear();

    ErrevaLPoly.clear();
    ErrLocationPoly.clear();
}

void QPC_Codeword::print_errInfo()
{
    std::cout << "errNum : " << errNum << std::endl;
    for(const auto& elem : errInfo) {
        std::cout << "pos : " << elem.first << ", val: " << elem.second << std::endl;
    }
    std::cout << std::endl;
}

void QPC_Codeword::print_codewordInfo()
{
    std::cout << "Codeword: " << std::endl;
    for(const auto& elem : codeword) {
        std::cout << elem << " ";
    }
    std::cout << std::endl;
}

void QPC_Codeword::print_decodedInfo()
{
    std::cout << "decodedErrNum : " << decodedErrNum << std::endl;
    for(const auto& elem : decodedErrInfo) {
        std::cout << "pos : " << elem.first << ", val : " << elem.second << std::endl;
    }
}

void QPC_Codeword::print_hmatrix()
{
    std::cout << "=== hmatrix ===" << std::endl;
    for(int i = 0; i < (QPC_N - QPC_K); i++) {
        for(int j = 0; j < QPC_N; j++) {
            std::cout << (*QPC_hmatrix)[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

void QPC_Codeword::print_syndromeInfo()
{
    std::cout << "syndrome : " << std::endl;
    for(int i = 0; i < (QPC_N - QPC_K); i++) {
        std::cout << syndrome[i].get_exp() << " ";
    }
    std::cout << std::endl;
}

void QPC_Codeword::basefiled_initialize()
{
    QPCBaseField = new GF_FIELD(QPC_PRIM_POLY, QPC_SYMBOL_SIZE);
}

void QPC_Codeword::QPC_hmatrix_init()
{
    for(int i = 0; i < (QPC_N - QPC_K); i++) {
        for(int j = 0; j < QPC_N; j++) {
            (*QPC_hmatrix)[i][j] = GF_ELEMENT(((i + 1) * j), QPCBaseField->alpha_to_binary((i + 1) * j), QPCBaseField, QPCBaseField->get_fieldsize());
        }
    }
}

void QPC_Codeword::QPC_dealloc()
{
    delete QPCBaseField;
    delete QPC_hmatrix;
}

void QPC_Codeword::apply_euclidean()
{
    GF_ELEMENT zero(-1, 0);
    GF_ELEMENT one(0, 1);

    ELP evaluator_1((QPC_N - QPC_K), QPCBaseField), evaluator_2((QPC_N - QPC_K), QPCBaseField); // r_(i - 2), r_(i - 1)
    evaluator_1.change_coeff((QPC_N - QPC_K), &one);
    for(int i = 0; i < (QPC_N - QPC_K); i++) evaluator_2.change_coeff(i, &syndrome[i]);

    ELP locator_1((QPC_N - QPC_K), QPCBaseField), locator_2((QPC_N - QPC_K), QPCBaseField); // t_(i - 2), t_(i - 1)
    locator_2.change_coeff(0, &one);
    
    while(true) {
        ELP quotient((QPC_N - QPC_K), QPCBaseField);

        // first calculate q_i and r_i by doing (r_(i - 2) / r_(i - 1)) ... (eval_1 / eval_2)
        GF_ELEMENT tempGFelem(-1, 0);
        ELP temp_quotient((QPC_N - QPC_K), QPCBaseField);
        while(evaluator_2.get_degree() <= evaluator_1.get_degree()) {
            int divident = (evaluator_1.get_highest_coeff_exp() - evaluator_2.get_highest_coeff_exp());
            while(divident < 0) divident += (QPC_FIELD_SIZE - 1);

            tempGFelem.set_value(divident, QPCBaseField->alpha_to_binary(divident));
            temp_quotient.change_coeff(evaluator_1.get_degree() - evaluator_2.get_degree(), &tempGFelem);

            evaluator_1 = (evaluator_1 + (temp_quotient * evaluator_2));
            quotient = (quotient + temp_quotient);

            tempGFelem.clear();
            temp_quotient.clear();
        } // now evaluator_1 is new remainder

        // swap
        ELP buffer = evaluator_1;
        evaluator_1 = evaluator_2;
        evaluator_2 = buffer;

        // calculate new t_i
        buffer = locator_2;
        locator_2 = locator_1 + (quotient * locator_2);
        locator_1 = buffer;

        if((evaluator_2.get_degree() < (QPC_N - QPC_K) / 2)) {
            ErrLocationPoly = locator_2;
            ErrevaLPoly = evaluator_2;
            return;
        }
    }
}

void QPC_Codeword::find_location_chien()
{
    GF_ELEMENT term(-1, 0, QPCBaseField, QPCBaseField->get_fieldsize());
    int s = 0;
    decodedErrNum = 0;

    for(int i = 1; i < QPCBaseField->get_fieldsize(); i++) {
        term.set_value(i, QPCBaseField->alpha_to_binary(i));
        if((ErrLocationPoly.assign(&term) == -1)) {
            int pos = (i - (QPCBaseField->get_fieldsize() - 1));
            decodedPos.push_back(-pos);
            decodedErrNum += 1;
        }
    }
}

void QPC_Codeword::find_err_value_forney()
{
    ELP derivative = ErrLocationPoly.derivative();    
    GF_ELEMENT term(-1, 0, QPCBaseField, QPCBaseField->get_fieldsize());
    for(auto elem : decodedPos) {
        int buffer = -elem;
        while(buffer < 0) buffer += (QPCBaseField->get_fieldsize() - 1);
        term.set_value(buffer, QPCBaseField->alpha_to_binary(buffer));

        int tempExp = (ErrevaLPoly.assign(&term) - derivative.assign(&term));
        while(tempExp < 0) tempExp += (QPCBaseField->get_fieldsize() - 1); 
        
        decodedErrInfo.insert(std::make_pair(elem, tempExp));
    }
}
