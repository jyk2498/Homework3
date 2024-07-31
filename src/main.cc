#include <iostream>
#include <stdio.h>
#include <vector>
#include <string.h>
#include <math.h>
#include <time.h>
#include <cstdlib>
#include <ctime>
#include <chrono>
#include <thread>
#include<random>
#include "GF_field.hh"
#include "ELP.hh"
#include "QPC_Codeword.hh"

void init();
void dealloc();

int main(void) {
    srand(time(NULL));
    init();

    QPC_Codeword mycodeword;
    mycodeword.print_hmatrix();

    mycodeword.set_err_value(8, 7);
    mycodeword.set_err_value(16, 21);
    mycodeword.set_err_value(25, 34);
    mycodeword.print_errInfo();

    mycodeword.calculate_syndrome();
    mycodeword.print_syndromeInfo();
    mycodeword.decode();
    mycodeword.print_decodedInfo();

    dealloc();
    return 0;
}

void init()
{
    QPC_Codeword::basefiled_initialize();
    QPC_Codeword::QPC_hmatrix_init();
}

void dealloc()
{
    QPC_Codeword::QPC_dealloc();
}