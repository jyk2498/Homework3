#ifndef QPC_CONFIG_HH
#define QPC_CONFIG_HH
    #define QPC_N                       40
    #define QPC_K                       32
    #define QPC_T                       4 // (QPC_N - QPC_K) / 2
    #define QPC_CODEWORD_NUM_PER_BLOCK  2  

    #define QPC_SYMBOL_SIZE             8
    #define QPC_CODEWORD_SIZE           255

    #define QPC_FIELD_SIZE              256
    #define QPC_PRIM_POLY               0x15F

    #define QPC_MAX_POLY_DEGREE         8
#endif