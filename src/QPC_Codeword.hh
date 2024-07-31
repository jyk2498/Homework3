#ifndef QPC_CODEWORD_HH
#define QPC_CODEWORD_HH
    #include "GF_field.hh"
    #include "QPC_config.hh"
    #include <map>
    #include "ELP.hh"
    #include <vector>

    class QPC_Codeword {
        public:
            void calculate_syndrome();
            void set_err_value(int symbolPos, int exp);
            void set_err_value_bin(int symbolPos, int exp);
            void add_err_value(int symbolPos, int binary);
            int get_errnum() {return errNum;}

            Decodingresult decode();
            int* get_decoded_positions();
            void clear();

            void print_errInfo();
            void print_codewordInfo();
            void print_decodedInfo();
            void print_hmatrix();
            void print_syndromeInfo();

            static void basefiled_initialize();
            static void QPC_hmatrix_init();
            static void QPC_dealloc();
            
            friend class MemTransferBlock;

            QPC_Codeword();
        private:
            static GF_FIELD*                                                    QPCBaseField;
            static std::array<std::array<GF_ELEMENT, QPC_N>, (QPC_N - QPC_K)>*  QPC_hmatrix;

            ELP                                                                 ErrevaLPoly;
            ELP                                                                 ErrLocationPoly;

            std::array<GF_ELEMENT, QPC_N>                                       codeword;
            std::array<GF_ELEMENT, QPC_N - QPC_K>                               syndrome;

            int                                                                 errNum;
            std::map<int, int>                                                  errInfo;

            int                                                                 decodedErrNum;
            std::map<int, int>                                                  decodedErrInfo;
            std::vector<int>                                                    decodedPos;

            void apply_euclidean();
            void find_location_chien();
            void find_err_value_forney();
    };
#endif