#ifndef ELP_HH
#define ELP_HH
    #include"GF_field.hh"
    #include<vector>

    class ELP {
        public:
            ELP();
            ELP(int maxDegree_, GF_FIELD* ELPBaseField);
            ~ELP();

            void print_poly(); 
            void change_coeff(int nthTerm, GF_ELEMENT* coefficient);
            int get_degree() const {return degree;}
            int get_coeff_exp(int nthTerm) const {return coefficients[nthTerm].get_exp();}
            int get_coeff_binary(int nthTerm) const {return coefficients[nthTerm].get_binary();}
            int get_highest_coeff_exp() const {return coefficients[degree].get_exp();}
            void set_degree(int degree_) {degree = degree_;}

            void set_maxdegree(int maxDegree_);
            void set_field(GF_FIELD* ELPBaseField_) {ELPBaseField = ELPBaseField_;}
            void clear();

            friend std::ostream& operator<<(std::ostream& os, ELP& otherelp);

            //ELP* add_elp(ELP* otherelp);
            //ELP* multiply_elp(ELP* otherelp);
            ELP operator+(const ELP& otherelp) const;
            ELP operator*(const ELP& otherelp) const;
            ELP derivative();
            int assign(GF_ELEMENT* term);
        private: 
            int                     degree;
            int                     maxDegree;
            GF_FIELD*               ELPBaseField;

            std::vector<GF_ELEMENT> coefficients;
    }; 
#endif