#ifndef GF_FIELD_HH
#define GF_FIELD_HH
    #include<array>
    #include<vector>
    #include"ECC.hh"

    class GF_FIELD {
        public:
            void print_alphaTobin();
            void print_binToalpha();

            int alpha_to_binary(int exp) const;
            int binary_to_alpha(int bin) const;

            int get_fieldsize() {return fieldSize;}

            GF_FIELD(int primitive_poly, int symbolSize);
            GF_FIELD() : primitivePolynomial(-1), symbolSize(-1) {} 
        private:
            int primitivePolynomial; 
            int symbolSize;
            int fieldSize;
            std::vector<int> alphaTobinary;
            std::vector<int> binaryToalpha;
    };

    class GF_ELEMENT {
        public:
            GF_ELEMENT(int exp, int binary, GF_FIELD* basefield, int fieldSize);
            GF_ELEMENT() : exp{-1}, binary{0}, basefield(nullptr) {}
            GF_ELEMENT(int exp, int binary) : exp{exp}, binary{binary}, basefield{nullptr}, fieldSize{-1} {}
            
            int get_binary() const {return binary;} 
            int get_exp() const {return exp;}
            void set_value(int exp_, int binary_);
            void set_field(GF_FIELD* basefield) {basefield = basefield;}
            int pow(int operandexp);
            GF_FIELD* get_basefiled() {return basefield;}
            void clear();

            GF_ELEMENT operator+(const GF_ELEMENT& b) const;
            GF_ELEMENT operator*(const GF_ELEMENT& b) const;

            friend std::ostream& operator<<(std::ostream& os, const GF_ELEMENT& b); 

        private:
            int binary;
            int exp;
            int fieldSize;
            GF_FIELD* basefield;
    }; 
#endif