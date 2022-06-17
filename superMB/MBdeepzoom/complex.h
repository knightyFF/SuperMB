//minimal implementation
#ifndef MY_COMPLEX_TYPE
#define MY_COMPLEX_TYPE
template<typename T>
class fcomplex{
public:
    fcomplex(){re=T(0.); im=T(0.);}
    fcomplex(T rv){re=rv; im=T(0.);}
    fcomplex(T rv, T iv){re=rv; im=iv;}
	template<typename U> fcomplex(U rv, U iv){re=T(rv); im=T(iv);}
	fcomplex(const fcomplex<T>& cv){re=cv.re; im=cv.im;}
    template<typename U> fcomplex(const fcomplex<U>& cv){re=T(cv.real()); im=T(cv.imag());}
    fcomplex<T>& operator=(const fcomplex<T>& cv){re=cv.re; im=cv.im; return *this;}
    fcomplex<T>& operator+=(const fcomplex<T>& cv){re+=cv.re; im+=cv.im; return *this;}
    fcomplex<T>& operator-=(const fcomplex<T>& cv){re-=cv.re; im-=cv.im; return *this;}
	fcomplex<T>& operator*=(const fcomplex<T>& cv){T re1 = re*cv.re-im*cv.im; im=re*cv.im+im*cv.re; re=re1; return *this;}

    T& real(){return re;}
    const T& real()const {return re;}
    T& imag(){return im;}
    const T& imag()const {return im;}
public:
    T re,im;
};
// Operators
template<typename T>
inline fcomplex<T> operator+(const fcomplex<T> &c1, const fcomplex<T>& c2){
    return fcomplex<T>(c1.re+c2.re, c1.im+c2.im);
}

template<typename T>
inline fcomplex<T> operator+(T c1, const fcomplex<T>& c2){
    return fcomplex<T>(c1+c2.re, c2.im);
}

template<typename T>
inline fcomplex<T> operator+(const fcomplex<T> &c1, T c2){
    return fcomplex<T>(c1.re+c2, c1.im);
}

template<typename T>
inline fcomplex<T> operator-(const fcomplex<T> &c1, const fcomplex<T>& c2){
    return fcomplex<T>(c1.re-c2.re, c1.im-c2.im);
}

template<typename T, typename U>
inline fcomplex<T> operator-(const fcomplex<T> &c1, const fcomplex<U>& c2){
    return fcomplex<T>(c1.re-c2.re, c1.im-c2.im);
}

template<typename T>
inline fcomplex<T> operator-(const fcomplex<T> &c){
    return fcomplex<T>(-c.re, -c.im);
}

template<typename T>
inline fcomplex<T> operator*(const fcomplex<T>& c1, const fcomplex<T>& c2){
    return fcomplex<T>(c1.re*c2.re-c1.im*c2.im, c1.re*c2.im+c1.im*c2.re);
}

template<typename T>
inline fcomplex<T> operator*(T c1, const fcomplex<T>& c2){
    return fcomplex<T>(c1*c2.re, c1*c2.im);
}

template<typename T>
inline fcomplex<T> operator*(const fcomplex<T>& c2, T c1){
    return fcomplex<T>(c1*c2.re, c1*c2.im);
}

template<typename T>
inline fcomplex<T> operator/(const fcomplex<T> &c1, const fcomplex<T>& c2){
    return T(1.)/norm(c2) * c1*conjugate(c2);
}

template<typename T>
inline std::ostream& operator<<(std::ostream& os, const fcomplex<T> &c){
    return os << "( " << c.real() << " , " << c.imag() << " )";
}

//Functions
template<typename T>
inline fcomplex<T> conjugate(const fcomplex<T> &c){
    return fcomplex<T>(c.re, -c.im);
}

template<typename T>
inline fcomplex<T> sqr(const fcomplex<T> &c){
    return fcomplex<T>(c.re*c.re-c.im*c.im, T(2.)*c.re*c.im);
}

template<typename T>
inline T norm(const fcomplex<T> &c){
    return c.re*c.re+c.im*c.im;
}

template<typename T>
inline T norm0(const fcomplex<T> &c){
    return std::abs(c.re)+std::abs(c.im);
}

template<typename T>
inline T normi(const fcomplex<T> &c){
    return std::max(std::abs(c.re), std::abs(c.im));
}

template<typename T>
inline T abs(const fcomplex<T> &c){
    return std::hypotl(c.re,c.im);//sqrt(norm(c));
}

template<typename T>
inline T real(const fcomplex<T> &c){
    return c.real();
}

template<typename T>
inline T imag(const fcomplex<T> &c){
    return c.imag();
}
#endif // COMPLEX

