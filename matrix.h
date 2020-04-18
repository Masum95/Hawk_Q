#ifndef DEF_MATRIX
#define DEF_MATRIX

#include <iostream>
#include <sstream>
#include <vector>
#include <stdexcept>

template <class T>
class Matrix{
    private:
        std::vector<std::vector<T> > array;
        int numRow;
        int numCol;

    public:
        Matrix<T>(int numRow, int numCol);
        Matrix<T>(std::vector<std::vector<T> > const &array);
        Matrix<T>();

        int getnumRow() const;
        int getnumCol() const;

        Matrix<T> add(const Matrix<T>& m) const;
        Matrix<T> subtract(const Matrix<T>& m) const;
        Matrix<T> multiply(const Matrix<T>& m) const;
        Matrix<T> dot(const Matrix<T>& m) const;
        Matrix<T> transpose() const;
        Matrix<T> multiply(const T& value) const;
        Matrix<T> divide(const T& value) const;

        Matrix<T> applyFunction(T (*function)(T)) const;
        Matrix<T> subMat(int startH, int startW, int h, int w) const;

        void fill(const T& value);
        void put(int h, int w, const T& value);
        T get(int h, int w) const;

        void print(std::ostream &flux) const;

        bool operator==(const Matrix<T>& m);
        bool operator!=(const Matrix<T>& m);
        Matrix<T> operator+=(const Matrix<T>& m);
        Matrix<T> operator-=(const Matrix<T>& m);
        Matrix<T> operator*=(const Matrix<T>& m);
        Matrix<T> operator*=(const T &m);
        Matrix<T> operator/=(const T &m);
        T& operator()(int y, int x);
};

template <class T> Matrix<T> operator+(const Matrix<T>& a, const Matrix<T>& b);
template <class T> Matrix<T> operator-(const Matrix<T>& a, const Matrix<T>& b);
template <class T> Matrix<T> operator*(const Matrix<T>& a, const Matrix<T>& b);
template <class T> Matrix<T> operator*(const T &b, const Matrix<T>& a);
template <class T> Matrix<T> operator/(const Matrix<T>& a, const T &b);
template <class T> std::ostream& operator<<(std::ostream &flux, const Matrix<T>& m);

#endif


template <class T>
Matrix<T>::Matrix(int numRow, int numCol){
    this->numRow = numRow;
    this->numCol = numCol;
    this->array = std::vector<std::vector<T> >(numRow, std::vector<T>(numCol));
}

template <class T>
Matrix<T>::Matrix(std::vector<std::vector<T> > const &array){
    if(array.size()==0)
        throw std::invalid_argument("Size of array must be greater than 0.");

    this->numRow = array.size();
    this->numCol = array[0].size();
    this->array = array;
}

template <class T>
Matrix<T>::Matrix(){
    numRow = 0;
    numCol = 0;
}

template <class T>
int Matrix<T>::getnumRow() const{
    return numRow;
}

template <class T>
int Matrix<T>::getnumCol() const{
    return numCol;
}

template <class T>
void Matrix<T>::fill(const T& value){
    for (int i=0 ; i<numRow ; i++){
        for (int j=0 ; j<numCol ; j++){
            array[i][j] = value;
        }
    }
}

template <class T>
void Matrix<T>::put(int h, int w, const T& value){
    if(!(h>=0 && h<numRow && w>=0 && w<numCol))
        throw std::invalid_argument("Index out of bounds.");

    array[h][w] = value;
}

template <class T>
T Matrix<T>::get(int h, int w) const{
    if(!(h>=0 && h<numRow && w>=0 && w<numCol))
        throw std::invalid_argument("Index out of bounds.");

    return array[h][w];
}

template <class T>
Matrix<T> Matrix<T>::multiply(const T& value) const{
    Matrix result(array);
    for (int i=0 ; i<numRow ; i++){
        for (int j=0 ; j<numCol ; j++){
            result.array[i][j] *= value;
        }
    }

    return result;
}

template <class T>
Matrix<T> Matrix<T>::divide(const T& value) const{
    Matrix result(array);
    for (int i=0 ; i<numRow ; i++){
        for (int j=0 ; j<numCol ; j++){
            result.array[i][j] /= value;
        }
    }

    return result;
}

template <class T>
Matrix<T> Matrix<T>::add(const Matrix& m) const{
    if(!(numRow==m.numRow && numCol==m.numCol))
        throw std::invalid_argument("Matrix dimension must be the same.");

    Matrix result(numRow, numCol);
    for (int i=0 ; i<numRow ; i++){
        for (int j=0 ; j<numCol ; j++){
            result.array[i][j] = array[i][j] + m.array[i][j];
        }
    }

    return result;
}

template <class T>
Matrix<T> Matrix<T>::subtract(const Matrix& m) const{
    if(!(numRow==m.numRow && numCol==m.numCol))
        throw std::invalid_argument("Matrix dimension must be the same.");

    Matrix result(numRow, numCol);
    for (int i=0 ; i<numRow ; i++){
        for (int j=0 ; j<numCol ; j++){
            result.array[i][j] = array[i][j] - m.array[i][j];
        }
    }
    return result;
}

template <class T>
Matrix<T> Matrix<T>::multiply(const Matrix& m) const{
    if(!(numCol==m.numRow))
        throw std::invalid_argument("Dimension Mismatch.");

    Matrix result(numRow, m.numCol);

    
    for (int i = 0; i < numRow; i++)
    {
        for (int j = 0; j < m.numCol; j++)
        {
            result.array[i][j] = 0;
            for (int k = 0; k < m.numRow; k++)
            {
                result.array[i][j] += array[i][k] * m.array[k][j];
            }
        }
    }
    return result;
}

template <class T>
Matrix<T> Matrix<T>::dot(const Matrix& m) const{
    if(!(numCol==m.numRow))
        throw std::invalid_argument("Dot product not compatible.");

    T w=0;
    int mnumCol = m.numCol;

    Matrix<T> result(numRow, mnumCol);
    for (int i=0 ; i<numRow ; i++){
        for (int j=0 ; j<mnumCol ; j++){
            for (int h=0 ; h<numCol ; h++){
                w += array[i][h]*m.array[h][j];
            }
            result.array[i][j] = w;
            w=0;
        }
    }

    return result;
}

template <class T>
Matrix<T> Matrix<T>::transpose() const{
    Matrix<T> result(numCol, numRow);

    for (int i=0 ; i<numCol ; i++){
        for (int j=0 ; j<numRow ; j++){
            result.array[i][j] = array[j][i];
        }
    }
    return result;
}

template <class T>
Matrix<T> Matrix<T>::applyFunction(T (*function)(T)) const{
    Matrix<T> result(numRow, numCol);

    for (int i=0 ; i<numRow ; i++){
        for (int j=0 ; j<numCol ; j++){
            result.array[i][j] = (*function)(array[i][j]);
        }
    }

    return result;
}

template <class T>
Matrix<T> Matrix<T>::subMat(int startH, int startW, int h, int w) const{
    if(!(startH>=0 && startH+h<=numRow && startW>=0 && startW+w<=numCol))
        throw std::invalid_argument("Index out of bounds");

    Matrix<T> result(h,w);
    for (int i=startH ; i<startH+h ; i++){
        for (int j=startW ; j<startW+w ; j++){
            result.array[i-startH][j-startW] = array[i][j];
        }
    }
    return result;
}

template <class T>
void Matrix<T>::print(std::ostream &flux) const{
    int maxLength[numCol] = {};
    std::stringstream ss;

    for (int i=0 ; i<numRow ; i++){
        for (int j=0 ; j<numCol ; j++){
            ss << array[i][j];
            if(maxLength[j] < ss.str().size()){
                maxLength[j] = ss.str().size();
            }
            ss.str(std::string());
        }
    }

    for (int i=0 ; i<numRow ; i++){
        for (int j=0 ; j<numCol ; j++){
            flux << array[i][j];
            ss << array[i][j];
            for (int k=0 ; k<maxLength[j]-ss.str().size()+1 ; k++){
                flux << " ";
            }
            ss.str(std::string());
        }
        flux << std::endl;
    }
}

template <class T>
bool Matrix<T>::operator==(const Matrix& m){
    if(numRow==m.numRow && numCol==m.numCol){
        for (int i=0 ; i<numRow ; i++){
            for (int j=0 ; j<numCol ; j++){
                if(array[i][j]!=m.array[i][j]){
                    return false;
                }
            }
        }
        return true;
    }
    return false;
}

template <class T>
bool Matrix<T>::operator!=(const Matrix& m){
    return !operator==(m);
}

template <class T>
Matrix<T> Matrix<T>::operator+=(const Matrix& m){
    this->array = add(m).array;
    return *this;
}

template <class T>
Matrix<T> Matrix<T>::operator-=(const Matrix& m){
    this->array = subtract(m).array;
    return *this;
}

template <class T>
Matrix<T> Matrix<T>::operator*=(const Matrix& m){
    this->array = multiply(m).array;
    return *this;
}

template <class T>
Matrix<T> Matrix<T>::operator*=(const T &m){
    *this = this->multiply(m);
    return *this;
}

template <class T>
Matrix<T> Matrix<T>::operator/=(const T &m){
    *this = this->divide(m);
    return *this;
}

template <class T>
T& Matrix<T>::operator()(int y, int x){
    if(!(y>=0 && y<numRow && x>=0 && x<numCol))
        throw std::invalid_argument("Index out of bounds.");
    return array[y][x];
}

template <class T>
Matrix<T> operator+(const Matrix<T>& a, const Matrix<T>& b){
    return a.add(b);
}

template <class T>
Matrix<T> operator-(const Matrix<T>& a, const Matrix<T>& b){
    return a.subtract(b);
}

template <class T>
Matrix<T> operator*(const Matrix<T>& a, const Matrix<T>& b){
    return a.multiply(b);
}

template <class T>
Matrix<T> operator*(const T &b, const Matrix<T>& a){
    return a.multiply(b);
}
template <class T>
Matrix<T> operator/(const Matrix<T>& a, const T &b){
    return a.divide(b);
}

template <class T>
std::ostream& operator<<(std::ostream &flux, const Matrix<T>& m){
    m.print(flux);
    return flux;
}