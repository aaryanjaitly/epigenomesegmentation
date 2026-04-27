#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <memory>
#include <iostream>
#include <stdexcept>
#include <initializer_list>

template <class T>
class Matrix 
{
    public:

    class col_iterator;
    class const_col_iterator;
    class row_iterator;
    class const_row_iterator;

    /**
     * @brief Construct empty matrix
     * 
     */
    Matrix() : Matrix(0, 0, T()){};

    /**
     * @brief Construct a empty matrix of size nrows x ncols
     * 
     * @param nrows 
     * @param ncols 
     */
    Matrix(std::size_t nrows, std::size_t ncols) : Matrix(nrows, ncols, T()) {}

    /**
     * @brief Construct a matrix of size nrows x ncols with inital value init
     * 
     * @param nrows 
     * @param ncols 
     * @param init 
     */
    Matrix(std::size_t nrows, std::size_t ncols, T init): nrows_(nrows), ncols_(ncols), data_(std::make_unique<T[]>(nrows * ncols))
    {
        for (std::size_t i = 0, end = nrows * ncols; i < end; ++i) 
        {
            data_[i] = init;
        }
    }

    /**
     * @brief Construct matrix from exisiting matrix (copy constructor)
     * 
     * @param m 
     */
    Matrix(const Matrix& m) : nrows_(m.nrows()), ncols_(m.ncols()), data_() 
    {
        data_ = std::make_unique<T[]>(nrows_ * ncols_);
        std::copy(m.data_.get(), m.data_.get() + (nrows_ * ncols_), data_.get());
    }
    
    /**
     * @brief Construct matrix from initializer list
     * 
     * @param elements 
     */
    Matrix(std::initializer_list<std::initializer_list<T>> elements)
    : nrows_(elements.size()), ncols_(0), data_() 
    {
        ncols_ = elements.begin()->size();
        data_ = std::make_unique<T[]>(nrows_ * ncols_);

        std::size_t i = 0;
        for (auto row : elements) 
        {
            for (auto e : row) 
            {
                data_[i] = e;
                ++i;
            }
        }
    }

    /**
     * @brief copy same elements from m to the matrix
     * 
     * @param m 
     * @return Matrix& 
     */
    Matrix& operator=(const Matrix& m) 
    {
        if (this == &m) 
        {
            return *this;
        }
        if (nrows_ * ncols_ != m.nrows() * m.ncols()) 
        {
            data_ = std::make_unique<T[]>(m.nrows() * m.ncols());
        }
        nrows_ = m.nrows();
        ncols_ = m.ncols();
        std::copy(m.data_.get(), m.data_.get() + (nrows_ * ncols_),
        data_.get());
        return *this;
    }

    /**
     * @brief get refernece to value at position (i,j)
     * 
     * @param i 
     * @param j 
     * @return value 
     */
    T& operator()(std::size_t i, std::size_t j) 
    {
        return data_[i * ncols_ + j];
    }

    /**
     * @brief get copy of value at position (i, j)
     * 
     * @param i 
     * @param j 
     * @return value
     */
    T operator()(std::size_t i, std::size_t j) const 
    {
        return data_[i * ncols_ + j];
    }

    /**
     * @brief return number of rows
     * 
     * @return std::size_t 
     */
    std::size_t nrows() const { return nrows_; }

    /**
     * @brief return number of columns
     * 
     * @return std::size_t 
     */
    std::size_t ncols() const { return ncols_; }

    /**
     * @brief return iterator to begin of column col_num
     * 
     * @param col_num 
     * @return col_iterator 
     */
    col_iterator col_begin(std::size_t col_num)
    {
        return col_iterator(&(this->operator()(0, col_num)), ncols_);     
    }

    /**
     * @brief return const iterator to begin of column col_num
     * 
     * @param col_num 
     * @return const_col_iterator 
     */
    const_col_iterator col_begin(std::size_t col_num) const
    {
        const T* p = &data_[0] + col_num;
        return const_col_iterator(p, ncols_);
    }

    /**
     * @brief return iterator to end of column col_num
     * 
     * @param col_num 
     * @return col_iterator 
     */
    col_iterator col_end(std::size_t col_num)
    {
        return col_iterator(&(this->operator()(nrows_-1, col_num))+ncols_, ncols_);
    }

    /**
     * @brief return const iterator to end of column col_num
     * 
     * @param col_num 
     * @return const_col_iterator 
     */
    const_col_iterator col_end(std::size_t col_num) const
    {
        const T* p = &data_[0] + (nrows_) * ncols_ + col_num;
        return const_col_iterator(p, ncols_);
    }

    /**
     * @brief return iterator to begin of row row_num
     * 
     * @param row_num 
     * @return row_iterator 
     */
    row_iterator row_begin(std::size_t row_num)
    {
        return row_iterator(&(this->operator()(row_num, 0)));
    }

    /**
     * @brief return const iterator to begin of row row_num
     * 
     * @param row_num 
     * @return const_row_iterator 
     */
    const_row_iterator row_begin(std::size_t row_num) const
    {
        const T* p = &data_[0] + row_num * ncols_;
        return const_row_iterator(p);
    }

    /**
     * @brief return iterator to end of row row_num
     * 
     * @param row_num 
     * @return row_iterator 
     */
    row_iterator row_end(std::size_t row_num)
    {
        return row_iterator(&(this->operator()(row_num, ncols_-1))+1);
    }

    /**
     * @brief return const iterator to end of row row_num
     * 
     * @param row_num 
     * @return const_row_iterator 
     */
    const_row_iterator row_end(std::size_t row_num) const
    {
        const T* p = &data_[0] + (row_num+1) * (ncols_);
        return const_row_iterator(p);
    }
    
    private:
    std::size_t nrows_;
    std::size_t ncols_;
    std::unique_ptr<T[]> data_;
};

/*
* Output matrix in row major ordering
*/
template <typename T>
std::ostream& operator<<(std::ostream& os, const Matrix<T>& m) 
{
    for (std::size_t i = 0; i < m.nrows(); ++i) 
    {
        for (std::size_t j = 0; j < m.ncols() - 1; ++j) 
        {
            os << m(i, j) << ' ';
        }
        os << m(i, m.ncols() - 1) << '\n';
    }
    return os;
}

template<typename T>
class Matrix<T>::col_iterator
{
    public:
        /**
         * define iterator traits
        */
        using iterator_category = std::forward_iterator_tag;
        using difference_type   = std::ptrdiff_t;
        using value_type        = T;
        using pointer           = T*;
        using reference         = T&;

        /**
         * Construct a new col iterator object
         */
        col_iterator(T* p, size_t n) : curr(p), ncols(n) {}

        /**
         * return pointer to next element
         */
        col_iterator& operator++() 
        {
            curr += ncols;
            return *this;
        }

        /**
         * return object the iterator is pointing to
         */
        T& operator*()
        {
            return *curr;
        }

        /**
         * compare if iterators are equal
         */
        bool operator==(const col_iterator& b) const
        {
            return curr == b.curr;
        }

        /**
         * compare if iterators are unequal
         */
        bool operator!=(const col_iterator& b) const
        {
            return curr != b.curr;
        }
    private:
        T* curr;
        size_t ncols;
};

template<typename T>
class Matrix<T>::const_col_iterator
{
    public:
        /**
         * define iterator traits
        */
        using iterator_category = std::forward_iterator_tag;
        using difference_type   = std::ptrdiff_t;
        using value_type        = T;
        using pointer           = T*;
        using reference         = T&;

        const_col_iterator(const T* p, size_t n) : curr(p), ncols(n){}

        /**
         * return pointer to next element
         */
        const_col_iterator& operator++() 
        {
            curr += ncols;
            return *this;
        }

        /**
         * return const reference to object the iterator is pointing to
         */
        const T& operator*() const
        {
            return *curr;
        }

        /**
         * compare if iterators are equal
         */
        bool operator==(const const_col_iterator& b) const
        {
            return curr == b.curr;
        }

        /**
         * compare if iterators are unequal
         */
        bool operator!=(const const_col_iterator& b) const
        {
            return curr != b.curr;
        }
    private:
        const T* curr;
        size_t ncols;
};

template<typename T>
class Matrix<T>::row_iterator
{
    public:
        /**
         * define iterator traits
        */
        using iterator_category = std::forward_iterator_tag;
        using difference_type   = std::ptrdiff_t;
        using value_type        = T;
        using pointer           = T*;
        using reference         = T&;

        /**
         * Construct a new row iterator object
         */
        row_iterator(T* p) : curr(p) {}

        /**
         * return pointer to next element
         */
        row_iterator& operator++() 
        {
            curr += 1;
            return *this;
        }

        /**
         * return object the iterator is pointing to
         */
        T& operator*()
        {
            return *curr;
        }

        /**
         * compare if iterators are equal
         */
        bool operator==(const row_iterator& b) const
        {
            return curr == b.curr;
        }

        /**
         * compare if iterators are unequal
         */
        bool operator!=(const row_iterator& b) const
        {
            return curr != b.curr;
        }
    private:
        T* curr;
};

template<typename T>
class Matrix<T>::const_row_iterator
{
    public:
        /**
         * define iterator traits
        */
        using iterator_category = std::forward_iterator_tag;
        using difference_type   = std::ptrdiff_t;
        using value_type        = T;
        using pointer           = T*;
        using reference         = T&;
        
        const_row_iterator(const T* p) : curr(p) {}

        /**
         * return pointer to next element
         */
        const_row_iterator& operator++() 
        {
            curr += 1;
            return *this;
        }

        /**
         * return const reference to object the iterator is pointing to
         */
        const T& operator*() const
        {
            return *curr;
        }

        /**
         * compare if iterators are equal
         */
        bool operator==(const const_row_iterator& b) const
        {
            return curr == b.curr;
        }

        /**
         * compare if iterators are unequal
         */
        bool operator!=(const const_row_iterator& b) const
        {
            return curr != b.curr;
        }
    private:
        const T* curr;
};

template <class T>
class Matrix3D
{
    public:

    enum dimension{FIRST, SECOND, THIRD};

    /**
     * @brief Construct empty matrix
     * 
     */
    Matrix3D() : Matrix3D(0, 0, 0, T()){};

    /**
     * @brief Construct a empty matrix of size nrows x ncols
     * 
     * @param nrows 
     * @param ncols 
     */
    Matrix3D(std::size_t dim1, std::size_t dim2, std::size_t dim3) : Matrix3D(dim1, dim2, dim3, T()) {}

    /**
     * @brief Construct a matrix of size nrows x ncols with inital value init
     * 
     * @param nrows 
     * @param ncols 
     * @param init 
     */
    Matrix3D(std::size_t dim1, std::size_t dim2, std::size_t dim3, T init): dim1_(dim1), dim2_(dim2), dim3_(dim3), data_(std::make_unique<T[]>(dim1 * dim2 * dim3))
    {
        std::size_t size = dim1 * dim2 * dim3;
        for (std::size_t i = 0; i < size; ++i) 
        {
            data_[i] = init;
        }
    }

    /**
     * @brief Construct matrix from exisiting matrix (copy constructor)
     * 
     * @param m 
     */
    Matrix3D(const Matrix3D& m) : dim1_(m.dim(FIRST)), dim2_(m.dim(SECOND)), dim3_(m.dim(THIRD)),data_() 
    {
        data_ = std::make_unique<T[]>(dim1_ * dim2_ * dim3_);
        std::copy(m.data_.get(), m.data_.get() + (dim1_ * dim2_ * dim3_), data_.get());
    }

    /**
     * @brief copy same elements from m to the matrix
     * 
     * @param m 
     * @return Matrix& 
     */
    Matrix3D& operator=(const Matrix3D& m) 
    {
        if (this == &m) 
        {
            return *this;
        }
        if (dim1_ * dim2_ * dim3_ != m.dim(FIRST) * m.dim(SECOND) * m.dim(THIRD)) 
        {
            data_ = std::make_unique<T[]>(m.dim(FIRST) * m.dim(SECOND) * m.dim(THIRD));
        }
        dim1_ = m.dim(FIRST);
        dim2_ = m.dim(SECOND);
        dim3_ = m.dim(THIRD);
        std::copy(m.data_.get(), m.data_.get() + (dim1_ * dim2_ * dim3_), data_.get());
        return *this;
    }

    /**
     * @brief get refernece to value at position (i,j,k)
     * 
     * @param i 
     * @param j 
     * @param k
     * @return value 
     */
    T& operator()(std::size_t i, std::size_t j, std::size_t k) 
    {
        return data_[(i * dim2_ + j) * dim3_+ k];
    }

    /**
     * @brief get copy of value at position (i, j, k)
     * 
     * @param i 
     * @param j 
     * @param k
     * @return value
     */
    T operator()(std::size_t i, std::size_t j, std::size_t k) const 
    {
        return data_[(i * dim2_ + j) * dim3_+ k];
    }

    /**
     * @brief return size of dimension dim 
     * 
     * @return std::size_t 
     */
    std::size_t dim(dimension d) const 
    { 
        switch(d)
        {
            case (FIRST):
                return dim1_;
            case (SECOND):
                return dim2_;
            case (THIRD):
                return dim3_;
            default:
                throw std::invalid_argument("Unknown dimension");
        }
     };
    
    private:
    std::size_t dim1_;
    std::size_t dim2_;
    std::size_t dim3_;
    std::unique_ptr<T[]> data_;
};

#endif  // MATRIX_HPP
