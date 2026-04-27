#ifndef _MATRIX_READER_H
#define _MATRIX_READER_H

#include <iostream>
#include "matrix.h"

/**
 * @brief Reading a count matrix.
 */
class Reader
{
    public:
    Reader(){};

    /**
     * @brief Parses a matrix from the input stream.
     * 
     * @return matrix<int> 
     */
    Matrix<int> parse_matrix(std::istream&);

    /**
     * @brief Parses a starting indices from the input stream.
     * 
     * @return matrix<int> 
     */
    std::vector<size_t> parse_regions(std::istream&, size_t);

    /**
     * @brief Parses a methylation matrix from the input stream (columns: Cov Meth). 
     * 
     * @return matrix<int> 
     */
    Matrix<int> parse_methylation_matrix(std::istream&);

    /**
     * @brief Returns error message.
     * 
     * @return std::string 
     */
    std::string get_message();

    private:
    void get_shape(std::istream&, size_t&, size_t&);
    std::string message;
};

#endif