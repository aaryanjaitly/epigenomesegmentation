#include <vector>
#include "boost/algorithm/string.hpp"

#include "matrix_reader.h"

void Reader::get_shape(std::istream& strm, size_t& rows,size_t& cols)
{
    std::string buffer;
    std::getline(strm, buffer);

    if (buffer.empty())
    {
        message = "Error: empty file.";
    }

    std::vector<std::string> row;
    boost::split(row, buffer, boost::is_any_of("\t "), boost::token_compress_on);
    cols = row.size();

    rows = 1;
    while(getline(strm, buffer))
    {
        ++rows;
    }
    strm.clear();
    strm.seekg(0, std::ios::beg);
}

Matrix<int> Reader::parse_matrix(std::istream& strm)
{
    size_t cols, rows;
    get_shape(strm, rows, cols);
    Matrix<int> m (rows, cols);

    std::string buffer;
    for (size_t i = 0; i < rows; ++i)
    {   
        std::getline(strm, buffer);
        std::vector<std::string> row;
        boost::split(row, buffer, boost::is_any_of("\t "), boost::token_compress_on);
        
        if (row.size() != cols)
        {
            message = "Error: not all rows have the same number of columns.";
            return m;
        }
        for (size_t j = 0; j < cols; ++j)
        {
            try
            {
                m(i, j) = std::stoi(row[j]);
            }
            catch (std::exception&)
            {
                message = "Error: matrix values cannot be converted to int.";
                return m;
            }
        }
    }

    return m;
}

std::vector<size_t> Reader::parse_regions(std::istream& strm, size_t max)
{
    std::istream_iterator<size_t> start(strm), end;
    std::vector<size_t> numbers(start, end);
    if (strm.eof())
    {
        if (numbers[0] != 0)
        {
            message = "Regions must start at 0.";
            return numbers;
        }

        if (!std::is_sorted(numbers.begin(), numbers.end()))
        {
            message = "Regions must be sorted.";
            return numbers;
        }

        if (max <= *numbers.rbegin())
        {
            message = "Last starting index must be smaller than the number of rows of the count matrix.";
            return numbers;
        }
    }
    else
    {
        message = "Wrong region format: cannot be converted to int (must be whitespace separated).";
    }
    return numbers;
}

Matrix<int> Reader::parse_methylation_matrix(std::istream& strm)
{
    size_t cols, rows;
    get_shape(strm, rows, cols);
    Matrix<int> m (rows, cols);

    if (cols != 2)
    {
        message = "Error: methylation matrix should have exactly two columns (coverage and methylation).";
        return m;
    }

    std::string buffer;
    for (size_t i = 0; i < rows; ++i)
    {   
        std::getline(strm, buffer);
        std::vector<std::string> row;
        boost::split(row, buffer, boost::is_any_of("\t "), boost::token_compress_on);
        
        if (row.size() != cols)
        {
            message = "Error: not all rows have the same number of columns.";
            return m;
        }
        try
        {
            m(i, 0) = std::stoi(row[0]);
            m(i, 1) = std::stoi(row[1]);
            if (m(i, 0) < m(i, 1))
            {
                message = "Error: coverage must be greater or equal to the methylation level.";
                return m;
            }
        }
        catch (std::exception&)
        {
            message = "Error: matrix values cannot be converted to int.";
            return m;
        }
    }
    return m;
}

std::string Reader::get_message()
{
    return message;
}