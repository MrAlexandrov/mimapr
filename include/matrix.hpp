#pragma once

#include <initializer_list>
#include <stdexcept>
#include <vector>
#include <iostream>

namespace NMatrix {

template <typename T = double>
class TMatrix {
private:
    struct ProxyRow {
        using iterator = typename std::vector<T>::iterator;
        using const_iterator = typename std::vector<T>::const_iterator;
        using reverse_iterator = typename std::vector<T>::reverse_iterator;
        using const_reverse_iterator = typename std::vector<T>::const_reverse_iterator;

        constexpr ProxyRow(size_t cols, const T& value = T()) : proxy_data_(cols, value) {}
        constexpr ProxyRow(const ProxyRow&) = default;
        constexpr ProxyRow(ProxyRow&&) = default;
        constexpr ProxyRow& operator=(const ProxyRow&) = default;
        constexpr ProxyRow& operator=(ProxyRow&&) = default;
        constexpr bool operator==(const ProxyRow& other) const;
        

        constexpr ProxyRow(const std::vector<T>& other) : proxy_data_(other) {}
        constexpr ProxyRow(std::vector<T>&& other) : proxy_data_(std::move(other)) {}
        constexpr ProxyRow& operator=(const std::vector<T>& other);
        constexpr ProxyRow& operator=(std::vector<T>&& other);
        constexpr bool operator==(const std::vector<T>& other) const;
        
        ProxyRow(std::initializer_list<T> list) : proxy_data_(list) {}
        
        operator std::vector<T>() const; 

        friend bool operator==(const std::vector<T>& vec, const ProxyRow& row) {
            return vec == row.proxy_data_;
        }

        constexpr size_t size() const noexcept;

        constexpr iterator begin() noexcept { return proxy_data_.begin(); }
        constexpr const_iterator begin() const noexcept { return proxy_data_.begin(); }
        constexpr iterator end() noexcept { return proxy_data_.end(); }
        constexpr const_iterator end() const noexcept { return proxy_data_.end(); }

        constexpr const_iterator cbegin() const noexcept { return proxy_data_.cbegin(); }
        constexpr const_iterator cend() const noexcept { return proxy_data_.cend(); }

        constexpr reverse_iterator rbegin() noexcept { return proxy_data_.rbegin(); }
        constexpr const_reverse_iterator rbegin() const noexcept { return proxy_data_.rbegin(); }
        constexpr reverse_iterator rend() noexcept { return proxy_data_.rend(); }
        constexpr const_reverse_iterator rend() const noexcept { return proxy_data_.rend(); }

        constexpr const_reverse_iterator crbegin() const noexcept { return proxy_data_.crbegin(); }
        constexpr const_reverse_iterator crend() const noexcept { return proxy_data_.crend(); }

        T& operator[](size_t col);
        const T& operator[](size_t col) const;

        std::vector<T> proxy_data_;
    };
public:
    using value_type = T;
    TMatrix() = default;
    constexpr TMatrix(size_t rows, size_t cols, const T& value = T()) : data_(rows, ProxyRow(cols, value)) {}
    constexpr TMatrix(const TMatrix&) = default;
    constexpr TMatrix(TMatrix&&) = default;
    constexpr TMatrix& operator=(const TMatrix&) = default;
    constexpr TMatrix& operator=(TMatrix&&) = default;

    constexpr TMatrix(const std::vector<ProxyRow>& data) : data_(data) {}
    constexpr TMatrix(std::initializer_list<std::initializer_list<T>> list);
    
    const size_t rows() const noexcept;
    const size_t cols() const noexcept;

	ProxyRow& operator[](size_t row);
	const ProxyRow& operator[](size_t row) const;

	TMatrix& operator+=(const TMatrix&);
	TMatrix& operator-=(const TMatrix&);
	TMatrix& operator*=(const TMatrix&);
	// TMatrix& operator/=(const TMatrix&);

    TMatrix operator+(const TMatrix&) const;
    TMatrix operator-(const TMatrix&) const;
    TMatrix operator*(const TMatrix&) const;

    TMatrix& operator+=(const T&);
    TMatrix& operator-=(const T&);
    TMatrix& operator*=(const T&);
    TMatrix& operator/=(const T&);

    TMatrix operator+(const T&) const;
    TMatrix operator-(const T&) const;
    TMatrix operator*(const T&) const;
    TMatrix operator/(const T&) const;
	
    TMatrix transpose() const;

    bool operator==(const TMatrix&) const;
    // TODO: add concepts
    template <typename U>
    friend std::ostream& operator<<(std::ostream&, const TMatrix<U>&);

private:
    // TODO: use pimpl
    std::vector<ProxyRow> data_;
};

template <typename T>
TMatrix<T>::ProxyRow::operator std::vector<T>() const {
    return proxy_data_;
} 

template <typename T>
constexpr typename TMatrix<T>::ProxyRow& TMatrix<T>::ProxyRow::operator=(const std::vector<T>& other) {
    proxy_data_ = other;
    return *this;
}

template <typename T>
constexpr typename TMatrix<T>::ProxyRow& TMatrix<T>::ProxyRow::operator=(std::vector<T>&& other) {
    proxy_data_ = std::move(other);
    return *this;
}

template <typename T>
constexpr bool TMatrix<T>::ProxyRow::operator==(const ProxyRow& other) const {
    return proxy_data_ == other.proxy_data_;
}

template <typename T>
constexpr bool TMatrix<T>::ProxyRow::operator==(const std::vector<T>& other) const {
    return proxy_data_ == other;
}

template <typename T>
constexpr size_t TMatrix<T>::ProxyRow::size() const noexcept {
    return proxy_data_.size();
} 

template <typename T>
T& TMatrix<T>::ProxyRow::operator[](size_t col) {
    if (col >= proxy_data_.size()) {
        throw std::out_of_range("Index out of range");
    }
    return proxy_data_[col];
}
template <typename T>
const T& TMatrix<T>::ProxyRow::operator[](size_t col) const {
    if (col >= proxy_data_.size()) {
        throw std::out_of_range("Index out of range");
    }
    return proxy_data_[col];
}


template <typename T>
constexpr TMatrix<T>::TMatrix(std::initializer_list<std::initializer_list<T>> list) {
    size_t rows_ = list.size();
    size_t cols_ = list.begin()->size();

    data_.reserve(rows_);
    for (const auto& row : list) {
        if (row.size() != cols_) {
            throw std::invalid_argument("All rows must have the same number of columns");
        }
        data_.emplace_back(row);
    }
}

template <typename T>
const size_t TMatrix<T>::rows() const noexcept { return data_.size(); }


template <typename T>
const size_t TMatrix<T>::cols() const noexcept { return data_.empty() ? 0 : data_.front().size(); }


template <typename T>
typename TMatrix<T>::ProxyRow& TMatrix<T>::operator[](size_t row) {
    if (row >= rows()) {
        throw std::out_of_range("Index out of range");
    }
    return data_[row];
}

template <typename T>
const typename TMatrix<T>::ProxyRow& TMatrix<T>::operator[](size_t row) const {
    if (row >= rows()) {
        throw std::out_of_range("Index out of range");
    }
    return data_[row];
}

template <typename T>
TMatrix<T>& TMatrix<T>::operator+=(const TMatrix<T>& other) {
    if (rows() != other.rows() || cols() != other.cols()) {
        throw std::invalid_argument("Matrix dimensions must match for addition");
    }
    for (size_t i = 0; i < rows(); ++i) {
        for (size_t j = 0; j < cols(); ++j) {
            (*this)[i][j] += other[i][j];
        }
    }
    return *this;
}

template <typename T>
TMatrix<T>& TMatrix<T>::operator-=(const TMatrix<T>& other) {
    if (rows() != other.rows() || cols() != other.cols()) {
        throw std::invalid_argument("Matrix dimensions must match for subtraction");
    }
    for (size_t i = 0; i < rows(); ++i) {
        for (size_t j = 0; j < cols(); ++j) {
            (*this)[i][j] -= other[i][j];
        }
    }
    return *this;
}

template <typename T>
TMatrix<T>& TMatrix<T>::operator*=(const TMatrix<T>& other) {
    if (cols() != other.rows()) {
        throw std::invalid_argument("Matrix dimensions are not suitable for multiplication");
    }

    TMatrix<T> result(rows(), other.cols(), T{});

    for (size_t i = 0; i < rows(); ++i) {
        for (size_t j = 0; j < other.cols(); ++j) {
            for (size_t k = 0; k < cols(); ++k) {
                result[i][j] += (*this)[i][k] * other[k][j];
            }
        }
    }

    *this = std::move(result);
    return *this;
}


template <typename T>
TMatrix<T> TMatrix<T>::operator+(const TMatrix<T>& other) const {
	TMatrix<T> result(*this);
	return result += other;
}

template <typename T>
TMatrix<T> TMatrix<T>::operator-(const TMatrix<T>& other) const {
	TMatrix<T> result(*this);
	return result -= other;
}

template <typename T>
TMatrix<T> TMatrix<T>::operator*(const TMatrix<T>& other) const {
	TMatrix<T> result(*this);
	return result *= other;
}


template <typename T>
TMatrix<T>& TMatrix<T>::operator+=(const T& scalar) {
    for (auto& rows : data_) {
        for (auto& element : rows) {
            element += scalar;
        }
    }
    return *this;
}

template <typename T>
TMatrix<T>& TMatrix<T>::operator-=(const T& scalar) {
    for (auto& rows : data_) {
        for (auto& element : rows) {
            element -= scalar;
        }
    }
    return *this;
}

template <typename T>
TMatrix<T>& TMatrix<T>::operator*=(const T& scalar) {
    for (auto& rows : data_) {
        for (auto& element : rows) {
            element *= scalar;
        }
    }
    return *this;
}

template <typename T>
TMatrix<T>& TMatrix<T>::operator/=(const T& scalar) {
    T zero = T(0.0);
    if (scalar == zero) {
        throw std::runtime_error("Division by zero!");
    }
    for (auto& rows : data_) {
        for (auto& element : rows) {
            element /= scalar;
        }
    }
    return *this;
}

template <typename T>
TMatrix<T> TMatrix<T>::operator+(const T& scalar) const {
    TMatrix<T> result(*this);
    return result += scalar;
}

template <typename T>
TMatrix<T> TMatrix<T>::operator-(const T& scalar) const {
    TMatrix<T> result(*this);
    return result -= scalar;
}

template <typename T>
TMatrix<T> TMatrix<T>::operator*(const T& scalar) const {
    TMatrix<T> result(*this);
    return result *= scalar;
}

template <typename T>
TMatrix<T> TMatrix<T>::operator/(const T& scalar) const {
    TMatrix<T> result(*this);
    return result /= scalar;
}

template <typename T>
TMatrix<T> TMatrix<T>::transpose() const {
	TMatrix<T> result(cols(), rows());
	for (size_t i = 0; i < rows(); ++i) {
		for (size_t j = 0; j < cols(); ++j) {
			result[j][i] = (*this)[i][j];
		}
	}
	return result;
}

template <typename T>
bool TMatrix<T>::operator==(const TMatrix<T>& other) const {
	return data_ == other.data_;
}

template <typename T>
std::ostream& operator<<(std::ostream& out, const TMatrix<T>& matrix) {
	for (const auto& row : matrix.data_) {
		for (const auto& elem : row) {
			out << elem << ' ';
		}
		out << '\n';
	}
	return out;
}

} // namespace NMatrix