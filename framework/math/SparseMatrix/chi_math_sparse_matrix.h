#pragma once

#include "framework/math/chi_math.h"

// ###################################################################
/**Sparse matrix utility. This is a basic CSR type sparse matrix
 * which allows efficient matrix storage and multiplication. It is
 * not intended for solving linear systems (use PETSc for that instead).
 * It was originally developed for the transfer matrices of transport
 * cross-sections.*/
class chi_math::SparseMatrix
{
private:
  size_t row_size_; ///< Maximum number of rows for this matrix
  size_t col_size_; ///< Maximum number of columns for this matrix

public:
  /**rowI_indices[i] is a vector indices j for the
   * non-zero columns.*/
  std::vector<std::vector<size_t>> rowI_indices_;
  /**rowI_values[i] corresponds to column indices and
   * contains the non-zero value.*/
  std::vector<std::vector<double>> rowI_values_;

public:
  SparseMatrix(size_t num_rows, size_t num_cols);
  SparseMatrix(const SparseMatrix& in_matrix);

  size_t NumRows() const { return row_size_; }
  size_t NumCols() const { return col_size_; }

  void Insert(size_t i, size_t j, double value);
  void InsertAdd(size_t i, size_t j, double value);
  double ValueIJ(size_t i, size_t j) const;
  void SetDiagonal(const std::vector<double>& diag);

  void Compress();

  std::string PrintStr() const;

private:
  void CheckInitialized() const;

public:
  virtual ~SparseMatrix() = default;

public:
  struct EntryReference
  {
    const size_t& row_index;
    const size_t& column_index;
    double& value;

    EntryReference(const size_t& row_id, const size_t& column_id, double& in_value)
      : row_index(row_id), column_index(column_id), value(in_value)
    {
    }
  };

  struct ConstEntryReference
  {
  public:
    const size_t& row_index;
    const size_t& column_index;
    const double& value;

    ConstEntryReference(const size_t& row_id, const size_t& column_id, const double& in_value)
      : row_index(row_id), column_index(column_id), value(in_value)
    {
    }
  };

  class RowIteratorContext
  {
  private:
    const std::vector<size_t>& ref_col_ids_;
    std::vector<double>& ref_col_vals_;
    const size_t ref_row_;

  public:
    RowIteratorContext(SparseMatrix& matrix, size_t ref_row)
      : ref_col_ids_(matrix.rowI_indices_[ref_row]),
        ref_col_vals_(matrix.rowI_values_[ref_row]),
        ref_row_(ref_row)
    {
    }

    class RowIterator
    {
    private:
      typedef RowIterator It;

    private:
      RowIteratorContext& context_;
      size_t ref_entry_;

    public:
      RowIterator(RowIteratorContext& context, size_t ref_entry)
        : context_{context}, ref_entry_{ref_entry}
      {
      }

      It operator++()
      {
        It i = *this;
        ref_entry_++;
        return i;
      }
      It operator++(int)
      {
        ref_entry_++;
        return *this;
      }

      EntryReference operator*()
      {
        return {
          context_.ref_row_, context_.ref_col_ids_[ref_entry_], context_.ref_col_vals_[ref_entry_]};
      }

      bool operator==(const It& rhs) const { return ref_entry_ == rhs.ref_entry_; }
      bool operator!=(const It& rhs) const { return ref_entry_ != rhs.ref_entry_; }
    };

    RowIterator begin() { return {*this, 0}; }
    RowIterator end() { return {*this, ref_col_vals_.size()}; }
  };

  RowIteratorContext Row(size_t row_id); // See .cc file

  class ConstRowIteratorContext
  {
  private:
    const std::vector<size_t>& ref_col_ids_;
    const std::vector<double>& ref_col_vals_;
    const size_t ref_row_;

  public:
    ConstRowIteratorContext(const SparseMatrix& matrix, size_t ref_row)
      : ref_col_ids_(matrix.rowI_indices_[ref_row]),
        ref_col_vals_(matrix.rowI_values_[ref_row]),
        ref_row_(ref_row)
    {
    }

    class ConstRowIterator
    {
    private:
      typedef ConstRowIterator It;

    private:
      const ConstRowIteratorContext& context_;
      size_t ref_entry_;

    public:
      ConstRowIterator(const ConstRowIteratorContext& context, size_t ref_entry)
        : context_(context), ref_entry_{ref_entry}
      {
      }

      It operator++()
      {
        It i = *this;
        ref_entry_++;
        return i;
      }
      It operator++(int)
      {
        ref_entry_++;
        return *this;
      }

      ConstEntryReference operator*()
      {
        return {
          context_.ref_row_, context_.ref_col_ids_[ref_entry_], context_.ref_col_vals_[ref_entry_]};
      }

      bool operator==(const It& rhs) const { return ref_entry_ == rhs.ref_entry_; }
      bool operator!=(const It& rhs) const { return ref_entry_ != rhs.ref_entry_; }
    };

    ConstRowIterator begin() const { return {*this, 0}; }
    ConstRowIterator end() const { return {*this, ref_col_vals_.size()}; }
  };

  ConstRowIteratorContext Row(size_t row_id) const;

  /**Iterator to loop over all matrix entries.*/
  class EntriesIterator
  {
  private:
    typedef EntriesIterator EIt;

  private:
    SparseMatrix& sp_matrix;
    size_t ref_row_;
    size_t ref_col_;

  public:
    explicit EntriesIterator(SparseMatrix& context, size_t row)
      : sp_matrix{context}, ref_row_{row}, ref_col_(0)
    {
    }

    void Advance()
    {
      ref_col_++;
      if (ref_col_ >= sp_matrix.rowI_indices_[ref_row_].size())
      {
        ref_row_++;
        ref_col_ = 0;
        while ((ref_row_ < sp_matrix.row_size_) and (sp_matrix.rowI_indices_[ref_row_].empty()))
          ref_row_++;
      }
    }

    EIt operator++()
    {
      EIt i = *this;
      Advance();
      return i;
    }
    EIt operator++(int)
    {
      Advance();
      return *this;
    }

    EntryReference operator*()
    {
      return {ref_row_,
              sp_matrix.rowI_indices_[ref_row_][ref_col_],
              sp_matrix.rowI_values_[ref_row_][ref_col_]};
    }
    bool operator==(const EIt& rhs) const
    {
      return (ref_row_ == rhs.ref_row_) and (ref_col_ == rhs.ref_col_);
    }
    bool operator!=(const EIt& rhs) const
    {
      return (ref_row_ != rhs.ref_row_) or (ref_col_ != rhs.ref_col_);
    }
  };

  EntriesIterator begin();
  EntriesIterator end();
};
