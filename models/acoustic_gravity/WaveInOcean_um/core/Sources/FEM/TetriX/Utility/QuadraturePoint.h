#pragma once 

#include <cstddef> // IWYU pragma: keep
#include <iosfwd>
#include <memory>
#include <vector>


namespace OndoMathX
{


/*!
    * \brief Quadrature point used to perform discrete integration.
    *
    */
class QuadraturePoint
{
    public:
    //! Alias to shared_ptr to constant object.
    using const_shared_ptr = std::shared_ptr<const QuadraturePoint>;

    //! Alias to a vector of const_shared_ptr.
    using vector_const_shared_ptr = std::vector<const_shared_ptr>;


    public:
    /// \name Constructors and destructor.
    ///@{

    /*!
        * \brief Constructor.
        *
        * \param[in] local_coords Coordinates of the quadrature point.
        * \param[in] weight Weight of the point in the integration.
        * \param[in] rule_name Name of the quadrature rule to which the point belongs.
        * \param[in] index Index within the quadrature rule.
        */
    explicit QuadraturePoint(RealVector&& local_coords,
                             double weight,
                             const std::string& rule_name,
                             Index index);

    //! Destructor.
    ~QuadraturePoint() = default;

    //! \copydoc doxygen_hide_copy_constructor
    QuadraturePoint(const QuadraturePoint& rhs) = delete;

    //! \copydoc doxygen_hide_move_constructor
    QuadraturePoint(QuadraturePoint&& rhs) = default;

    //! \copydoc doxygen_hide_copy_affectation
    QuadraturePoint& operator=(const QuadraturePoint& rhs) = delete;

    //! \copydoc doxygen_hide_move_affectation
    QuadraturePoint& operator=(QuadraturePoint&& rhs) = delete;

    ///@}

    //! Get the local coordinates of the quadrature point.
    const RealVector& GetCoordinates() const noexcept;

    //! Get the weight.
    double GetWeight() const noexcept;

    //! Print function (used also for operator<< overload).
    //! \copydoc doxygen_hide_stream_inout
    void Print(std::ostream& stream = std::cout) const;

    //! Get the name of the quadrature rule to which the point belongs.
    const std::string& GetRuleName() const;

    //! Get the index of the quadrature point within the quadrature rule.
    Index GetIndex() const noexcept;


    private:
    //! Coordinates of the quadrature point.
    const RealVector local_coords_;

    //! Weight of the point.
    const double weight_;

    //! Reference to the name of the quadrature rule to which the point belongs.
    const std::string& rule_name_;

    //! Index within the quadrature rule (from 0 to Nquadrature_point - 1).
    const Index index_;
};


//! \copydoc doxygen_hide_std_stream_out_overload
std::ostream& operator<<(std::ostream& stream, const QuadraturePoint& rhs);


} // namespace OndoMathX



#include "QuadraturePoint.hxx" // IWYU pragma: export


