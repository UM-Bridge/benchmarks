#pragma once

// -----------------------------------------------------------------------------------------//
#include <vector>
#include <cmath>
#include <cassert>
#include <array>
#include <memory>
#include <cstddef> // IWYU pragma: keep
#include <iosfwd>

// -----------------------------------------------------------------------------------------//

#include "QuadraturePoint.h" 

// -----------------------------------------------------------------------------------------//
namespace OndoMathX
{


/*!
    * \brief Defines a quadrature rule.
    */
class QuadratureRule
{

    public:
    //! Alias to shared_ptr to constant object.
    using const_shared_ptr = std::shared_ptr<const QuadratureRule>;

    //! Alias to unique_ptr to constant object.
    using const_unique_ptr = std::unique_ptr<const QuadratureRule>;

    //! Alias to a vector of const_shared_ptr.
    using vector_const_shared_ptr = std::vector<const_shared_ptr>;


    public:
    /*!
        * \class doxygen_hide_quadrature_rule_constructor_args
        *
        * \param[in] topology_id Topology identifier.
        * \param[in] degree_of_exactness Degree of exactness of the rule. Choose NumericNS::UninitiazedIndex() for
        * rules for which it is pointless.
        *
        * \todo This should be handled more properly with inheritance, so that rules for which it is pointless do not
        * even have such an attribute, but convention above will do for the time being.
        */

    /*!
        * \class doxygen_hide_quadrature_rule_constructor_name_arg
        *
        * \param[in] name Name of the quadrature rule.
        */


    /// \name Constructors and destructor.
    ///@{

    /*!
        * Constructor.
        *
        * \copydetails doxygen_hide_quadrature_rule_constructor_args
        * \copydetails doxygen_hide_quadrature_rule_constructor_name_arg
        * \param[in] point_list List of quadrature points.
        */
    explicit QuadratureRule(std::string&& name,
                            QuadraturePoint::vector_const_shared_ptr&& point_list,
                            RefElement topology_id,
                            Index degree_of_exactness);

    /*!
        * Constructor.
        *
        * In this overload, quadrature points aren't yet defined and must be added afterwards through
        * \a AddQuadraturePoint() method.
        *
        * \copydetails doxygen_hide_quadrature_rule_constructor_args
        * \copydetails doxygen_hide_quadrature_rule_constructor_name_arg
        */
    explicit QuadratureRule(std::string&& name,
                            RefElement topology_id,
                            Index degree_of_exactness = 0);


    //! Destructor.
    ~QuadratureRule() = default;

    //! \copydoc doxygen_hide_copy_constructor
    QuadratureRule(const QuadratureRule& rhs) = delete;

    //! \copydoc doxygen_hide_move_constructor
    QuadratureRule(QuadratureRule&& rhs) = default;

    //! \copydoc doxygen_hide_copy_affectation
    QuadratureRule& operator=(const QuadratureRule& rhs) = delete;

    //! \copydoc doxygen_hide_move_affectation
    QuadratureRule& operator=(QuadratureRule&& rhs) = delete;


    ///@}


    /*!
        * \brief Add a new quadrature point.
        *
        * The quadrature point is created inside the function from \a local_coords and \a weight.
        *
        * \param[in] local_coords Local coordinates of the quadrature point to create.
        * \param[in] weight Weight of the quadrature point to create.
        */
    void AddQuadraturePoint(RealVector&& local_coords, double weight);

    //! Number of quadrature points.
    Index NquadraturePoint() const noexcept;

    //! Access to one quadrature point.
    //! \param[in] index Position in the vector of the sought \a QuadraturePoint.
    const QuadraturePoint& GetPoint(Index index) const noexcept;

    //! Return the degree of exactness if relevant (and will assert in debug if not).
    Index DegreeOfExactness() const noexcept;

    //! Identifier of the topology upon which the rule is defined.
    RefElement GetTopologyIdentifier() const noexcept;

    //! Returns the name of the quadrature rule.
    const std::string& GetName() const noexcept;

    //! Returns the list of quadrature points.
    const QuadraturePoint::vector_const_shared_ptr& GetQuadraturePointList() const noexcept;

    private:
    //! Name of the quadrature rule.
    std::string name_;

    //! List of quadrature points.
    QuadraturePoint::vector_const_shared_ptr point_list_;

    //! Identifier of the geometric element upon which the rule is defined.
    const RefElement topology_id_;

    /*!
        * \brief Degree of exactness.
        *
        * Equal to NumericNS::UninitiazedIndex() for rules for which it is pointless.
        *
        * \todo This should be handled more properly with inheritance, so that rules for which it is pointless do not
        * even have such an attribute, but convention above will do for the time being.
        */
    const Index degree_of_exactness_;
};



} // OndoMathX


// #include "Tet15.hxx"
#include "QuadratureRule.hxx" // IWYU pragma: export
#include "Instances/Tetrahedron.h"


