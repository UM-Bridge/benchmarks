#pragma once

#include<memory>

#include "FEM/TetriX/Utility/QuadratureRule.h"


namespace OndoMathX::QuadratureNS
{

enum Type
{
    Tet5 = 0,
    Tet15 = 1
};


static std::array<QuadratureRule::const_shared_ptr, 2> GetTetrahedronQuadratureRules()
{
    auto five_points_ptr =
        std::make_shared<QuadratureRule>("tetrahedron_5_points", OndoMathX::RefElement::Tetrahedron, 2);
    {
        auto& five_points = *five_points_ptr;

        static constexpr Real weight_vertex = 1./ 20.;
        five_points.AddQuadraturePoint(RealVector{0., 0., 0.}, weight_vertex);
        five_points.AddQuadraturePoint(RealVector{1., 0., 0.}, weight_vertex);
        five_points.AddQuadraturePoint(RealVector{0., 1., 0.}, weight_vertex);
        five_points.AddQuadraturePoint(RealVector{0., 0., 1.}, weight_vertex);
        
        static constexpr Real weight_barycenter = 4. / 5.;
        five_points.AddQuadraturePoint(RealVector{0.25, 0.25, 0.25}, weight_barycenter);
    }

    auto fifteen_points_ptr =
        std::make_shared<QuadratureRule>("tetrahedron_15_points", OndoMathX::RefElement::Tetrahedron, 3);

    {
        auto& fifteen_points = *fifteen_points_ptr;

        static constexpr Real weight_vertex = 17. / 5040.;
        fifteen_points.AddQuadraturePoint(RealVector{0., 0., 0.}, weight_vertex);
        fifteen_points.AddQuadraturePoint(RealVector{1., 0., 0.}, weight_vertex);
        fifteen_points.AddQuadraturePoint(RealVector{0., 1., 0.}, weight_vertex);
        fifteen_points.AddQuadraturePoint(RealVector{0., 0., 1.}, weight_vertex);

        static constexpr Real weight_edge = 2. / 315.;
        fifteen_points.AddQuadraturePoint(RealVector{0.5, 0.0, 0.0},weight_edge);
        fifteen_points.AddQuadraturePoint(RealVector{0.5, 0.5, 0.0},weight_edge);
        fifteen_points.AddQuadraturePoint(RealVector{0.0, 0.5, 0.0},weight_edge);
        fifteen_points.AddQuadraturePoint(RealVector{0.0, 0.0, 0.5},weight_edge);
        fifteen_points.AddQuadraturePoint(RealVector{0.5, 0.0, 0.5},weight_edge);
        fifteen_points.AddQuadraturePoint(RealVector{0.0, 0.5, 0.5},weight_edge);

        static constexpr Real one_third = 1. / 3.;
        static constexpr Real weight_face = 9. / 560.;
        fifteen_points.AddQuadraturePoint(RealVector{one_third, one_third, 0.}, weight_face);
        fifteen_points.AddQuadraturePoint(RealVector{one_third, 0., one_third}, weight_face);
        fifteen_points.AddQuadraturePoint(RealVector{0., one_third, one_third}, weight_face);
        fifteen_points.AddQuadraturePoint(RealVector{one_third, one_third, one_third}, weight_face);

        static constexpr Real weight_barycenter = 16. / 315.;
        fifteen_points.AddQuadraturePoint(RealVector{0.25, 0.25, 0.25}, weight_barycenter);
    }
    
    assert(five_points_ptr->NquadraturePoint() == 5);
    assert(fifteen_points_ptr->NquadraturePoint() == 15);
    return { { five_points_ptr, fifteen_points_ptr } };

}

} // namespace OndoMathX::QuadratureNS
