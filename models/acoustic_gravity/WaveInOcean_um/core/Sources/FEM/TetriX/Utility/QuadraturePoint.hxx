#pragma once

#include <cstddef> // IWYU pragma: keep
#include <iosfwd>

namespace OndoMathX
{


QuadraturePoint::QuadraturePoint(RealVector&& local_coords,
                                    double weight,
                                    const std::string& rule_name,
                                    const std::size_t index)
: local_coords_(std::move(local_coords)), weight_(weight), rule_name_(rule_name), index_(index)
{ }


void QuadraturePoint::Print(std::ostream& out) const
{
    out << "\nQuadrature point #"<< GetIndex() << " at local coords: [";
    std::string separator("");
    for (auto value : GetCoordinates())
    {
        out << separator << value;
        separator = ", ";
    }
    out << "], (weight = " << GetWeight() << ')';
}


std::ostream& operator<<(std::ostream& stream, const QuadraturePoint& point)
{
    point.Print(stream);
    return stream;
}

inline const RealVector& QuadraturePoint::GetCoordinates() const noexcept
{
    return local_coords_;
}


inline double QuadraturePoint::GetWeight() const noexcept
{
    return weight_;
}


inline const std::string& QuadraturePoint::GetRuleName() const
{
    return rule_name_;
}


inline Index QuadraturePoint::GetIndex() const noexcept
{
    return index_;
}


} // namespace OndoMathX



