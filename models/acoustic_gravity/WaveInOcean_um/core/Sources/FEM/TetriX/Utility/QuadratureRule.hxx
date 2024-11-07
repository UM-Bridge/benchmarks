#pragma once


// -----------------------------------------------------------------------------------------//
namespace OndoMathX
{


QuadratureRule::QuadratureRule(std::string&& name,
                                QuadraturePoint::vector_const_shared_ptr&& point_list,
                                RefElement topology_id,
                                std::size_t degree_of_exactness)
: name_(std::move(name)), point_list_(std::move(point_list)), topology_id_(topology_id),
    degree_of_exactness_(degree_of_exactness)
{ }


QuadratureRule::QuadratureRule(std::string&& name, RefElement topology_id, std::size_t degree_of_exactness)
: name_(std::move(name)), topology_id_(topology_id), degree_of_exactness_(degree_of_exactness)
{ }


void QuadratureRule::AddQuadraturePoint(RealVector&& local_coords, double weight)
{
    auto&& quad_pt =
        std::make_shared<const QuadraturePoint>(std::move(local_coords), weight, GetName(), point_list_.size());

    point_list_.emplace_back(std::move(quad_pt));
}


inline Index QuadratureRule::NquadraturePoint() const noexcept
{
    assert(!point_list_.empty());
    return point_list_.size();
}


inline const QuadraturePoint& QuadratureRule::GetPoint(Index index) const noexcept
{
    assert(!point_list_.empty());
    assert(index < point_list_.size());
    assert(!(!point_list_[index]));
    return *(point_list_[index]);
}


inline Index QuadratureRule::DegreeOfExactness() const noexcept
{
    assert(degree_of_exactness_ != 0
            && "Should be called only if it is relevant for the law!");
    return degree_of_exactness_;
}


inline RefElement QuadratureRule::GetTopologyIdentifier() const noexcept
{
    return topology_id_;
}


inline const std::string& QuadratureRule::GetName() const noexcept
{
    assert(!name_.empty());
    return name_;
}


inline const QuadraturePoint::vector_const_shared_ptr& QuadratureRule::GetQuadraturePointList() const noexcept
{
    assert(!point_list_.empty());
    return point_list_;
}


} // namespace OndoMathX