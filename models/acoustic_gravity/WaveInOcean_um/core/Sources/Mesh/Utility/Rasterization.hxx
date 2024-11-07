#include <algorithm>

namespace OndoMathX {


struct Rasterization
{
    std::vector<std::vector<Real>> grid;
    Real Ox;
    Real Oy;
    Index Nx;
    Index Ny;
    Real dx;
    Real dy;
};

inline Rasterization
CreateRasterization(const Mesh &mesh, Real Ox, Real Oy, Index Nx, Index Ny, Real dx, Real dy)
{
    Rasterization rasterization;
    
    rasterization.Ox = Ox;
    rasterization.Oy = Oy;
    rasterization.Nx = Nx;
    rasterization.Ny = Ny;
    rasterization.dx = dx;
    rasterization.dy = dy;
    
    
    
    rasterization.grid.resize(Nx);
    for (Index nx = 0; nx<Nx; nx++) rasterization.grid[nx].resize(Ny);
    
    Mesh planeMesh = mesh;
    
    planeMesh.makePlane();
    
    for (Index nx = 0; nx<Nx; nx++)
    {
        for (Index ny = 0; ny< Ny; ny++)
        {
            RealVector xy = { Ox+nx*dx, Oy+ny*dy, 0.0 };
            RealVector uv;
            
            bool pos_found = false;
            
            for (Index el = 0; el< planeMesh.getNumElements(); ++el)
            {
                const Element &       elem =      mesh.getElement(el);
                const Element & plane_elem = planeMesh.getElement(el);
                
                assert(elem.getType() == GeoTriangle3);
                
                plane_elem.xyz2uvw(xy, uv);
                
                RealVector uvw ={uv[0], uv[1], 0.0};
                
                if (RTriangle::isInside(uvw))
                {
                    plane_elem.xyz2uvw(xy, uv);
                    
                    Real z = (1-uv[0]-uv[1])*elem.getVertex(0)[2]
                    +          uv[0]*elem.getVertex(1)[2]
                    +          uv[1]*elem.getVertex(2)[2];
                    
                    
                    rasterization.grid[nx][ny] = z;
                    
                    pos_found = true;
                    break;
                }
            }
            
            assert(pos_found == true);
            
        }
    }
    
    return rasterization;
}




inline Mesh CreatMesh(Rasterization & rasterization)
{
    Index Nx = rasterization.Nx;
    Index Ny = rasterization.Ny;
    Real dx = rasterization.dx;
    Real dy = rasterization.dy;
    Real Ox = rasterization.Ox;
    Real Oy = rasterization.Oy;
    
    std::vector<std::shared_ptr<Point> > coords(Nx*Ny);
    
    std::vector<std::shared_ptr<Element> > elements;
    
    for (Index ny = 0;ny<Ny;ny++)
    for (Index nx = 0;nx<Nx;nx++)
    {
        coords[nx+Nx*ny] =  std::make_shared<Point>(Ox+nx*dx, Oy+ny*dy, rasterization.grid[nx][ny],nx+Nx*ny);
    }
    
    for (Index ny = 1;ny<Ny;ny++)
    for (Index nx = 1;nx<Nx;nx++)
    {
        std::vector<std::shared_ptr<Point> > local_nodes(4);
        
        local_nodes[0] = coords[(nx-1) + Nx*(ny-1)];
        local_nodes[1] = coords[(nx  ) + Nx*(ny-1)];
        local_nodes[2] = coords[(nx  ) + Nx*(ny )];
        local_nodes[3] = coords[(nx-1) + Nx*(ny )];
        
        elements.push_back(std::make_shared<Quadrangle4>(local_nodes,0));
    }
    
    return Mesh(coords, elements);
}


inline Real Interpolate(const Rasterization & rasterization, Real x, Real y)
{
    // Find the vertices of the quad containing the point
    Index nx = std::max((int)floor((x-rasterization.Ox)/rasterization.dx),0);
    Index ny = std::max((int)floor((y-rasterization.Oy)/rasterization.dy),0);
    
    if (nx == rasterization.Nx-1) nx--;
    if (ny == rasterization.Ny-1) ny--;
    
    Real u = (x - rasterization.Ox - nx * rasterization.dx)/rasterization.dx;
    Real v = (y - rasterization.Oy - ny * rasterization.dy)/rasterization.dy;
    
    // Find the z coordinate
    Real z1 = rasterization.grid[nx][ny];
    Real z2 = rasterization.grid[nx+1][ny];
    Real z3 = rasterization.grid[nx+1][ny+1];
    Real z4 = rasterization.grid[nx][ny+1];
    
    Real z = (1-u)*(1-v)*z1
             +   u*(1-v)*z2
             +   u*   v *z3
             +(1-u)*  v *z4;
    
    return z;
    
}

//Since the gradient is discontinuous in each element we need a direction (dx,dy) to define the element in which to compute the gradient
inline std::array<Real,2> InterpolateGradient(const Rasterization & rasterization, Real x, Real y, Real dx, Real dy)
{
    std::array<Real,2> grad;
    
    // Find the vertices of the quad containing the point
    Index nx = std::max((int)floor((x-rasterization.Ox)/ rasterization.dx),0);
    Index ny = std::max((int)floor((y-rasterization.Oy) / rasterization.dy),0);
    
    if (nx == rasterization.Nx) nx--;
    if (ny == rasterization.Ny) ny--;
    
    Real u = (x - rasterization.Ox - nx * rasterization.dx)/rasterization.dx;
    Real v = (y - rasterization.Oy - ny * rasterization.dy)/rasterization.dy;
    
    // Find the zb coordinate
    // Know in which cell the point is: necessary only on the boundary
    // If the point is on one of the boundaries:
    // If direction in the x coord is positive: the x-barycenter is smaller than the x-point hence the point is in the cell nx-1
    // If direction in the y coord is positive: the y-barycenter is smaller than the y-point hence the point is in the cell ny-1
    if (u < REF_COORD_TOL && dx>0) nx = nx-1;
    if (v < REF_COORD_TOL && dy>0) ny = ny-1;
        
    Real z1 = rasterization.grid[nx][ny];
    Real z2 = rasterization.grid[nx+1][ny];
    Real z3 = rasterization.grid[nx+1][ny+1];
    Real z4 = rasterization.grid[nx][ny+1];

    grad[0] = (1-v)*(z2-z1) + v*(z3-z4);
    grad[1] = (1-u)*(z4-z1) + u*(z3-z2);
    
    return grad;
}



} // namespace OndoMathX




