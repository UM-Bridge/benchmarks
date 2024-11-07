
namespace OndoMathX {
    
inline void  getBoundingBox(const Mesh & mesh,
                     Real &x_min,Real &x_max,
                     Real &y_min,Real &y_max,
                     Real &z_min,Real &z_max, bool useLabel = false, Index label = 0) 
{
    x_min= std::numeric_limits<Real>::max();
    x_max=-std::numeric_limits<Real>::max();
    y_min= std::numeric_limits<Real>::max();
    y_max=-std::numeric_limits<Real>::max();
    z_min= std::numeric_limits<Real>::max();
    z_max=-std::numeric_limits<Real>::max();
    
    bool atLeastOne = false;
    
    for (Index el=0;el<mesh.getNumElements();++el)
    {
        const Element & elem = mesh.getElement(el);
        
        if( (useLabel && elem.getLabel() == label) || !useLabel)
        {
            for (Index k=0;k<elem.getNumNodes();++k)
            {
                const Point &xyz = elem.getNode(k);
                
                if (xyz[0]<x_min) x_min=xyz[0];
                if (xyz[0]>x_max) x_max=xyz[0];
                if (xyz[1]<y_min) y_min=xyz[1];
                if (xyz[1]>y_max) y_max=xyz[1];
                if (xyz[2]<z_min) z_min=xyz[2];
                if (xyz[2]>z_max) z_max=xyz[2];
                
                atLeastOne = true;
            }
        }
    }
    
    if (atLeastOne == false)
        x_min=x_max=y_min=y_max=z_min=z_max=0.0;
}

     
    
} // namespace OndoMathX




