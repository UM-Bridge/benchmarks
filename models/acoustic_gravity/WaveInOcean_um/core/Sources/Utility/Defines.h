#pragma once



 
namespace OndoMathX {
            
        // IO file formats
        enum MeshFormat
        {
            FormatMsh,
            FormatMesh,
            FormatVtk
        };
        
        // Element enum (following msh standard)
        enum GeoElement
        {
            GeoNone         =0,
            GeoSegment2     =1,
            GeoTriangle3    =2,
            GeoQuadrangle4  =3,
            GeoTetrahedron4 =4,
            GeoHexahedron8  =5,
            //GeoPrism6       =6,
            //GeoPyramid5     =7,
            GeoSegment3     =8,
            GeoTriangle6    =9,
            GeoQuadrangle9  =10,
            GeoHexahedron27 =12,
            GeoPoint        =15 
        };
        
        //Enum for the reference elements
        enum RefElement
        {
            RefNone,
            RefPoint,
            Segment,
            Triangle,
            Quadrangle,
            Tetrahedron,
            Hexahedron,
            Prism,
            Pyramid
        };
        
        //Tolerance on nodes position (on the reference element and in absolute)
        constexpr Real REF_COORD_TOL = 1e-12;
        constexpr Real COORD_TOL = 1e-12;
        
} //OndoMathX

 
