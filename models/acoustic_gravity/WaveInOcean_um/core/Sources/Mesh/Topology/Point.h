#pragma once

 

# include <array>
# include <cassert>
# include <iostream>
# include <math.h>
# include <vector>
#include "../../Utility/Defines.h"

namespace OndoMathX {
    
 
        
        /*!
         * \brief Define a point in an euclidian frame.
         */
        class Point
        {
            
        private:
            
            //! Coordinates.
            RealVector _coords;
            Index _index;
            
        public:
            
            
            //! \name Constructors and destructor.
            //@{
            
            // Default construcotr, set everything to zero
            explicit Point() {zero();};
            
            //! From three scalars and an integer.
            explicit Point(Real x, Real y = 0.0, Real z = 0.0, Index index = 0)
                : _coords({{x, y, z}}), _index(index) { };
            
            explicit Point(RealVector p, Index index = 0)
                : _coords(p), _index(index) { };
            
            //! Copy constructor.
            Point(const Point& rhs) = default;
            
            //! Move constructor.
            explicit Point(Point&& rhs) = default;
            
            //! Destructor.
            virtual ~Point();
            
            //@}
            
            //! Assignment
            Point& operator=(const Point& rhs) = default;
            
            
            //! Return the first component of the point.
            Real x() const {return _coords[0];};
            
            //! Return the second component of the point (if relevant; if not an exception is thrown).
            Real y() const {return _coords[1];};
            
            //! Return the third component of the point (if relevant; if not an exception is thrown).
            Real z() const {return _coords[2];};

                        
            //! Return the first component of the point.
            Real& x() {return _coords[0];};
            
            //! Return the second component of the point (if relevant; if not an exception is thrown).
            Real& y() {return _coords[1];};
            
            //! Return the third component of the point (if relevant; if not an exception is thrown).
            Real& z() {return _coords[2];};
            
            
            //! Return the index
            Index getIndex() const {return _index;}
            
            //! Set the index
            void setIndex(Index index)  {_index=index;}
            
            const RealVector & coord() const {return _coords;}
            
            /*!
             * \brief Subscript operator, non-const version.
             *
             * \param[in] i Component i. x, y and z coordinates may be accessed respectively with [0], [1], [2].
             *
             * \return Coordinate at the index-th component.
             */
            inline Real& operator[](Index i) {assert(i < 3); return _coords[i];};
            
            /*!
             * \brief Subscript operator, const version.
             *
             * \param[in] i Component i. x, y and z coordinates may be accessed respectively with [0], [1], [2].
             *
             * \return Coordinate at the index-th component.
             */
            inline const Real& operator[](Index i) const {assert(i < 3); return _coords[i];};
            
            //! Set all coordinates and index at 0.
            inline void zero() {std::fill(_coords.begin(), _coords.end(), 0.);_index = 0;};
            
            //! Print function.
            void print(std::ostream& out) const;
            
            
        };
        
        
        /*!
         * \brief Calculates the distance between two points, using an euclidian norm.
         */
        Real distance(const Point&, const Point&);
        
        
        Real norm(const Point&);
        
        Point getBarycenter(const std::vector<Point>&);
        
        struct sameCoordinatesXYZ
        {
            bool operator()(const Point &n1, const Point &n2) const
            {
                if( (fabs(n1[0] - n2[0]) < REF_COORD_TOL)
                   && (fabs(n1[1] - n2[1]) < REF_COORD_TOL)
                   && (fabs(n1[2] - n2[2]) < REF_COORD_TOL)) return true;
                
                return false;
            }
        };
        
        struct lessCoordinatesXYZ
        {
            
            bool operator()(const Point &n1, const Point &n2) const
            {
                if(n1[2] < n2[2] - REF_COORD_TOL) return true;
                if(n2[2] < n1[2] - REF_COORD_TOL) return false;
                
                if(n1[1] < n2[1] - REF_COORD_TOL) return true;
                if(n2[1] < n1[1] - REF_COORD_TOL) return false;
                
                if(n1[0] < n2[0] - REF_COORD_TOL) return true;
                if(n2[0] < n1[0] - REF_COORD_TOL) return false;
                
                
                return false;
            }
        };
    
    
        inline bool operator<(const Point &n1, const Point &n2)
        {
            lessCoordinatesXYZ cmp;
            return cmp(n1,n2);
        }

        
        struct lessCoordinatesAngle2D
        {
            Point Barycenter;
            
            bool operator()(const Point &n1, const Point &n2) const
            {
                Point d1(n1[0]-Barycenter[0],n1[1]-Barycenter[1],0);
                Point d2(n2[0]-Barycenter[0],n2[1]-Barycenter[1],0);
                
                Real nd1 = norm(d1);
                Real nd2 = norm(d2);
                
                if (nd1 < REF_COORD_TOL)
                {
                    if (nd2 < REF_COORD_TOL) return false;
                    else return true;
                }
                if (nd2 < REF_COORD_TOL) return false;
                
                //Renormalization
                d1[0]/=nd1; d1[1]/=nd1;
                d2[0]/=nd2; d2[1]/=nd2;
                
                //Compute the angles
                Real theta1 = acos(d1[0]);
                if (d1[1]<0) theta1=2*M_PI - theta1;
                
                Real theta2 = acos(d2[0]);
                if (d2[1]<0) theta2=2*M_PI - theta2;
                
                //Test the angles
                if (theta1<theta2) return true;
                
                return false;
            }
        };


inline Point::~Point() = default;
 

inline Real distance(const Point& Point1, const Point& Point2)
{
    
    Real ret = 0.0;
    
    for (Index i = 0u; i < 3; ++i)
        ret += (Point1[i] - Point2[i])*(Point1[i] - Point2[i]);
    
    return sqrt(ret);
}

inline Real norm(const Point& Point1)
{
    
    Real ret = 0.0;
    
    for (Index i = 0u; i < 3; ++i)
        ret +=  Point1[i]*Point1[i];
    
    return sqrt(ret);
}

inline Point getBarycenter(const std::vector<Point>& list)
{
    Point barycenter;
    
    for (auto coord:list)
    {
        for (Index i=0;i<3;++i)
            barycenter[i]+=coord[i];
    }
    
    for (Index i=0;i<3;++i)
        barycenter[i]/=list.size();
    
    return barycenter;
}

inline void Point::print(std::ostream& stream) const
{
    stream << '(';
    
    for (Index i = 0u; i < 2; ++i)
        stream << _coords[i] << ", ";
    
    stream << _coords[2] << "; " << _index << ')';
}
     
} //OndoMathX

namespace std
{
    inline std::ostream& operator<<(std::ostream& stream, const OndoMathX::Point& Point)
    {
        Point.print(stream);
        return stream;
    }
}
 
