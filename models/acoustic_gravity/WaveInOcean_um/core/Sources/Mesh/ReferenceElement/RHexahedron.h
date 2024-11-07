#pragma once

#include "RElementT.h"
#include "RQuadrangle.h"

namespace OndoMathX
{

    class RHexahedron : public RElementT<RHexahedron>
    {
    private:
        static inline void _pos2PosRefFace(const RealVector &pos, RealVector &posRef, Index numFace)
        {
            assert(numFace < 6);
            assert(isInside(pos));

            if (numFace == 0)
            {
                posRef[0] = pos[0];
                posRef[1] = pos[1];
            }
            if (numFace == 1)
            {
                posRef[0] = pos[0];
                posRef[1] = pos[2];
            }
            if (numFace == 2)
            {
                posRef[0] = pos[1];
                posRef[1] = pos[2];
            }
            if (numFace == 3)
            {
                posRef[0] = pos[1];
                posRef[1] = pos[2];
            }
            if (numFace == 4)
            {
                posRef[0] = pos[0];
                posRef[1] = pos[2];
            }
            if (numFace == 5)
            {
                posRef[0] = pos[0];
                posRef[1] = pos[1];
            }

            posRef[2] = 0;
        }

        static inline void _posRefFace2Pos(const RealVector &posRef, RealVector &pos, Index numFace)
        {
            assert(numFace < 6);
            assert(isInside(posRef));

            switch (numFace)
            {
            case 0:
                pos[0] = posRef[0];
                pos[1] = posRef[1];
                pos[2] = 0.0;
                break;
            case 1:
                pos[0] = posRef[0];
                pos[1] = 0.0;
                pos[2] = posRef[1];
                break;
            case 2:
                pos[0] = 0.0;
                pos[1] = posRef[0];
                pos[2] = posRef[1];
                break;
            case 3:
                pos[0] = 1.0;
                pos[1] = posRef[0];
                pos[2] = posRef[1];
                break;
            case 4:
                pos[0] = posRef[0];
                pos[1] = 1.0;
                pos[2] = posRef[1];
                break;
            case 5:
                pos[0] = posRef[0];
                pos[1] = posRef[1];
                pos[2] = 1.0;
                break;
            }
        }

    public:
        static inline Index getNumFaces() { return 6; }
        static inline Index getNumEdges() { return 12; }
        static inline Index getNumVertices() { return 8; }
        static inline Index getDim() { return 3; }
        static inline RefElement getRefType() { return Hexahedron; }
        static inline Index getNumInterfaces() { return 6; }

        static inline RefElement getInterfaceRefType(Index numInterface)
        {
            assert(numInterface < 6);
            return Quadrangle;
        }

        static inline RefElement getFaceRefType(Index numInterface)
        {
            assert(numInterface < 6);
            return Quadrangle;
        }

        static inline void getNormal(Index numInterface, RealVector &normal)
        {
            assert((numInterface >= 0) && (numInterface < 6));

            switch (numInterface)
            {
            case 0:
                normal[0] = 0.0;
                normal[1] = 0.0;
                normal[2] = -1.0;
                break;
            case 1:
                normal[0] = 0.0;
                normal[1] = -1.0;
                normal[2] = 0.0;
                break;
            case 2:
                normal[0] = -1.0;
                normal[1] = 0.0;
                normal[2] = 0.0;
                break;
            case 3:
                normal[0] = 1.0;
                normal[1] = 0.0;
                normal[2] = 0.0;
                break;
            case 4:
                normal[0] = 0.0;
                normal[1] = 1.0;
                normal[2] = 0.0;
                break;
            case 5:
                normal[0] = 0.0;
                normal[1] = 0.0;
                normal[2] = 1.0;
                break;
            default:
                break;
            }
        }

        static inline void getPosVertices(RealVector &pos, Index numRealVector)
        {
            assert(numRealVector < 8);

            switch (numRealVector)
            {
                {
                case 0:
                    pos[0] = 0.0;
                    pos[1] = 0.0;
                    pos[2] = 0.0;
                    break;
                case 1:
                    pos[0] = 1.0;
                    pos[1] = 0.0;
                    pos[2] = 0.0;
                    break;
                case 2:
                    pos[0] = 1.0;
                    pos[1] = 1.0;
                    pos[2] = 0.0;
                    break;
                case 3:
                    pos[0] = 0.0;
                    pos[1] = 1.0;
                    pos[2] = 0.0;
                    break;
                case 4:
                    pos[0] = 0.0;
                    pos[1] = 0.0;
                    pos[2] = 1.0;
                    break;
                case 5:
                    pos[0] = 1.0;
                    pos[1] = 0.0;
                    pos[2] = 1.0;
                    break;
                case 6:
                    pos[0] = 1.0;
                    pos[1] = 1.0;
                    pos[2] = 1.0;
                    break;
                case 7:
                    pos[0] = 0.0;
                    pos[1] = 1.0;
                    pos[2] = 1.0;
                    break;
                }
            }
        }

        static inline bool isInside(const RealVector &pos)
        {
            if (pos[0] < 1.0 + REF_COORD_TOL && pos[0] > -REF_COORD_TOL && pos[1] < 1.0 + REF_COORD_TOL && pos[1] > -REF_COORD_TOL && pos[2] < 1.0 + REF_COORD_TOL && pos[2] > -REF_COORD_TOL)
                return true;

            return false;
        }

        static inline bool isOnEdge(const RealVector &pos, Index numInterface)
        {
            assert(numInterface < 12);
            if (!isInside(pos))
                return false;

            switch (numInterface)
            {
            case 0:
                if (pos[1] < REF_COORD_TOL && pos[2] < REF_COORD_TOL)
                    return true;
                break;
            case 1:
                if (pos[0] > 1.0 - REF_COORD_TOL && pos[2] < REF_COORD_TOL)
                    return true;
                break;
            case 2:
                if (pos[1] > 1.0 - REF_COORD_TOL && pos[2] < REF_COORD_TOL)
                    return true;
                break;
            case 3:
                if (pos[0] < REF_COORD_TOL && pos[2] < REF_COORD_TOL)
                    return true;
                break;

            case 4:
                if (pos[0] < REF_COORD_TOL && pos[1] < REF_COORD_TOL)
                    return true;
                break;
            case 5:
                if (pos[0] > 1.0 - REF_COORD_TOL && pos[1] < REF_COORD_TOL)
                    return true;
                break;
            case 6:
                if (pos[0] > 1.0 - REF_COORD_TOL && pos[1] > 1.0 - REF_COORD_TOL)
                    return true;
                break;
            case 7:
                if (pos[0] < REF_COORD_TOL && pos[1] > 1.0 - REF_COORD_TOL)
                    return true;
                break;

            case 8:
                if (pos[1] < REF_COORD_TOL && pos[2] > 1.0 - REF_COORD_TOL)
                    return true;
                break;
            case 9:
                if (pos[0] > 1.0 - REF_COORD_TOL && pos[2] > 1.0 - REF_COORD_TOL)
                    return true;
                break;
            case 10:
                if (pos[1] > 1.0 - REF_COORD_TOL && pos[2] > 1.0 - REF_COORD_TOL)
                    return true;
                break;
            case 11:
                if (pos[0] < REF_COORD_TOL && pos[2] > 1.0 - REF_COORD_TOL)
                    return true;
                break;
            }
            return false;
        }

        static inline bool isOnVertex(const RealVector &pos, Index numInterface)
        {
            assert(numInterface < 8);
            switch (numInterface)
            {
            case 0:
                if (pos[0] < REF_COORD_TOL && pos[1] < REF_COORD_TOL && pos[2] < REF_COORD_TOL)
                    return true;
                break;
            case 1:
                if (pos[0] > 1.0 - REF_COORD_TOL && pos[1] < REF_COORD_TOL && pos[2] < REF_COORD_TOL)
                    return true;
                break;
            case 2:
                if (pos[0] > 1.0 - REF_COORD_TOL && pos[1] > 1.0 - REF_COORD_TOL && pos[2] < REF_COORD_TOL)
                    return true;
                break;
            case 3:
                if (pos[0] < REF_COORD_TOL && pos[1] > 1.0 - REF_COORD_TOL && pos[2] < REF_COORD_TOL)
                    return true;
                break;
            case 4:
                if (pos[0] < REF_COORD_TOL && pos[1] < REF_COORD_TOL && pos[2] > 1.0 - REF_COORD_TOL)
                    return true;
                break;
            case 5:
                if (pos[0] > 1.0 - REF_COORD_TOL && pos[1] < REF_COORD_TOL && pos[2] > 1.0 - REF_COORD_TOL)
                    return true;
                break;
            case 6:
                if (pos[0] > 1.0 - REF_COORD_TOL && pos[1] > 1.0 - REF_COORD_TOL && pos[2] > 1.0 - REF_COORD_TOL)
                    return true;
                break;
            case 7:
                if (pos[0] < REF_COORD_TOL && pos[1] > 1.0 - REF_COORD_TOL && pos[2] > 1.0 - REF_COORD_TOL)
                    return true;
                break;
            }
            return false;
        }

        static inline bool isOnFace(const RealVector &pos, Index numInterface)
        {
            assert(numInterface < 6);
            switch (numInterface)
            {
            case 0:
                if (pos[2] < REF_COORD_TOL)
                    return true;
                break;
            case 1:
                if (pos[1] < REF_COORD_TOL)
                    return true;
                break;
            case 2:
                if (pos[0] < REF_COORD_TOL)
                    return true;
                break;
            case 3:
                if (pos[0] > 1.0 - REF_COORD_TOL)
                    return true;
                break;
            case 4:
                if (pos[1] > 1.0 - REF_COORD_TOL)
                    return true;
                break;
            case 5:
                if (pos[2] > 1.0 - REF_COORD_TOL)
                    return true;
                break;
            }
            return false;
        }

        static inline Index edges(Index edge, Index vert)
        {
            assert(edge < 12);
            assert((vert == 0) || (vert == 1));

            static const Index e[12][2] = {
                {0, 1},
                {1, 2},
                {3, 2},
                {0, 3},
                {0, 4},
                {1, 5},
                {2, 6},
                {3, 7},
                {4, 5},
                {5, 6},
                {7, 6},
                {4, 7}};

            return e[edge][vert];
        }

        static inline Index faces(Index face, Index vert)
        {
            assert((face < 6));
            assert((vert < 4));

            static const Index f[6][4] = {
                {0, 1, 2, 3},
                {0, 1, 5, 4},
                {0, 3, 7, 4},
                {1, 2, 6, 5},
                {3, 2, 6, 7},
                {4, 5, 6, 7}};

            return f[face][vert];
        }

        static inline void posEdgeTransformed(const RealVector &posRes, RealVector &pos,
                                              Index numEdge, Index numOrientation)
        {
            assert(isOnEdge(posRes, numEdge));
            assert(numOrientation < 2);

            pos = posRes;

            if (numOrientation == 1)
                switch (numEdge)
                {
                case 0:
                    pos[0] = 1 - pos[0];
                    break;
                case 1:
                    pos[1] = 1 - pos[1];
                    break;
                case 2:
                    pos[0] = 1 - pos[0];
                    break;
                case 3:
                    pos[1] = 1 - pos[1];
                    break;
                case 4:
                    pos[2] = 1 - pos[2];
                    break;
                case 5:
                    pos[2] = 1 - pos[2];
                    break;
                case 6:
                    pos[2] = 1 - pos[2];
                    break;
                case 7:
                    pos[2] = 1 - pos[2];
                    break;
                case 8:
                    pos[0] = 1 - pos[0];
                    break;
                case 9:
                    pos[1] = 1 - pos[1];
                    break;
                case 10:
                    pos[0] = 1 - pos[0];
                    break;
                case 11:
                    pos[1] = 1 - pos[1];
                    break;
                }
        }

        static inline void posFaceTransformed(const RealVector &posRef, RealVector &pos,
                                              Index numFace, Index numOrientation)
        {
            assert(isOnFace(posRef, numFace));
            assert(numOrientation < 8);

            RealVector posTmp1;
            RealVector posTmp2;

            _pos2PosRefFace(posRef, posTmp1, numFace);
            RQuadrangle::posFaceTransformed(posTmp1, posTmp2, 0, numOrientation);
            _posRefFace2Pos(posTmp2, pos, numFace);
        }
    };

} // OndoMathX
