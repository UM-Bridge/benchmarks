#pragma once

// -----------------------------------------------------------------------------------------//
#include <memory>
#include <vector>
#include "LinearAlgebraLibrary.h"
// -----------------------------------------------------------------------------------------//

// -----------------------------------------------------------------------------------------//
// Class definition
namespace OndoMathX
{

	namespace LAL
	{

		class VectorSequence
		{

		public:
			// ---------------------------------------------------------------------------------//
			/*! \brief Default constructor.
			 */
			VectorSequence(Index NIterates)
			{
				// Resizing.
				_pData.resize(NIterates);
				for (Index i = 0; i < _pData.size(); i++)
					_pData[i] = nullptr;
			}

			/*! \brief Copy constructor is deleted
			 */
			VectorSequence(const VectorSequence &aVectorSequence) = delete;
			// ---------------------------------------------------------------------------------//

			// ---------------------------------------------------------------------------------//
			/*! \brief Accessing the number of time step.
			 */
			Index getNumIterates() const
			{
				return _pData.size();
			}

			/*! \brief Accessing the size of an unknown.
			 \param iUnknown is an index of an unknown.
			 */
			Index getDimension() const
			{
				if (_pData[0] == nullptr)
					return 0;
				else
					return LAL::getDimension(*_pData[0]);
			}

			/*! \brief Allocating an unknown.
				\param Size is the size of the unknown after allocation.
				\param iUnknown is an index of an unknown.
				\param allocateBuffer indicates if a buffer of the same size has to be allocated for intermediate computation step.
			*/
			void Allocate(Index Dimension)
			{
				for (Index iStep = 0; iStep < _pData.size(); iStep++)
				{
					_pData[iStep] = std::make_unique<Vector>();

					LAL::Allocate(*_pData[iStep], Dimension);
				}
			}
			// ---------------------------------------------------------------------------------//

			// ---------------------------------------------------------------------------------//

			/*! \brief Extracting the raw pointer to a vector of DoF of the solution.
				\param Step is a specific time step.
				\param Unknown is the corresponding unknown index.
				\return A raw pointer to a vector of DoF of the solution.
			*/
			LAL::Vector &getVector(Index i = 0)
			{
				return *_pData[i].get();
			}

			// ---------------------------------------------------------------------------------//

			// ---------------------------------------------------------------------------------//
			/*! \brief Switches between the solutions between time steps
			 */
			void Swap()
			{
				// Swaping solutions.
				for (Index iStep = _pData.size() - 1; iStep > 0; iStep--)
					std::swap(_pData[iStep], _pData[iStep - 1]);
			}

			/*! \brief Switches between the solutions between time steps
			 */
			void Reverse_Swap()
			{
				// Swaping solutions.
				for (Index iStep = 1; iStep < _pData.size(); iStep++)
					std::swap(_pData[iStep], _pData[iStep - 1]);
			}
			// ---------------------------------------------------------------------------------//

			// ---------------------------------------------------------------------------------//
			/*! \brief Reversing the last operations performed on time steps.
			 */
			void Reverse()
			{
				if (_pData.size() >= 2)
					std::swap(_pData[_pData.size() - 1], _pData[_pData.size() - 2]);
			}
			// ---------------------------------------------------------------------------------//

		private:
			/*! \brief A solution is defined as a vector of vector of pointers towards raw data containers.
			The first vector index the time steps and the second the unknown of each time steps.
			*/
			std::vector<std::unique_ptr<Vector>> _pData;
		};

	} // LAL

} // OndoMathX
