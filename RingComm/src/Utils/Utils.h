/** @addtogroup RingComm
 * @{
 */
/**
 * \namespace UtilsRC
 * \ingroup RingComm
 * \brief Different utility functions used in the implementation. Specific for Ring exchange approach.
 */
#ifndef __UTILS_H
#define __UTILS_H

#include <boost/dynamic_bitset.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <cmath>

#include "../Environment/Environment.h"
#include "../Basis/Basis.h"

namespace UtilsRC
{
  /** \brief A different definition for modulus %, accounts for negative integers.
    * \param a An integer.
    * \param b An integer.
    * \return The modulus of the two integers.
    */
  LLInt mod(LLInt a, LLInt b);
  /** \brief Computes the integer representation of a bitset object.
    * \param bs The bitset object.
    * \param l The number of sites in the system.
    * \return An integer value for the binary representation.
    */
  LLInt binary_to_int(boost::dynamic_bitset<> bs,
                      unsigned int l);
  /** \brief Binary search algorithm.
    * \param array Sorted array of integers.
    * \param len Number of integers in the array.
    * \param value Integer value to locate.
    * \return The index of the found value, 0 if unfound.
    */
  LLInt binsearch(const LLInt *array, 
                  LLInt len, 
                  LLInt value);
  /** \brief Returns the position of the Neel state of the system in computational basis.
    * \param env An instance of class Environment.
    * \param bas An instance of class Basis.
    * \return The index of the Neel state.
    */
}
#endif
/** @}*/
