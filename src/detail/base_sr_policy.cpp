/* Copyright 2017-2018 PaGMO development team

This file is part of the PaGMO library.

The PaGMO library is free software; you can redistribute it and/or modify
it under the terms of either:

  * the GNU Lesser General Public License as published by the Free
    Software Foundation; either version 3 of the License, or (at your
    option) any later version.

or

  * the GNU General Public License as published by the Free Software
    Foundation; either version 3 of the License, or (at your option) any
    later version.

or both in parallel, as here.

The PaGMO library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received copies of the GNU General Public License and the
GNU Lesser General Public License along with the PaGMO library.  If not,
see https://www.gnu.org/licenses/. */

#include <cassert>
#include <cmath>
#include <stdexcept>
#include <string>

#include <boost/serialization/variant.hpp>
#include <boost/variant/get.hpp>

#include <pagmo/detail/base_sr_policy.hpp>
#include <pagmo/exceptions.hpp>
#include <pagmo/s11n.hpp>

namespace pagmo
{

namespace detail
{

// Helper to verify the ctor from a fractional rate.
void base_sr_policy::verify_fp_ctor() const
{
    assert(m_migr_rate.which() == 1);

    const auto rate = boost::get<double>(m_migr_rate);

    if (!std::isfinite(rate) || rate < 0. || rate > 1.) {
        pagmo_throw(std::invalid_argument,
                    "Invalid fractional migration rate specified in the constructor of a replacement/selection "
                    "policy: the rate must be in the [0., 1.] range, but it is "
                        + std::to_string(rate) + " instead");
    }
}

template <typename Archive>
void base_sr_policy::serialize(Archive &ar, unsigned)
{
    detail::archive(ar, m_migr_rate);
}

// NOTE: explicitly instantiate the specialisation of the serialize() member
// function that we will be using. Note that if we add further archive types
// in s11n.hpp, we'll have to add the specialisations here as well.
template void base_sr_policy::serialize(boost::archive::binary_iarchive &, unsigned);
template void base_sr_policy::serialize(boost::archive::binary_oarchive &, unsigned);
template void base_sr_policy::serialize(boost::archive::text_iarchive &, unsigned);
template void base_sr_policy::serialize(boost::archive::text_oarchive &, unsigned);

} // namespace detail

} // namespace pagmo
