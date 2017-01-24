/* -*- c-basic-offset:2; tab-width:2; indent-tabs-mode:nil -*-
 *
 *
 */

#ifndef _BLOCK_STENCIL_H
#define _BLOCK_STENCIL_H

#include "stenseal/stencil.h"

namespace stenseal
{

  template <int width, int height>
  struct BlockStencil
  {
    const std::array<Stencil<width>,height> stencils;

    template <typename... Ss>
    constexpr  BlockStencil(const Ss... s)
      : stencils{s...}
    {
      static_assert(sizeof...(Ss) == height, "wrong number of arguments");
    }

    constexpr const Stencil<width>& operator[](int row) const {
      return stencils[row];
    }
  };
}

#endif /* _BLOCK_STENCIL_H */
