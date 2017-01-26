/* -*- c-basic-offset:2; tab-width:2; indent-tabs-mode:nil -*-
 *
 *
 */

#ifndef _BLOCK_STENCIL_H
#define _BLOCK_STENCIL_H

#include "stenseal/stencil.h"

namespace stenseal
{

  /**
   * Class representing a block of `Stencil`s with `height` rows, each `width`
   * wide.
   */
  template <int width, int height>
  class StencilArray
  {
  private:
    const std::array<Stencil<width>,height> stencils;

  public:

    /**
     * Constructor. Called with a matching number of `height` `Stencil`s as
     * input.
     */
    template <typename... Ss>
    constexpr  StencilArray(const Ss... s);

    /**
     * Access a specific row of the block.
     */
    inline constexpr const Stencil<width>& operator[](int row) const;
  };


  //---------------------------------------------------------------------------
  // Implementations
  //---------------------------------------------------------------------------

  template <int width, int height>
  template <typename... Ss>
  constexpr StencilArray<width,height>::StencilArray(const Ss... s)
    : stencils{s...}
  {
    static_assert(sizeof...(Ss) == height, "wrong number of arguments");
  }

  template <int width, int height>
  inline constexpr const Stencil<width>& StencilArray<width,height>::operator[](int row) const {
    return stencils[row];
  }

}

#endif /* _BLOCK_STENCIL_H */
