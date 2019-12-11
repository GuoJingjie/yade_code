/*************************************************************************
*  2019 Janek Kozicki                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

// Note: include this file only for ThinRealWrapper<UnderlyingReal>

#ifndef YADE_REAL_NUMERIC_LIMITS_HPP
#define YADE_REAL_NUMERIC_LIMITS_HPP

namespace std {
template <> struct numeric_limits<ThinRealWrapper<UnderlyingReal>> {
	constexpr static const auto& is_specialized    = std::numeric_limits<UnderlyingReal>::is_specialized;
	constexpr static const auto& is_signed         = std::numeric_limits<UnderlyingReal>::is_signed;
	constexpr static const auto& is_integer        = std::numeric_limits<UnderlyingReal>::is_integer;
	constexpr static const auto& is_exact          = std::numeric_limits<UnderlyingReal>::is_exact;
	constexpr static const auto& has_infinity      = std::numeric_limits<UnderlyingReal>::has_infinity;
	constexpr static const auto& has_quiet_NaN     = std::numeric_limits<UnderlyingReal>::has_quiet_NaN;
	constexpr static const auto& has_signaling_NaN = std::numeric_limits<UnderlyingReal>::has_signaling_NaN;
	constexpr static const auto& has_denorm        = std::numeric_limits<UnderlyingReal>::has_denorm;
	constexpr static const auto& has_denorm_loss   = std::numeric_limits<UnderlyingReal>::has_denorm_loss;
	constexpr static const auto& round_style       = std::numeric_limits<UnderlyingReal>::round_style;
	constexpr static const auto& is_iec559         = std::numeric_limits<UnderlyingReal>::is_iec559;
	constexpr static const auto& is_bounded        = std::numeric_limits<UnderlyingReal>::is_bounded;
	constexpr static const auto& is_modulo         = std::numeric_limits<UnderlyingReal>::is_modulo;
	constexpr static const auto& digits            = std::numeric_limits<UnderlyingReal>::digits;
	constexpr static const auto& digits10          = std::numeric_limits<UnderlyingReal>::digits10;
	constexpr static const auto& max_digits10      = std::numeric_limits<UnderlyingReal>::max_digits10;
	constexpr static const auto& radix             = std::numeric_limits<UnderlyingReal>::radix;
	constexpr static const auto& min_exponent      = std::numeric_limits<UnderlyingReal>::min_exponent;
	constexpr static const auto& min_exponent10    = std::numeric_limits<UnderlyingReal>::min_exponent10;
	constexpr static const auto& max_exponent      = std::numeric_limits<UnderlyingReal>::max_exponent;
	constexpr static const auto& max_exponent10    = std::numeric_limits<UnderlyingReal>::max_exponent10;
	constexpr static const auto& traps             = std::numeric_limits<UnderlyingReal>::traps;
	constexpr static const auto& tinyness_before   = std::numeric_limits<UnderlyingReal>::tinyness_before;
	static inline auto           min() { return static_cast<ThinRealWrapper<UnderlyingReal>>(std::numeric_limits<UnderlyingReal>::min()); }
	static inline auto           lowest() { return static_cast<ThinRealWrapper<UnderlyingReal>>(std::numeric_limits<UnderlyingReal>::lowest()); }
	static inline auto           max() { return static_cast<ThinRealWrapper<UnderlyingReal>>(std::numeric_limits<UnderlyingReal>::max()); }
	static inline auto           epsilon() { return static_cast<ThinRealWrapper<UnderlyingReal>>(std::numeric_limits<UnderlyingReal>::epsilon()); }
	static inline auto           round_error() { return static_cast<ThinRealWrapper<UnderlyingReal>>(std::numeric_limits<UnderlyingReal>::round_error()); }
	static inline auto           infinity() { return static_cast<ThinRealWrapper<UnderlyingReal>>(std::numeric_limits<UnderlyingReal>::infinity()); }
	static inline auto           quiet_NaN() { return static_cast<ThinRealWrapper<UnderlyingReal>>(std::numeric_limits<UnderlyingReal>::quiet_NaN()); }
	static inline auto           signaling_NaN() { return static_cast<ThinRealWrapper<UnderlyingReal>>(std::numeric_limits<UnderlyingReal>::signaling_NaN()); }
	static inline auto           denorm_min() { return static_cast<ThinRealWrapper<UnderlyingReal>>(std::numeric_limits<UnderlyingReal>::denorm_min()); }
	// constexpr static const auto& float_round_style = std::numeric_limits<UnderlyingReal>::float_round_style ;
	// constexpr static const auto& float_denorm_style= std::numeric_limits<UnderlyingReal>::float_denorm_style;
};
}

#endif
