/*************************************************************************
*  2022 Tóth János                                                       *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

// Token is a lexical unit that the lexer outputs for the parser from the input characters.
// E.g.: "format      ascii;" -> WORD WORD ';'
//            characters           tokens
//
// Valid tokens are:
//   + EMPTY (default)
//   + WORD (one or more characters)
//   + CHAR (one character)
//   + INT (integer number)
//   + DOUBLE (double precision floating point number)
//   + END (last token of the file)

#pragma once

#include "VerifyMacros.hpp"

#include <string>

namespace yade {
namespace ymport {
	namespace foamfile {

		class Token {
		public:
			enum Type { EMPTY = 0, WORD, CHAR, INT, DOUBLE, END };

			// NOTE: It wolud be nice to include the "string (word)" value of the token as well
			// in form of a "const char * s", create it with strdup() and release it with free().
			// I tried to do that, but not sure when or where to use "free()".
			union Value {
				double d;
				int    i;
				char   c;
			};

			Token(Type type = Type::EMPTY)
			{
				if (type == Type::END) {
					m_type = Type::END;
				} else {
					m_type = Type::EMPTY;
				}
			}

			Token(const std::string& s)
			        : m_type(Type::WORD)
			        , m_str_value(s)
			{
			}

			Token(char c)
			        : m_type(Type::CHAR)
			{
				m_value.c = c;
			}

			Token(int i)
			        : m_type(Type::INT)
			{
				m_value.i = i;
			}

			Token(double d)
			        : m_type(Type::DOUBLE)
			{
				m_value.d = d;
			}

			Type type() const { return m_type; }

			bool isEmpty() const { return m_type == Type::EMPTY; }

			bool isEnd() const { return m_type == Type::END; }

			bool isWord() const { return m_type == Type::WORD; }

			bool isChar() const { return m_type == Type::CHAR; }

			bool isInt() const { return m_type == Type::INT; }

			bool isDouble() const { return m_type == Type::DOUBLE; }

			bool isNumber() const { return isInt() || isDouble(); }

			const std::string& asWord() const
			{
				VERIFY(isWord());
				return m_str_value;
			};

			char asChar() const
			{
				VERIFY(isChar());
				return m_value.c;
			};

			int asInt() const
			{
				VERIFY(isInt());
				return m_value.i;
			};

			double asDouble() const
			{
				VERIFY(isDouble());
				return m_value.d;
			};

			double asNumber() const
			{
				VERIFY(isNumber());
				if (isInt()) {
					return (double)asInt();
				} else {
					return asDouble();
				}
			};

			bool isEqual(const Token& t)
			{
				if (t.isWord()) {
					return asWordEquals(t.asWord());
				} else if (t.isInt()) {
					return asIntEquals(t.asInt());
				} else if (t.isChar()) {
					return asCharEquals(t.asChar());
				} else if (t.isDouble()) {
					return asDoubleEquals(t.asDouble());
				} else if (t.type() == Type::EMPTY || t.type() == Type::END) {
					return type() == t.type();
				} else {
					UNREACHABLE();
					return false;
				}
			}

			bool asWordEquals(const std::string& s) { return isWord() && asWord() == s; }

			bool asCharEquals(char c) { return isChar() && asChar() == c; }

			bool asIntEquals(int i) { return isInt() && asInt() == i; }

			bool asDoubleEquals(double d) { return isDouble() && asDouble() == d; }

			bool asNumberEquals(double d) { return isNumber() && asNumber() == d; }

			std::string debugString() const
			{
				switch (m_type) {
					case Type::EMPTY: return "EMPTY";

					case Type::WORD: return "WORD(" + asWord() + ")";

					case Type::CHAR: return "CHAR(" + std::string(1, asChar()) + ")";

					case Type::INT: return "INT(" + std::to_string(asInt()) + ")";

					case Type::DOUBLE: return "DOUBLE(" + std::to_string(asDouble()) + ")";

					case Type::END: return "END";

					default: UNREACHABLE(); return "";
				}
			}

		private:
			Type        m_type { EMPTY };
			Value       m_value {};
			std::string m_str_value;
		};

	} // namespace foamfile
} // namespace ymport
} // namespace yade
