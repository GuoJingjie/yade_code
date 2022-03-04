/*************************************************************************
*  2022 Tóth János                                                       *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

//	The parser interprets the syntax of the file and parses the values of the tokens.
// 	This is the base class and implements all the common functions.

#pragma once

#include "Lexer.hpp"
#include "Token.hpp"

#include <memory>
#include <string>

namespace yade {
namespace ymport {
	namespace foamfile {

		class Parser {
		public:
			Parser(const boost::filesystem::path& path, const std::string& klass, const std::string& object);
			virtual ~Parser() {};

			virtual void parse() = 0;

		protected:
			void expect(Token token);

			void skip(int n);
			void skipUntil(Token token);

			std::string getWord();
			int         getInt();
			double      getDouble();
			double      getNumber();

			void parseHeader();

			void error(const char* fmt, ...);

			std::unique_ptr<Lexer> m_lexer { nullptr };
			std::string            m_klass;
			std::string            m_object;
		};

	} // namespace foamfile
} // namespace ymport
} // namespace yade
