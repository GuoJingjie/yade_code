/*************************************************************************
*  2022 Tóth János                                                       *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

// The lexer converts the character stream into a token stream.
// The output of the lexer goes to the input of the parser.
//
// The "getNextToken()" will return "END" after the last non-end token and
// it can be called more than one time after "END", it will return with "END".

#pragma once

#include "Token.hpp"
#include "VerifyMacros.hpp"

#include <boost/filesystem.hpp>

namespace yade {
namespace ymport {
	namespace foamfile {

		class Lexer {
		public:
			Lexer(const boost::filesystem::path& path);
			~Lexer();

			Token getNextToken();

			const boost::filesystem::path& path() const { return m_path; }
			std::uint32_t                  lineNumber() const { return m_lineNumber; }

		private:
			void readFullFile();

			void removeSingleLineComments();
			void removeMultiLineComments();

			void consumeWhiteSpaces();

			bool isWhiteSpace(char c);

			void error(const char* fmt, ...);

			std::uint32_t           m_currentCharIndex { 0 };
			std::uint32_t           m_lineNumber { 1 };
			boost::filesystem::path m_path;
			std::string             m_data;
			const char*             m_oldLocale;
		};

	} // namespace foamfile
} // namespace ymport
} // namespace yade
