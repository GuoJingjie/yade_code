/*************************************************************************
*  2022 Tóth János                                                       *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include "Parser.hpp"
#include "VerifyMacros.hpp"

#include <boost/filesystem.hpp>

#include <stdarg.h>
#include <stdio.h>

namespace yade {
namespace ymport {
	namespace foamfile {

		Parser::Parser(const boost::filesystem::path& path, const std::string& klass, const std::string& object)
		        : m_klass(klass)
		        , m_object(object)
		{
			m_lexer = std::make_unique<Lexer>(path);

			parseHeader();
		}

		void Parser::parseHeader()
		{
			expect({ "FoamFile" });
			expect({ '{' });

			enum MustCheck { VERSION = 0, FORMAT, OBJECT, CLASS };

			bool isOk[4] = { false };

			// Token::END will eventually appear.
			while (true) {
				Token token = m_lexer->getNextToken();

				if (token.isEnd()) { error("Unexpected END."); }

				if (token.asWordEquals("version")) {
					// NOTE: Maybe compare as a word?
					expect({ 2.0 });
					expect({ ';' });
					isOk[VERSION] = true;
				} else if (token.asWordEquals("format")) {
					expect({ "ascii" });
					expect({ ';' });
					isOk[FORMAT] = true;
				} else if (token.asWordEquals("class")) {
					expect({ m_klass });
					expect({ ';' });
					isOk[CLASS] = true;
				} else if (token.asWordEquals("object")) {
					expect({ m_object });
					expect({ ';' });
					isOk[OBJECT] = true;
				} else if (token.asWordEquals("location")) {
					skipUntil({ ';' });
				} else if (token.asCharEquals('}')) {
					break;
				} else {
					error("Unexpected token: %s", token.debugString().c_str());
				}
			}

			if (!isOk[VERSION]) { error("No version was provided in the header."); }

			if (!isOk[FORMAT]) { error("No format was provided in the header."); }

			if (!isOk[CLASS]) { error("No class was provided in the header."); }

			if (!isOk[OBJECT]) { error("No object was provided in the header."); }
		}

		void Parser::expect(Token token)
		{
			Token current = m_lexer->getNextToken();
			if (!current.isEqual(token)) { error("Expected: %s, got: %s", token.debugString().c_str(), current.debugString().c_str()); }
		}

		void Parser::skip(int n)
		{
			for (int i = 0; i < n; ++i) {
				m_lexer->getNextToken();
			}
		}

		void Parser::skipUntil(Token token)
		{
			Token current;
			while (!current.isEqual(token)) {
				current = m_lexer->getNextToken();

				if (token.isEnd()) { error("Unexpected END."); }
			}
		}

		std::string Parser::getWord()
		{
			Token token = m_lexer->getNextToken();

			if (!token.isWord()) { error("Expected 'WODD'', got: %s", token.debugString().c_str()); }

			return token.asWord();
		}

		int Parser::getInt()
		{
			Token token = m_lexer->getNextToken();

			if (!token.isInt()) { error("Expected 'INT'', got: %s", token.debugString().c_str()); }

			return token.asInt();
		}

		double Parser::getDouble()
		{
			Token token = m_lexer->getNextToken();

			if (!token.isDouble()) { error("Expected 'DOUBLE'', got: %s", token.debugString().c_str()); }

			return token.asDouble();
		}

		double Parser::getNumber()
		{
			Token token = m_lexer->getNextToken();

			if (!token.isNumber()) { error("Expected 'NUMBER'', got: %s", token.debugString().c_str()); }

			return token.asNumber();
		}

		void Parser::error(const char* fmt, ...)
		{
			va_list args;
			va_start(args, fmt);
			size_t size = vsnprintf(NULL, 0, fmt, args);
			va_end(args);

			std::vector<char> data(size + 1);

			va_start(args, fmt);
			vsnprintf(data.data(), data.size(), fmt, args);
			va_end(args);

			std::string message = m_lexer->path().filename().string() + ":" + std::to_string(m_lexer->lineNumber()) + ": ";
			message.append(std::string(data.begin(), data.end()));

			throw std::runtime_error(message);
		}

	} // namespace foamfile
} // namespace ymport
} // namespace yade
