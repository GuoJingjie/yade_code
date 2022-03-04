/*************************************************************************
*  2022 Tóth János                                                       *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include "Lexer.hpp"
#include "VerifyMacros.hpp"

#include <boost/regex.hpp>

#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include <locale.h>
#include <stdarg.h>
#include <stdio.h>

namespace yade {
namespace ymport {
	namespace foamfile {

		Lexer::Lexer(const boost::filesystem::path& path)
		        : m_path(path)
		{
			m_oldLocale = setlocale(LC_NUMERIC, nullptr);

			// Force '.' as the radix point.
			// This is neccessary for the string->double conversion,
			// e.g.: my system uses ',' as a decimal separator.
			setlocale(LC_NUMERIC, "C");

			// NOTE: Read in chunks? If the file is too large loading it in the memory can be a problem.
			readFullFile();

			removeSingleLineComments();
			removeMultiLineComments();

			m_currentCharIndex = 0;
			m_lineNumber       = 1;
		}

		Lexer::~Lexer() { setlocale(LC_NUMERIC, m_oldLocale); }

		// Do it in a POSIX way, this must work regardelss of distro or version.
		void Lexer::readFullFile()
		{
			int fd = open(m_path.string().c_str(), O_RDONLY);
			if (fd == -1) { error("Unable to open file."); }

			struct stat st;
			fstat(fd, &st);

			VERIFY(st.st_size > 0);

			try {
				m_data.resize(st.st_size);
			} catch (...) {
				close(fd);
				error("Unable to allocate buffer.");
			}

			ssize_t n = read(fd, &m_data[0], m_data.size());

			if (n == -1 || (size_t)n != m_data.size()) {
				close(fd);
				error("Unable to read from file.");
			}

			close(fd);
		}

		bool Lexer::isWhiteSpace(char c) { return c == ' ' || c == '\t' || c == '\n'; }

		void Lexer::consumeWhiteSpaces()
		{
			for (; m_currentCharIndex < m_data.length(); m_currentCharIndex++) {
				char curr = m_data[m_currentCharIndex];
				if (curr == '\n') { m_lineNumber++; }

				if (!isWhiteSpace(curr)) { break; }
			}
		}

		Token Lexer::getNextToken()
		{
			consumeWhiteSpaces();
			if (m_currentCharIndex >= m_data.length() - 1) { return { Token::Type::END }; }

			// Blob consists one or more characters, they will become 'word's, 'int's or 'double's.
			auto isBlobChar = [](char c) {
				return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || (c >= '0' && c <= '9') || c == '_' || c == '-' || c == '+'
				        || c == '.';
			};

			// Char is exactly one character.
			auto isCharChar = [](char c) { return c == ';' || c == '{' || c == '}' || c == '(' || c == ')' || c == '<' || c == '>'; };

			enum State { EMPTY = 0, BLOB, STRING };

			State       state = EMPTY;
			std::string currentTokenData;
			for (; m_currentCharIndex < m_data.length(); m_currentCharIndex++) {
				char curr = m_data[m_currentCharIndex];

				if (isWhiteSpace(curr)) { break; }

				if (state == EMPTY) {
					if (isBlobChar(curr)) {
						state = BLOB;
						currentTokenData += curr;
					} else if (isCharChar(curr)) {
						m_currentCharIndex++;
						return { curr };
					} else if (curr == '"') {
						currentTokenData += curr;
						state = STRING;
					} else {
						error("Unknown character: '%c'", curr);
					}
				} else if (state == BLOB) {
					if (isBlobChar(curr)) {
						currentTokenData += curr;
					} else {
						break;
					}
					// string will "eat" everything
				} else if (state == STRING) {
					if (curr == '"') {
						currentTokenData += curr;
						m_currentCharIndex++;
						break;
					} else {
						currentTokenData += curr;
					}
				} else {
					UNREACHABLE();
				}
			}

			VERIFY(currentTokenData.length() > 0);

			boost::regex re_int("[-+]?[0-9]+");
			boost::regex re_dbl("[-+]?[0-9]?[.][0-9]+");

			if (boost::regex_match(currentTokenData, re_int)) {
				int i = 0;
				try {
					i = std::stoi(currentTokenData);
				} catch (std::exception& ex) {
					error("%s", ex.what());
				}

				return { i }; // int
			} else if (boost::regex_match(currentTokenData, re_dbl)) {
				double d = 0;
				try {
					d = std::stod(currentTokenData);
				} catch (std::exception& ex) {
					error("%s", ex.what());
				}

				return { d }; // double
			} else {
				return { currentTokenData }; // word
			}
		}

		void Lexer::removeSingleLineComments()
		{
			size_t len       = 0;
			bool   inComment = false;
			size_t start     = 0;
			m_lineNumber     = 1;
			for (size_t i = 0; i < m_data.length() - 1; ++i) {
				char c1 = m_data[i];
				char c2 = m_data[i + 1];

				if (c1 == '\n') { m_lineNumber++; }

				if (c1 == '/' && c2 == '/') {
					inComment = true;
					start     = i;
					len       = 0;
				}

				if (inComment) {
					len++;
					if (c1 == '\n') {
						m_data.replace(start, len, "\n");
						inComment = false;
						i         = 0;
					}
				}
			}

			if (inComment && m_data[m_data.length() - 1] == '\n') {
				len++;
				m_data.replace(start, len, "\n");
				inComment = false;
			}

			if (inComment) { error("Non-terminated single line comment."); }
		}

		void Lexer::removeMultiLineComments()
		{
			size_t len              = 0;
			bool   inComment        = false;
			size_t numberOfNewLines = 0;
			size_t start            = 0;
			m_lineNumber            = 1;
			for (size_t i = 0; i < m_data.length() - 1; ++i) {
				char c1 = m_data[i];
				char c2 = m_data[i + 1];

				if (c1 == '\n') { m_lineNumber++; }

				if (c1 == '/' && c2 == '*') {
					numberOfNewLines = 0;
					inComment        = true;
					start            = i;
					len              = 0;
				}

				if (inComment) {
					len++;

					if (c1 == '\n') {
						numberOfNewLines++;
					} else if (c1 == '*' && c2 == '/') {
						len++;
						std::string newLines = "";

						for (size_t j = 0; j < numberOfNewLines; ++j) {
							newLines += '\n';
						}

						m_data.replace(start, len, newLines);
						inComment = false;
						i         = 0;
					}
				}
			}

			if (inComment) { error("Non-terminated multi line comment."); }
		}

		void Lexer::error(const char* fmt, ...)
		{
			va_list args;

			va_start(args, fmt);
			size_t size = vsnprintf(NULL, 0, fmt, args);
			va_end(args);

			std::vector<char> data(size + 1);

			va_start(args, fmt);
			vsnprintf(data.data(), data.size(), fmt, args);
			va_end(args);

			std::string message = m_path.filename().string() + ":" + std::to_string(m_lineNumber) + ": ";
			message.append(std::string(data.begin(), data.end()));

			throw std::runtime_error(message);
		}

	} // namespace foamfile
} // namespace ymport
} // namespace yade
