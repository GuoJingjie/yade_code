/*************************************************************************
*  2022 Janek Kozicki                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

/*
 * Each class is getting 4 or 5 extra functions. The goal is to show in documentation automatically retreived information from code:
 *
 *  1. getSupportStatus() - the color in the inheritance graph says: green - mature and supported, yellow - in development, red - abandoned
 *  2. we can call getAuthors() to check who wrote this class
 *  3. we can call getPublications() to see what publications mention or use this class
 *     TODO: use this function in documentation to print this info in references next to each class
 *  4. we can call getCreationYearMonth() to know when the class was written for the first time
 *  5. getBibtex() sounds useful, but we will manage to put it everywhere ? We could go through bibliography
 *     and move to here from there, then generate bibliography. Not sure.
 *
 */

#ifndef SUPPORT_HPP
#define SUPPORT_HPP

class AuthorsPublicationsStatus : public Factorable {
public:
	AuthorsPublicationsStatus() { }
	virtual ~AuthorsPublicationsStatus() { }

	enum class SupportStatus { MATURE, DEVELOPMENT, ABANDONED, UNKNOWN };

	// The first returned argument is always getClassName() so that self inspection in yade --test can discover an error when the
	// parent virtual function returns its info i.e. the function are not overridden in some particular class.
	// Meaning that someone forgot to add the info.
	virtual std::pair<std::string,SupportStatus> getSupportStatus    const { return {getClassName(), SupportStatus::DEVELOPMENT }; }
	virtual std::vector<std::string> getAuthors() const { return {getClassName(), "Janek Kozicki"}; }
	virtual std::pair<std::string, double>    getCreationYearMonth() const { return {getClassName(), 2022.06}; }
	virtual std::vector<std::string> getPublications() const { return {getClassName(), "TODO: new yade-computer-physics paper will mention this thing"}; };
	virtual std::vector<std::string> getBibtex()       const { return {getClassName(), "???"}; }
};

#endif
