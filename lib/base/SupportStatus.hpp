/*************************************************************************
*  2022 Janek Kozicki                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

/*
 * Each class is getting 4 or 5 extra functions. The goal is to show in documentation automatically retreived information from code:
 *  1. the color in the inheritance graph says: green - mature and supported, yellow - in development, red - abandoned
 *  2. 
 *
 */

#ifndef SUPPORT_STATUS_HPP
#define SUPPORT_STATUS_HPP

class SupportStatus : public Factorable {
public:
	SupportStatus() { }
	virtual ~SupportStatus() { }

	// The first returned argument is getClassName() so that self inspection in yade --test can discover an error when the
	// parent virtual function returns its info i.e. the function are not overridden in some particular class.
	// Meaning that someone forgot to add the info.
	virtual std::vector<std::string> getAuthors() const { return {getClassName(), "Janek Kozicki"}; }
	virtual std::pair<std::string, double>    getCreationYearMonth() const { return {getClassName(), 2022.06}; }
	virtual std::vector<std::string> getPublications() const { return {getClassName(), "TODO: new yade-computer-physics paper will mention this thing"}; };
	virtual std::vector<std::string> getBibtex()       const { return {getClassName(), "???"}; }
	virtual std::pair<std::string,Status> getSupportStatus    const { return {getClassName(), Status::UNKNOWN }; }
};

#endif
