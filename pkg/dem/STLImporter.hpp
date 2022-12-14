/*************************************************************************
*  Copyright (C) 2008 by Sergei Dorofeenko				 *
*  sega@users.berlios.de                                                 *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/
#pragma once
#include <core/Body.hpp>
#include <core/BodyContainer.hpp>
#include <pkg/common/Facet.hpp>

namespace yade { // Cannot have #include directive inside.

class STLImporter {
public:
	vector<shared_ptr<Body>> import(const char*);
	DECLARE_LOGGER;
};

template <class T> std::pair<T, T> minmax(const T& a, const T& b) { return (a < b) ? std::pair<T, T>(a, b) : std::pair<T, T>(b, a); }


class STLReader {
public:
	// if it is binary there are 80 char of comment, the number fn of faces and then exactly fn*4*3 bytes.
	enum { STL_LABEL_SIZE = 80 };

	template <class OutV, class OutE, class OutF, class OutN> bool open(const char* filename, OutV vertices, OutE edges, OutF facets, OutN normals);

	float tolerance;

protected:
	template <class OutV, class OutE, class OutF, class OutN> bool open_binary(const char* filename, OutV vertices, OutE edges, OutF facets, OutN normals);

	template <class OutV, class OutE, class OutF, class OutN> bool open_ascii(const char* filename, OutV vertices, OutE edges, OutF facets, OutN normals);

	struct Vrtx {
		float pos[3];
		Vrtx() { }
		Vrtx(float x, float y, float z)
		{
			pos[0] = x;
			pos[1] = y;
			pos[2] = z;
		}
		bool operator<(const Vrtx& v) const
		{
#if (YADE_REAL_BIT <= 64)
			return memcmp(pos, v.pos, 3 * sizeof(float)) < 0;
#else
			return std::tie(pos[0], pos[1], pos[2]) < std::tie(v.pos[0], v.pos[1], v.pos[2]);
#endif
		}
		float  operator[](int id) const { return pos[id]; }
		float& operator[](int id) { return pos[id]; }
	};

	class IsDifferent {
		float tolerance;

	public:
		IsDifferent(float _tolerance)
		        : tolerance(_tolerance)
		{
		}
		bool operator()(const Vrtx& v1, const Vrtx& v2)
		{
			if (std::abs(v1[0] - v2[0]) < tolerance && std::abs(v1[1] - v2[1]) < tolerance && std::abs(v1[2] - v2[2]) < tolerance) return false;
			return true;
		}
	};
};

template <class OutV, class OutE, class OutF, class OutN> bool STLReader::open(const char* filename, OutV vertices, OutE edges, OutF facets, OutN normals)
{
	FILE* fp;
	bool  binary = false;
	fp           = fopen(filename, "r");
	if (fp == NULL) return false;

	/* Find size of file */
	fseek(fp, 0, SEEK_END);
	int file_size = ftell(fp);
	int facenum;
	/* Check for binary or ASCII file */
	fseek(fp, STL_LABEL_SIZE, SEEK_SET);
	size_t res                = fread(&facenum, sizeof(int), 1, fp);
	int    expected_file_size = STL_LABEL_SIZE + 4 + (sizeof(short) + 4 * sizeof(float)) * facenum;
	if (file_size == expected_file_size) binary = true;
	unsigned char tmpbuf[128];
	res = fread(tmpbuf, sizeof(tmpbuf), 1, fp);
	if (res) {
		for (size_t i = 0; i < sizeof(tmpbuf); i++) {
			if (tmpbuf[i] > 127) {
				binary = true;
				break;
			}
		}
	}
	// Now we know if the stl file is ascii or binary.
	fclose(fp);
	if (binary) return open_binary(filename, vertices, edges, facets, normals);
	else
		return open_ascii(filename, vertices, edges, facets, normals);
}

template <class OutV, class OutE, class OutF, class OutN> bool STLReader::open_ascii(const char* filename, OutV vertices, OutE edges, OutF facets, OutN normals)
{
	FILE* fp;
	fp = fopen(filename, "r");
	if (fp == NULL) return false;

	/* Skip the first line of the file */
	while (getc(fp) != '\n')
		;

	vector<Vrtx>             vcs;
	std::set<pair<int, int>> egs;

	/* Read a single facet from an ASCII .STL file 
   * ASCII file looks like:
   *   solid name
   *    facet normal n1 n2 n3
   *     outer loop
   *     vertex p1x p1y p1z
   *     vertex p2x p2y p2z
   *     vertex p3x p3y p3z
   *     endloop
   *    endfacet
   *   endsolid name
   */
	int r = 0; //just to escape "warning: ignoring return value of 'int fscanf(FILE*, const char*, ...)',"
	while (!feof(fp)) {
		float n[3];
		Vrtx  v[3];
		r = fscanf(fp, "%*s %*s %f %f %f\n", &n[0], &n[1], &n[2]);
		r += fscanf(fp, "%*s %*s");
		r += fscanf(fp, "%*s %f %f %f\n", &v[0][0], &v[0][1], &v[0][2]);
		r += fscanf(fp, "%*s %f %f %f\n", &v[1][0], &v[1][1], &v[1][2]);
		r += fscanf(fp, "%*s %f %f %f\n", &v[2][0], &v[2][1], &v[2][2]);
		r += fscanf(fp, "%*s"); // end loop
		r += fscanf(fp, "%*s"); // end facet
		if (feof(fp)) {
			r = 1; // last time reading procedures from above are called, the end-of-file was reached so r got negative; assigning some value > 0 let's succed this function
			break;
		}

		int vid[3];
		for (int i = 0; i < 3; ++i) {
			(normals++)              = n[i];
			bool        is_different = true;
			IsDifferent isd(tolerance);
			int         j = 0;
			for (int ej = vcs.size(); j < ej; ++j)
				if (!(is_different = isd(v[i], vcs[j]))) break;
			if (is_different) {
				vid[i] = vcs.size();
				vcs.push_back(v[i]);
			} else
				vid[i] = j;
			(facets++) = vid[i];
		}
		egs.insert(minmax(vid[0], vid[1]));
		egs.insert(minmax(vid[1], vid[2]));
		egs.insert(minmax(vid[2], vid[0]));
	}
	fclose(fp);

	for (vector<Vrtx>::iterator it = vcs.begin(), end = vcs.end(); it != end; ++it) {
		(vertices++) = (*it)[0];
		(vertices++) = (*it)[1];
		(vertices++) = (*it)[2];
	}
	for (std::set<pair<int, int>>::iterator it = egs.begin(), end = egs.end(); it != end; ++it) {
		(edges++) = it->first;
		(edges++) = it->second;
	}
	return (r > 0);
}

template <class OutV, class OutE, class OutF, class OutN>
bool STLReader::open_binary(const char* filename, OutV vertices, OutE edges, OutF facets, OutN normals)
{
	FILE* fp;
	fp = fopen(filename, "rb");
	if (fp == NULL) { return false; }

	int facenum;
	fseek(fp, STL_LABEL_SIZE, SEEK_SET);
	int res = fread(&facenum, sizeof(int), 1, fp);

	vector<Vrtx>             vcs;
	std::set<pair<int, int>> egs;

	if (res) {
		// For each triangle read the normal, the three coords and a short set to zero
		for (int i = 0; i < facenum; ++i) {
			short attr;
			float n[3];
			Vrtx  v[3];
			res = fread(&n, 3 * sizeof(float), 1, fp);
			res = fread(&v, sizeof(Vrtx), 3, fp);
			res = fread(&attr, sizeof(short), 1, fp);

			//FIXME: ???????????? ???????????????????????? ???????? ?? open_ascii
			int vid[3];
			for (int l = 0; l < 3; ++l) {
				(normals++)              = n[l];
				bool        is_different = true;
				IsDifferent isd(tolerance);
				int         j = 0;
				for (int ej = vcs.size(); j < ej; ++j)
					if (!(is_different = isd(v[l], vcs[j]))) break;
				if (is_different) {
					vid[l] = vcs.size();
					vcs.push_back(v[l]);
				} else
					vid[l] = j;
				(facets++) = vid[l];
			}
			egs.insert(minmax(vid[0], vid[1]));
			egs.insert(minmax(vid[1], vid[2]));
			egs.insert(minmax(vid[2], vid[0]));
		}
	}

	fclose(fp);

	for (vector<Vrtx>::iterator it = vcs.begin(), end = vcs.end(); it != end; ++it) {
		(vertices++) = (*it)[0];
		(vertices++) = (*it)[1];
		(vertices++) = (*it)[2];
	}
	for (std::set<pair<int, int>>::iterator it = egs.begin(), end = egs.end(); it != end; ++it) {
		(edges++) = it->first;
		(edges++) = it->second;
	}
	return true;
}

} // namespace yade
