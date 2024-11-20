/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2016.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#pragma once

#include <cassert>
#include <cstdarg>
#include <stdexcept>
#include <map>
#include <algorithm>
#include <iomanip>
#include <memory>
#include <cfloat>
#include <netcdf>

using namespace netCDF;
using namespace netCDF::exceptions;

#include "vector_utils.h"
#include "string_utils.h"
#include "logger.h"
#include "file_formats.h"

// GDAL functionality must be explicitly enabled
#ifdef ENABLE_GDAL
#include "crs.h"
#endif

// CGAL functionality must be explicitly enabled
#ifdef ENABLE_CGAL
#include "cgal_utils.h"
#endif

#include "marray.hxx"

namespace GeophysicsNetCDF {

constexpr auto DN_POINT = "point";
constexpr auto DN_LINE = "line";

constexpr auto VN_LI_START = "line_index_start";
constexpr auto VN_LI_COUNT = "line_index_count";
constexpr auto VN_LINE_INDEX = "line_index";

constexpr auto LN_LINE_INDEX = "line_index";
constexpr auto LN_LINE_NUMBER = "line_number";
constexpr auto LN_SAMPLE_NUMBER = "point_number";
constexpr auto LN_FLIGHT_NUMBER = "flight_number";

constexpr auto AN_STANDARD_NAME = "standard_name";
constexpr auto AN_LONG_NAME = "long_name";
constexpr auto AN_UNITS = "units";
constexpr auto AN_DESCRIPTION = "description";
constexpr auto AN_MISSINGVALUE = "_FillValue";
constexpr auto AN_ORIGINAL_DATASET_NAME = "original_dataset_name";
constexpr auto AN_ORIGINAL_DATASET_FIELDNAME = "original_dataset_fieldname";

inline NcType getnctype(const uint8_t) { return ncUbyte; }
inline NcType getnctype(const int8_t) { return ncByte; }
inline NcType getnctype(const short) { return ncShort; }
inline NcType getnctype(const int) { return ncInt; }
inline NcType getnctype(const unsigned int) { return ncUint; }
inline NcType getnctype(const float) { return ncFloat; }
inline NcType getnctype(const double) { return ncDouble; }
inline NcType getnctype(const std::string) { return ncString; }

inline uint8_t defaultmissingvalue(const NcUbyte&) { return static_cast<uint8_t>NC_FILL_UBYTE; }
inline int8_t defaultmissingvalue(const NcByte&) { return static_cast<int8_t>NC_FILL_BYTE; }
inline short defaultmissingvalue(const NcShort&) { return static_cast<short>NC_FILL_SHORT; }
inline int defaultmissingvalue(const NcInt&) { return static_cast<int>NC_FILL_INT; }
inline unsigned int defaultmissingvalue(const NcUint&) { return static_cast<unsigned int>NC_FILL_UINT; }
inline float defaultmissingvalue(const NcFloat&) { return static_cast<float>NC_FILL_FLOAT; }
inline double defaultmissingvalue(const NcDouble&) { return static_cast<double>NC_FILL_DOUBLE; }
inline std::string defaultmissingvalue(const NcString&) { return std::string(NC_FILL_STRING); }

class cExportFormat {

public:
	char   form = '\0';
	int    width = 0;
	int    decimals = 0;
	double nullvalue = 0;

	cExportFormat() {};

	cExportFormat(char _form, int _width, int _decimals, double _nullvalue) {
		form = _form;
		width = _width;
		decimals = _decimals;
		nullvalue = _nullvalue;
	}
};

class GFile;

class GVar : public NcVar {

private:
	const GFile& Parent;

public:

	//Do not allow implicit definition of copy constructor
	//cGeophysicsVar(const cGeophysicsVar& rhs) = delete;
	
	//Do not allow implicit definition of assignment operator
	GVar& operator=(const GVar& rhs) = delete;
	GVar& operator=(GVar& rhs) = delete;

	GVar(const GFile& parent, const NcVar& var) 
		: NcVar(var), Parent(parent)
	{};

	size_t line_index_start(const size_t& index) const;
	
	size_t line_index_count(const size_t& index) const;

	size_t length() {
		std::vector<NcDim> dims = getDims();
		size_t len = dims[0].getSize();
		for (size_t di = 1; di < dims.size(); di++) {
			len *= dims[di].getSize();
		}
		return len;
	}

	size_t sizeBytes() {
		return length() * getType().getSize();
	}

	size_t elementspersample() const {
		std::vector<NcDim> dims = getDims();
		size_t len = 1;
		if (dims.size() == 1) {
			return len;
		}
		else {
			len = dims[1].getSize();
			for (size_t di = 2; di < dims.size(); di++) {
				len *= dims[di].getSize();
			}
		}
		return len;
	}

	size_t lineelements(const size_t& lineindex) const {
		return elementspersample() * line_index_count(lineindex);
	}

	size_t nbands() const {
		std::vector<NcDim> dims = getDims();
		if (dims.size() == 1) return 1;
		else if (dims.size() == 2) return dims[1].getSize();
		return 0;
		//this need fixing for case of extra dims
	}

	bool isLineVar() const {
		if (getDimCount() == 0) return false;
		if (getDim(0).getName() == DN_LINE) return true;
		return false;
	}

	bool isSampleVar() const {
		if (getDimCount() == 0) return false;
		if (getDim(0).getName() == DN_POINT) return true;
		return false;
	}

	NcVarAtt add_attribute(const std::string& att, std::string value) {
		return putAtt(att, value);
	}

	NcVarAtt add_standard_name(const std::string& value) {
		return putAtt(AN_STANDARD_NAME, value);
	}

	NcVarAtt add_long_name(const std::string& value) {
		return putAtt(AN_LONG_NAME, value);
	}

	NcVarAtt add_original_dataset_fieldname(const std::string& value) {
		return putAtt(AN_ORIGINAL_DATASET_FIELDNAME, value);
	}

	NcVarAtt add_units(const std::string& value) {
		return putAtt(AN_UNITS, value);
	}

	NcVarAtt add_description(const std::string& value) {
		return putAtt(AN_DESCRIPTION, value);
	}

	static bool hasAtt(const NcVar& var, const std::string& name) {
		nc_type vr_type;
		size_t  vr_len;
		int status = nc_inq_att(var.getParentGroup().getId(), var.getId(), name.c_str(), &vr_type, &vr_len);
		if (status == NC_NOERR) return true;
		else return false;
	};

	bool hasAtt(const std::string& name) const {
		return GVar::hasAtt(*this, name);
	}

	std::string getStringAtt(const std::string& attname) const {
		std::string attvalue;
		if (hasAtt(attname)) {
			NcVarAtt a = getAtt(attname);
			a.getValues(attvalue);
		}
		return attvalue;
	}

	std::string getUnits() const {
		return getStringAtt(AN_UNITS);
	}

	std::string getDescription() const {
		return getStringAtt(AN_DESCRIPTION);
	}

	double lowest_possible_value() const {
		const NcType type = getType();
		if (type == ncShort) return std::numeric_limits<short>::lowest();
		else if (type == ncUint) return std::numeric_limits<unsigned int>::lowest();
		else if (type == ncInt) return (double)std::numeric_limits<int>::lowest();
		else if (type == ncFloat) return std::numeric_limits<float>::lowest();
		else return std::numeric_limits<double>::lowest();
	}

	double highest_possible_value() {
		const NcType type = getType();
		if (type == ncShort) return std::numeric_limits<short>::max();
		else if (type == ncUint) return std::numeric_limits<unsigned int>::max();
		else if (type == ncInt) return std::numeric_limits<int>::max();
		else if (type == ncFloat) return std::numeric_limits<float>::max();
		else return std::numeric_limits<double>::max();
	}

	void set_default_missingvalue() {
		nc_type t = getType().getId();
		switch (t) {
		case NC_UBYTE: putAtt(AN_MISSINGVALUE, ncUbyte, defaultmissingvalue(ncUbyte)); break;
		case NC_BYTE: putAtt(AN_MISSINGVALUE, ncByte, defaultmissingvalue(ncByte)); break;
		case NC_SHORT: putAtt(AN_MISSINGVALUE, ncShort, defaultmissingvalue(ncShort)); break;
		case NC_INT: putAtt(AN_MISSINGVALUE, ncInt, defaultmissingvalue(ncInt)); break;
		case NC_UINT: putAtt(AN_MISSINGVALUE, ncUint, defaultmissingvalue(ncUint)); break;
		case NC_FLOAT: putAtt(AN_MISSINGVALUE, ncFloat, defaultmissingvalue(ncFloat)); break;
		case NC_DOUBLE: putAtt(AN_MISSINGVALUE, ncDouble, defaultmissingvalue(ncDouble)); break;
		case NC_STRING: {
			//std::string n = defaultmissingvalue(ncString);
			//putAtt(AN_MISSINGVALUE, ncString, n.length(), n.c_str()); break;
			break;
		}
		default: {
			std::string msg = _SRC_ + strprint("\nAttempt to set default missing value of unsupported datatype\n");
			throw(std::exception(msg.c_str()));
		}
		}
		return;
	}

	template<typename T>
	bool add_missing_value(const T& value) {
		const NcType type = getType();
		if (type == ncShort) putAtt(AN_MISSINGVALUE, ncShort, value);
		else if (type == ncInt) putAtt(AN_MISSINGVALUE, ncInt, value);
		else if (type == ncFloat) putAtt(AN_MISSINGVALUE, ncFloat, value);
		else if (type == ncDouble) putAtt(AN_MISSINGVALUE, ncDouble, value);
		else return false;
		return true;
	}

	template<typename T>
	T missingvalue(const T&) const {
		if (hasAtt(AN_MISSINGVALUE)) {
			T v;
			getAtt(AN_MISSINGVALUE).getValues(&v);
			return v;
		}
		else {
			const NcType type = getType();
			if (type == ncShort) return (T)NC_FILL_SHORT;
			else if (type == ncInt)  return (T)NC_FILL_INT;
			else if (type == ncFloat) return (T)NC_FILL_FLOAT;
			else if (type == ncDouble) return (T)NC_FILL_DOUBLE;
			else return (T)NC_FILL_SHORT;
		}
	}

	template<typename T>
	bool getAll(std::vector<T>& vals) {
		vals.resize(length());
		getVar(vals.data());
		return true;
	}

	template<typename T>
	bool minmax(T& minval, T& maxval) {
		std::vector<T> vals;
		getAll(vals);
		T nullv;
		nullv = missingvalue(nullv);
		minval = highest_possible_value();
		maxval = lowest_possible_value();
		for (size_t i = 0; i < vals.size(); i++) {
			if (vals[i] == nullv) continue;
			if (vals[i] < minval) minval = vals[i];
			if (vals[i] > maxval) maxval = vals[i];
		}
		return true;
	}

	template<typename T>
	void getLine(const size_t& lineindex, andres::Marray<T>& A) const {
		if (isNull()) {
			std::string msg = _SRC_ + strprint("\nAttempt to read from a Null variable\n");
			throw(std::exception(msg.c_str()));
		}
		std::vector<NcDim>  dims = getDims();
		std::vector<size_t> start((size_t)getDimCount());
		std::vector<size_t> count((size_t)getDimCount());


		if (isLineVar()) {
			start[0] = lineindex;
			count[0] = 1;
		}
		else {
			start[0] = line_index_start(lineindex);
			count[0] = line_index_count(lineindex);
		}

		for (size_t i = 1; i < dims.size(); i++) {
			start[i] = 0;
			count[i] = dims[i].getSize();
		}

		A.resize(count.data(), count.data() + count.size());
		getVar(start, count, &(A(0)));
	}

	template<typename T>
	bool getRecord(const size_t& record, std::vector<T>& v) const
	{
		if (isNull()) {
			std::string msg = _SRC_ + strprint("\nAttempt to read from a Null variable\n");
			throw(std::exception(msg.c_str()));
			return false;
		}
		std::vector<size_t> startp = { record, 0 };
		std::vector<size_t> countp = { 1,      nbands() };
		v.resize(nbands());
		getVar(startp, countp, v.data());
		return true;
	};

	template<typename T>
	bool putRecord(const size_t& record, const T& v) const
	{
		if (isNull()) {
			std::string msg = _SRC_ + strprint("\nAttempt to write to a Null variable\n");
			throw(std::exception(msg.c_str()));
			return false;
		}
		std::vector<size_t> startp = { record, 0 };
		std::vector<size_t> countp = { 1,      1 };
		putVar(startp, countp, &v);
		return true;
	};

	template<typename T>
	bool putRecord(const size_t& record, const std::vector<T>& v) const
	{
		if (isNull()) {
			std::string msg = _SRC_ + strprint("\nAttempt to write to a Null variable\n");
			throw(std::exception(msg.c_str()));
			return false;
		}
		std::vector<size_t> startp = { record, 0 };
		std::vector<size_t> countp = { 1,      nbands() };
		assert(countp[1] == v.size());
		putVar(startp, countp, v.data());
		return true;
	};

	bool donotexport() {
		std::vector<std::string> s = {
			DN_POINT, VN_LI_START, VN_LI_COUNT,
			"longitude_first", "longitude_last",
			"latitude_first", "latitude_last",
			"easting_first", "easting_last",
			"northing_first", "northing_last" };

		auto it = std::find(s.begin(), s.end(), getName());
		if (it == s.end()) return false;
		else return true;
	};

	cExportFormat defaultexportformat() const {
		cExportFormat e;
		nc_type t = getType().getId();
		switch (t) {
		case NC_SHORT: e = cExportFormat('I', 8, 0, -999); break;
		case NC_INT: e = cExportFormat('I', 12, 0, -999); break;
		case NC_UINT: e = cExportFormat('I', 12, 0, 999); break;
		case NC_FLOAT: e = cExportFormat('F', 10, 4, -999); break;
		case NC_DOUBLE: e = cExportFormat('F', 16, 6, -999); break;
		default: e = cExportFormat('F', 16, 6, -999);  break;
		}
		return e;
	}
};

class GLineVar : public GVar {

private:

public:

	GLineVar(const GFile& parent, const NcVar& var)
		: GVar(parent, var)
	{}

	template<typename T>
	bool putAll(std::vector<T> vals) {
		if (isNull()) {
			std::string msg = _SRC_ + strprint("\nAttempt to write to a Null variable\n");
			throw(std::exception(msg.c_str()));
		}

		if (vals.size() != length()) {
			std::string msg = _SRC_ + strprint("\nAttempt to write variable (%s) with non-matching size\n", getName().c_str());
			throw(std::exception(msg.c_str()));
		}

		putVar(vals.data());
		return true;
	}

	template<typename T>
	bool getLine(const size_t& lineindex, const size_t& bandindex, T& val) {
		if (isNull()) { return false; }
		std::vector<size_t> startp = { lineindex, bandindex };
		std::vector<size_t> countp = { 1, 1 };
		getVar(startp, countp, &val);
		return true;
	}

	template<typename T>
	bool getLine(const size_t& lineindex, T& val) const {
		return getLine(lineindex, 0, val);
	}

	template<typename T>
	bool getLine(const size_t& lineindex, const size_t& bandindex, std::vector<T>& vals) const {
		if (isNull()) { return false; }
		std::vector<size_t> startp = { lineindex, bandindex };
		std::vector<size_t> countp = { 1, 1 };
		vals.resize(countp[0] * countp[1]);
		getVar(startp, countp, vals.data());
		return true;
	}

	template<typename T>
	bool getLine(const size_t& lineindex, std::vector<T>& vals) const {
		return getLine(lineindex, 0, vals);
	}

	template<typename T>
	T getSample(const size_t& lineindex, const size_t& sampleindex, const size_t& bandindex, T& val) const {
		if (isNull()) { return false; }
		std::vector<size_t> startp = { lineindex, bandindex };
		std::vector<size_t> countp = { 1, 1 };
		getVar(startp, countp, &val);
		return val;
	};

};

class GSampleVar : public GVar {

private:

public:

	GSampleVar(const GFile& parent, const NcVar& var)
		: GVar(parent, var)
	{}

	template<typename T>
	bool putAll(const std::vector<T>& vals) {
		if (isNull()) {
			std::string msg = _SRC_ + strprint("\nAttempt to write to a Null variable\n");
			throw(std::exception(msg.c_str()));
		}

		if (vals.size() != length()) {
			std::string msg = _SRC_ + strprint("\nAttempt to write variable (%s) with non-matching size\n", getName().c_str());
			throw(std::exception(msg.c_str()));
		}

		putVar(vals.data());
		return true;
	}

	template<typename T>
	bool putLineBand(const size_t& lineindex, const size_t& bandindex, const std::vector<T>& vals) {
		if (isNull()) {
			std::string msg = _SRC_ + strprint("\nAttempt to write to a Null variable\n");
			throw(std::exception(msg.c_str()));
		}

		std::vector<NcDim>  dims = getDims();
		if (dims.size() > 2) {
			std::string msg = _SRC_ + strprint("\nAttempt to use putLineBand() to write to a variable with more than 2 dimensions\n");
			throw(std::exception(msg.c_str()));
		}

		if (vals.size() != line_index_count(lineindex)) {
			std::string msg = _SRC_ + strprint("\nAttempt to write line/band of variable (%s) with non-matching size\n", getName().c_str());
			throw(std::exception(msg.c_str()));
		}

		std::vector<size_t> start(getDimCount());
		std::vector<size_t> count(getDimCount());
		start[0] = line_index_start(lineindex);
		count[0] = line_index_count(lineindex);
		start[1] = bandindex;
		count[1] = 1;
		putVar(start, count, vals.data());
		return true;
	}

	template<typename T>
	bool putLine(const size_t& lineindex, const std::vector<T>& vals) {
		if (isNull()) {
			std::string msg = _SRC_ + strprint("\nAttempt to write to a Null variable\n");
			throw(std::exception(msg.c_str()));
		}

		size_t sz = lineelements(lineindex);
		if (vals.size() != sz) {
			std::string msg = _SRC_ + strprint("\nAttempt to write line/band of variable (%s) with non-matching size\n", getName().c_str());
			throw(std::exception(msg.c_str()));
		}

		std::vector<NcDim>  dims = getDims();
		std::vector<size_t> start(getDimCount());
		std::vector<size_t> count(getDimCount());
		start[0] = line_index_start(lineindex);
		count[0] = line_index_count(lineindex);
		for (size_t i = 1; i < dims.size(); i++) {
			start[i] = 0;
			count[i] = dims[i].getSize();
		}
		putVar(start, count, vals.data());
		return true;
	}

	template<typename T>
	bool getLine(const size_t& lineindex, std::vector<T>& vals) {
		if (isNull()) {
			std::string msg = _SRC_ + strprint("\nAttempt to read from a Null variable\n");
			throw(std::exception(msg.c_str()));
		}

		std::vector<NcDim>  dims = getDims();
		std::vector<size_t> start((size_t)getDimCount());
		std::vector<size_t> count((size_t)getDimCount());
		start[0] = line_index_start(lineindex);
		count[0] = line_index_count(lineindex);
		for (size_t i = 1; i < dims.size(); i++) {
			start[i] = 0;
			count[i] = dims[i].getSize();
		}
		size_t sz = lineelements(lineindex);
		vals.resize(sz);
		getVar(start, count, vals.data());
		return true;
	}

	template<typename T>
	void getLine_temp(const size_t& lineindex, andres::Marray<T>& A) const {
		if (isNull()) {
			std::string msg = _SRC_ + strprint("\nAttempt to read from a Null variable\n");
			throw(std::exception(msg.c_str()));
		}
		std::vector<NcDim>  dims = getDims();
		std::vector<size_t> start((size_t)getDimCount());
		std::vector<size_t> count((size_t)getDimCount());
		start[0] = line_index_start(lineindex);
		count[0] = line_index_count(lineindex);
		for (size_t i = 1; i < dims.size(); i++) {
			start[i] = 0;
			count[i] = dims[i].getSize();
		}
		A.resize(count.data(), count.data() + count.size());
		size_t sz = lineelements(lineindex);
		getVar(start, count, &(A(0)));
	}

	template<typename T>
	bool getSample(const size_t& lineindex, const size_t& sampleindex, const size_t& bandindex, T& val) const {
		if (isNull()) { return false; }
		std::vector<size_t> startp = { line_index_start(lineindex) + sampleindex, bandindex };
		std::vector<size_t> countp = { 1, 1 };
		getVar(startp, countp, &val);
		return true;
	};

	template<typename T>
	bool getPoint(const size_t& pointindex, std::vector<T>& v) const
	{
		if (isNull()) return false;
		getRecord(pointindex, v);
		return true;
	};

	template<typename T>
	bool putPoint(const size_t& pointindex, const std::vector<T>& v) const
	{
		if (isNull()) return false;
		putRecord(pointindex, v);
		return true;
	};

};

class GFile : public NcFile {

private:

	std::vector<unsigned int> line_index_start;
	std::vector<unsigned int> line_index_count;
	std::vector<unsigned int> line_number;

	NcDim dim_sample() { return getDim(DN_POINT); }

	NcDim dim_line() { return getDim(DN_LINE); }

	bool InitialiseExisting() {
		if (readLineIndex() == false) return false;
		if (getLineNumbers(line_number) == false) return false;
		return true;
	}

	bool readLineIndex() {
		if (hasVar(VN_LI_COUNT)) {
			GLineVar vc = getLineVar(VN_LI_COUNT);
			if (vc.getAll(line_index_count) == false)return false;

			GLineVar vs = getLineVar(VN_LI_START);
			if (vs.getAll(line_index_start) == false)return false;
		}
		else if (hasVar(VN_LINE_INDEX)) {
			GLineVar vl = getLineVar(VN_LINE_INDEX);
			std::vector<unsigned int> line_index;
			if (vl.getAll(line_index) == false)return false;
			set_start_count(line_index);
		}

		GLineVar v = getLineVar(DN_LINE);
		v.getAll(line_number);
		return true;
	}

	static std::vector<unsigned int> compute_line_index(std::vector<unsigned int>& count)
	{
		std::vector<unsigned int> index(sum(count));
		size_t k = 0;
		for (size_t i = 0; i < count.size(); i++) {
			for (size_t j = 0; j < count[i]; j++) {
				index[k] = (unsigned int)i;
				k++;
			}
		}
		return index;
	}

	void set_start_count(const std::vector<unsigned int>& line_index)
	{
		line_index_start = compute_start_from_line_index(line_index);
		line_index_count = compute_count_from_start(line_index_start, line_index.size());
		return;
	}

	static std::vector<unsigned int> compute_start_from_line_index(const std::vector<unsigned int>& line_index)
	{
		size_t ns = line_index.size();
		std::vector<unsigned int> start;
		start.push_back(0);
		for (size_t i = 1; i < ns; i++) {
			if (line_index[i] != line_index[i - 1]) {
				start.push_back((unsigned int)i);
			}
		}
		return start;
	}

	static std::vector<unsigned int> compute_start_from_count(const std::vector<unsigned int>& count)
	{
		const size_t nl = count.size();
		std::vector<unsigned int> start(nl);
		start[0] = 0;
		for (size_t i = 1; i < nl; i++) {
			start[i] = start[i - 1] + count[i - 1];
		}
		return start;
	}

	static std::vector<unsigned int> compute_count_from_start(const std::vector<unsigned int>& start, size_t nsamplestotal)
	{
		size_t nl = start.size();
		std::vector<unsigned int> count(nl);
		for (size_t i = 0; i < nl - 1; i++) {
			count[i] = start[i + 1] - start[i];
		}
		count[nl - 1] = (unsigned int)(nsamplestotal - start[nl - 1]);
		return count;
	}

	static size_t nelements(const NcVar& v)
	{
		std::vector<NcDim> dims = v.getDims();
		size_t n = 0;
		for (size_t di = 0; di < dims.size(); di++) {
			if (di == 0) n = 1;
			n *= dims[di].getSize();
		}
		return n;
	}

	static size_t nbytes(const NcVar& v)
	{
		size_t n = nelements(v);
		n *= v.getType().getSize();
		return n;
	}

	static size_t nbytes(const NcAtt& a)
	{
		size_t n = a.getAttLength() * a.getType().getSize();
		return n;
	}



public:

	//Do not allow implicit definition of copy constructor or assignment operators	
	GFile& operator=(const NcGroup& rhs) = delete;
	GFile& operator=(const NcFile& rhs) = delete;
	GFile& operator=(const GFile& rhs) = delete;

	GFile(const NcGroup& rhs) = delete;
	GFile(const NcFile& rhs) = delete;
	GFile(const GFile& rhs) = delete;

	// Constructor generates a null object.
	GFile() : NcFile() {} // invoke base class constructor	

	//Open existing file constructor
	GFile(const std::string& ncpath, const netCDF::NcFile::FileMode& filemode = netCDF::NcFile::FileMode::read)
		: netCDF::NcFile(ncpath, filemode)
	{
		open(ncpath, filemode);
	};

	//Destructor
	~GFile() {};

	size_t get_line_index_start(const size_t& li) const { return line_index_start[li]; }
	size_t get_line_index_count(const size_t& li) const { return line_index_count[li]; }

	std::string pathname() const {
		size_t len;
		nc_inq_path(getId(), &len, nullptr);
		char* p = new char[len + 1];
		nc_inq_path(getId(), nullptr, p);
		std::string path(p);
		delete[] p;
		return path;
	}

	bool isopen() const {
		if (line_index_start.size() > 0)	return true;
		return false;
	}

	std::vector<unsigned int> get_line_index() const {
		std::vector<unsigned int> index(ntotalsamples());
		getVar(VN_LINE_INDEX).getVar(index.data());
		return index;
	}

	void open(const std::string& ncpath, const FileMode& filemode = NcFile::FileMode::read)
	{
		if (filemode == NcFile::read) {
			NcFile::open(ncpath, filemode);
			InitialiseExisting();
		}
		else if (filemode == NcFile::write) {
			NcFile::open(ncpath, filemode);
			InitialiseExisting();
		}
		else if (filemode == NcFile::replace) {

		}
		else if (filemode == NcFile::newFile) {

		}
		else {

		}
	}

	bool InitialiseNew(const std::vector<size_t>& linenumbers, const std::vector<size_t>& linesamplecount) {
		size_t n = linenumbers.size();
		std::vector<unsigned int> uint_linenumbers(n);
		std::vector<unsigned int> uint_linesamplecount(n);
		for (size_t i = 0; i < n; i++) {
			uint_linenumbers[i] = (unsigned int)linenumbers[i];
			uint_linesamplecount[i] = (unsigned int)linesamplecount[i];
		}
		return InitialiseNew(uint_linenumbers, uint_linesamplecount);
	}

	bool InitialiseNew(const std::vector<unsigned int>& linenumbers, const std::vector<unsigned int>& linesamplecount) {
		const size_t nl = linenumbers.size();
		line_number = linenumbers;
		line_index_count = linesamplecount;
		line_index_start = compute_start_from_count(line_index_count);

		size_t nsamples = sum(line_index_count);
		NcDim ds = addDim(DN_POINT, nsamples);
		NcDim dl = addDim(DN_LINE, nl);

		std::vector<unsigned int> line_index = compute_line_index(line_index_count);
		bool status = add_line_index(line_index);
		if (status == false) {
			std::string msg = _SRC_ + strprint("\nCould no add line_index variable (%s)\n", VN_LINE_INDEX);
			throw(std::exception(msg.c_str()));
		}

		status = add_line_number_variable(linenumbers);
		if (status == false) {
			std::string msg = _SRC_ + strprint("\nCould no add line_number variable (%s)\n", DN_LINE);
			throw(std::exception(msg.c_str()));
		}

		return true;
	}

	bool add_line_index(const std::vector<unsigned int> line_index)
	{
		bool status = addSampleVar(VN_LINE_INDEX, ncUint);
		if (status) {
			GSampleVar v = getSampleVar(VN_LINE_INDEX);
			v.putVar(line_index.data());
			v.add_long_name(LN_LINE_INDEX);
			v.add_description("zero-based index of line associated with point");
			//v.add_units("1");
		}
		return status;
	}

	//Deprecated
	bool add_line_index_start(const std::vector<unsigned int> line_index_start)
	{
		bool status = addLineVar(VN_LI_START, ncUint);
		if (status) {
			GLineVar vstart = getLineVar(VN_LI_START);
			vstart.putVar(line_index_start.data());
			vstart.add_long_name(VN_LI_START);
			vstart.add_description("zero-based index of the first sample in the line");
			//vstart.add_units("1");
		}
		return status;
	}

	//Deprecated
	bool add_line_index_count(const std::vector<unsigned int> line_index_count)
	{
		bool status = addLineVar(VN_LI_COUNT, ncUint);
		if (status) {
			GLineVar vcount = getLineVar(VN_LI_COUNT);
			vcount.putVar(line_index_count.data());
			vcount.add_long_name(VN_LI_COUNT);
			vcount.add_description("number of samples in the line");
			//vcount.add_units("1");
		}
		return status;
	}

	//Deprecated
	bool add_point_variable()
	{
		bool status = addSampleVar(DN_POINT, ncUint);
		if (status) {
			//std::vector<unsigned int> sample = increment((unsigned int)ntotalsamples(), (unsigned int)0, (unsigned int)1);
			std::vector<unsigned int> sample = increment<unsigned int>(ntotalsamples(), 0, 1);
			GSampleVar v = getSampleVar(DN_POINT);
			v.putVar(sample.data());
			v.add_long_name(LN_SAMPLE_NUMBER);
			v.add_description("sequential point number");
			//v.add_units("1");
		}
		return status;
	}

	bool add_line_number_variable(const std::vector<unsigned int> linenumbers)
	{
		bool status = addLineVar(DN_LINE, ncUint);
		if (status) {
			GLineVar vline = getLineVar(DN_LINE);
			vline.putVar(linenumbers.data());
			vline.add_long_name(LN_LINE_NUMBER);
			vline.add_description("flight line number");
			//vline.add_units("1");
		}
		return status;
	}

	template<typename T>
	bool isinlist(const std::vector<T>& vec, const T& val) {
		if (std::find(vec.begin(), vec.end(), val) == vec.end()) return false;
		return true;
	}

	bool subsample(const GFile& srcfile, const size_t subsamplerate, std::vector<std::string> include_varnames, std::vector<std::string> exclude_varnames)
	{
		unsigned int nsnew = (unsigned int)std::ceil((double)srcfile.ntotalsamples() / (double)subsamplerate);

		const auto srclineindex = srcfile.get_line_index();
		std::vector<unsigned int> ss_srclineindex(nsnew);
		size_t k = 0;
		for (size_t pi = 0; pi < nsnew; pi++) {
			ss_srclineindex[pi] = srclineindex[k];
			k += subsamplerate;
		}
		auto start = compute_start_from_line_index(ss_srclineindex);
		auto count = compute_count_from_start(start, nsnew);

		std::vector<unsigned int> linenumber(start.size());
		for (size_t li = 0; li < start.size(); li++) {
			linenumber[li] = srcfile.line_number[ss_srclineindex[start[li]]];
		}

		InitialiseNew(linenumber, count);

		auto dm = srcfile.getDims();
		for (auto dit = dm.begin(); dit != dm.end(); dit++) {
			NcDim& srcdim = dit->second;
			const std::string dimname = srcdim.getName();
			if (hasDim(dimname)) continue;
			else {
				addDim(dimname, srcdim.getSize());
			}
		}

		copy_global_atts(srcfile);
		auto vm = srcfile.getVars();
		for (auto vit = vm.begin(); vit != vm.end(); vit++) {
			NcVar& srcvar = vit->second;
			const std::string vname = srcvar.getName();
			if (vname == VN_LI_START) continue;
			if (vname == VN_LI_COUNT) continue;
			if (vname == DN_POINT) continue;

			bool incstatus = isinlist(include_varnames, vname);
			bool excstatus = isinlist(exclude_varnames, vname);

			bool status = true;
			if (include_varnames.size() > 0 && incstatus == false) {
				status = false;
			}

			if (exclude_varnames.size() > 0 && excstatus == true) {
				status = false;
			}

			if (status) {
				copy_var(subsamplerate, srcvar);
			}
		}
		return true;
	}

#ifdef ENABLE_GDAL
	//Convert legacy file
	bool convert_legacy(const cGeophysicsNcFile& srcfile)
	{
		InitialiseNew(srcfile.line_number, srcfile.line_index_count);
		copy_global_atts(srcfile);
		copy_dims(srcfile);
		auto vm = srcfile.getVars();
		for (auto vit = vm.begin(); vit != vm.end(); vit++) {
			NcVar& srcvar = vit->second;
			if (srcvar.getName() == VN_LI_START) continue;
			if (srcvar.getName() == VN_LI_COUNT) continue;
			if (srcvar.getName() == DN_POINT) continue;

			if (srcvar.getName() == "crs") {
				NcVarAtt a = srcvar.getAtt("epsg_code");
				std::string epsg_string;
				a.getValues(epsg_string);
				int epsgcode = cCRS::epsgcode(epsg_string);
				addCRS(epsgcode);
			}
			else {
				copy_var(1, srcvar);
				NcVar v = getVar(srcvar.getName());
				nc_rename_att(getId(), v.getId(), "standard_name", "long_name");

				if (v.getName() == "latitude") {
					v.putAtt("standard_name", "latitude");
				}

				if (v.getName() == "longitude") {
					v.putAtt("standard_name", "longitude");
				}

				//Remove units = "1"				
				if (cGeophysicsVar::hasAtt(v, AN_UNITS)) {
					NcVarAtt a = v.getAtt(AN_UNITS);
					if (!a.isNull()) {
						std::string units;
						a.getValues(units);
						if (units == "1") {
							nc_del_att(getId(), v.getId(), AN_UNITS);
						}
					}
				}
			}
		}
		return true;
	}
#endif

	bool copy_global_atts(const GFile& srcfile)
	{
		auto m = srcfile.getAtts();
		for (auto it = m.begin(); it != m.end(); it++) {
			NcGroupAtt& srcatt = it->second;
			size_t len = nbytes(srcatt);
			if (len > 0) {
				const std::vector<uint8_t> buf(len);
				srcatt.getValues((void*)buf.data());
				putAtt(srcatt.getName(), srcatt.getType(), srcatt.getAttLength(), (void*)buf.data());
			}

		}
		return true;
	}

	bool copy_dims(const GFile& srcfile)
	{
		auto dm = srcfile.getDims();
		for (auto dit = dm.begin(); dit != dm.end(); dit++) {
			NcDim& srcdim = dit->second;
			if (hasDim(srcdim.getName())) continue;
			addDim(srcdim.getName(), srcdim.getSize());
		}
		return true;
	}


	bool copy_var(const size_t& subsample, const NcVar& srcvar, std::string newname = std::string())
	{
		std::string name = newname;
		if (name.size() == 0) name = srcvar.getName();
		if (hasVar(name)) return false;

		std::vector<NcDim> srcdims = srcvar.getDims();
		size_t nd = srcdims.size();
		if (nd == 0) return true;

		std::vector<size_t> start(nd);
		std::vector<size_t> count(nd);
		std::vector<ptrdiff_t> stride(nd);
		std::vector<ptrdiff_t> stride_out(nd);

		std::vector<NcDim> dims = srcvar.getDims();
		size_t nelements = 1;
		for (size_t i = 0; i < srcdims.size(); i++) {
			start[i] = 0;
			count[i] = srcdims[i].getSize();
			stride[i] = 1;
			stride_out[i] = 1;
			if (srcdims[i].getName() == DN_POINT) {
				count[i] = (size_t)std::ceil(srcdims[i].getSize() / (double)subsample);
				stride[i] = subsample;
				dims[i] = getDim(DN_POINT);
			}
			nelements *= count[i];
		}

		NcVar v = addVar(name, srcvar.getType(), srcvar.getDims());
		v.setCompression(true, true, 9);
		copy_varatts(srcvar, v);

		size_t elementsize = v.getType().getSize();
		size_t bufsize = nelements * elementsize;
		if (bufsize > 0) {
			std::vector<uint8_t> buf(bufsize);
			srcvar.getVar(start, count, stride, (void*)buf.data());
			v.putVar(start, count, stride_out, (void*)buf.data());
			//v.putVar((void*)buf.data());
		}
		return true;
	}

	bool copy_varatts(const NcVar& srcvar, const NcVar& dstvar)
	{
		auto am = srcvar.getAtts();
		for (auto ait = am.begin(); ait != am.end(); ait++) {
			NcVarAtt& srcatt = ait->second;
			size_t len = nbytes(srcatt);
			if (len > 0) {
				const std::vector<uint8_t> buf(len);
				srcatt.getValues((void*)buf.data());
				dstvar.putAtt(srcatt.getName(), srcatt.getType(), srcatt.getAttLength(), (void*)buf.data());
			}
		}
		return true;
	}

	size_t nlines() const { return line_index_start.size(); }
	size_t ntotalsamples() const { return sum(line_index_count); }
	size_t nlinesamples(const size_t lineindex) const { return line_index_count[lineindex]; }

	size_t getLineIndexByPointIndex(const int& pointindex) const
	{
		const size_t n = line_index_start.size() - 1;
		for (size_t k = 0; k < n; k++) {
			if ((int)line_index_start[k + 1] > pointindex) return k;
		}
		if ((int)(line_index_start[n] + line_index_count[n]) > pointindex) return n;
		return undefinedvalue<size_t>();
	};

	size_t getLineIndex(const int& linenumber) {
		auto it = std::find(line_number.begin(), line_number.end(), linenumber);
		return (size_t)(it - line_number.begin());
	}

	bool hasVar(const std::string& varname) {
		NcVar v = getVar(varname);
		if (v.isNull()) return false;
		return true;
	}

	bool hasVarCaseInsensitive(const std::string& invarname) {
		std::string realname;
		NcVar v = getVarByNameCaseInsensitive(invarname,realname);
		if (v.isNull()) return false;
		return true;
	}

	NcVar getVarByNameCaseInsensitive(const std::string& invarname, std::string& varname) {
		std::multimap<std::string, NcVar> vars = getVars();
		for (auto vit = vars.begin(); vit != vars.end(); vit++) {
			if (strcasecmp(invarname, vit->first) == 0) {
				varname = vit->first;
				return vit->second;
			}
		}
		return NcVar();
	}

	bool hasDim(const std::string& dimname) {
		NcDim d = getDim(dimname);
		if (d.isNull()) return false;
		return true;
	}

	NcDim addDim(const std::string& dimname, const size_t& dimsize) {
		if (hasDim(dimname)) {
			NcDim d = getDim(dimname);
			if (d.getSize() == dimsize) return d;
			else {
				std::string msg = _SRC_ + strprint("\nAttempt to add dimension (%s) with size (%zu) unequal to existing dimension (%zu) with same name\n", dimname.c_str(), dimsize, d.getSize());
				throw(std::exception(msg.c_str()));
			}
		}
		else {
			NcDim d = NcGroup::addDim(dimname, dimsize);
			return d;
		}
	}

	NcVar getVarByAtt(const std::string& att_name, const std::string& att_value) {
		std::multimap<std::string, NcVar> vars = getVars();
		for (auto vit = vars.begin(); vit != vars.end(); vit++) {
			std::map<std::string, NcVarAtt> atts = vit->second.getAtts();
			for (auto ait = atts.begin(); ait != atts.end(); ait++) {
				std::string name = ait->second.getName();
				if (name == att_name) {
					std::string value;
					ait->second.getValues(value);
					if (value == att_value) {
						return vit->second;
					}
				}
			}
		}
		return NcVar();
	}

	NcVar getVarByLongName(const std::string& att_value) {
		return getVarByAtt(AN_LONG_NAME, att_value);
	}

	std::string getVarNameByAtt(const std::string& att_name, const std::string& att_value) {
		NcVar var = getVarByAtt(att_name, att_value);
		if (var.isNull()) {
			return std::string();
		}
		else return var.getName();
	}

	std::string getVarNameByLongName(const std::string& att_value) {
		return getVarNameByAtt(AN_LONG_NAME, att_value);
	}

	std::string getVarNameByStandardName(const std::string& att_value) {
		return getVarNameByAtt(AN_STANDARD_NAME, att_value);
	}

	template<typename T>
	bool getLineNumbers(std::vector<T>& vals) {
		NcVar v = getVarByLongName(LN_LINE_NUMBER);
		if (v.isNull() == false) {
			GLineVar var(*this, v);
			return var.getAll(vals);
		}
		return false;
	}

	std::vector<int> getLineNumbers() {
		std::vector<int> num;
		getLineNumbers(num);
		return num;
	}

	template<typename T>
	bool getFlightNumbers(std::vector<T>& vals) {
		NcVar v = getVarByLongName(LN_FLIGHT_NUMBER);
		if (v.isNull() == false) {
			GLineVar var(*this, v);
			return var.getAll(vals);
		}
		return false;
	}

	std::vector<int> getFlightNumbers() {
		std::vector<int> num;
		getFlightNumbers(num);
		return num;
	}

	template<typename T>
	bool getDataByLineIndex(const GSampleVar& var, const size_t& lineindex, std::vector<T>& vals) {
		std::vector<NcDim> dims = var.getDims();
		if (dims.size() != 1) return false;
		std::vector<size_t> start(1);
		std::vector<size_t> count(1);
		start[0] = line_index_start[lineindex];
		count[0] = line_index_count[lineindex];
		vals.resize(count[0]);
		var.getVar(start, count, vals.data());
		return true;
	}

	template<typename T>
	bool getDataByLineNumber(const std::string& varname, const size_t& linenumber, std::vector<T>& vals) {
		size_t index = getLineIndex(linenumber);
		if (index >= nlines())return false;
		return getDataByLineIndex(varname, index, vals);
	}

	template<typename T>
	bool getDataByLineIndex(const std::string& varname, const size_t& lineindex, std::vector<T>& vals) {
		GSampleVar var = getSampleVar(varname);
		return getDataByLineIndex(var, lineindex, vals);
	}

	template<typename T>
	bool getDataByLineIndex(const std::string& varname, const size_t& lineindex, std::vector<std::vector<T>>& vals) {
		NcVar var = NcFile::getVar(varname);
		std::vector<NcDim> dims = var.getDims();
		size_t nd = dims.size();
		if (nd != 2) return false;

		size_t nsamples = line_index_count[lineindex];
		size_t nbands = dims[1].getSize();

		std::vector<size_t> start(2);
		std::vector<size_t> count(2);
		start[0] = line_index_start[lineindex];
		count[0] = nsamples;

		vals.resize(nbands);
		for (size_t bi = 0; bi < nbands; bi++) {
			start[1] = bi;
			count[1] = 1;
			vals[bi].resize(nsamples);
			var.getVar(start, count, vals[bi].data());
		}
		return true;
	}

	template<typename T>
	bool getDataByPointIndex(const std::string& varname, const size_t& pointindex, std::vector<T>& vals) {
		GVar var(*this, getVar(varname));
		if (var.isNull()) {
			std::string msg = _SRC_ + strprint("\nAttempt to read variable (%s)\n", varname.c_str());
			throw(std::exception(msg.c_str()));
		}

		size_t record;
		if (isScalarVar(var)) {
			vals.resize(1);
			NcVar v = getVar(varname);
			v.getVar(vals.data());
			return true;
		}
		else if (isSampleVar(var)) record = pointindex;
		else if (isLineVar(var))   record = getLineIndexByPointIndex((int)pointindex);
		else {
			std::string msg = _SRC_ + strprint("\nVariable %s is neither a \"point\" or a \"line\" variable\n", varname.c_str());
			throw(std::exception(msg.c_str()));
		}
		return var.getRecord(record, vals);
	}

	template<typename T>
	NcDim addDimVar(const std::string& dimname, const std::vector<T>& dimvals) {
		size_t dimsize = dimvals.size();
		NcDim dim = getDim(dimname);
		if (dim.isNull()) {
			dim = addDim(dimname, dimsize);
		}
		else {
			if (dim.getSize() != dimsize) {
				std::string msg = _SRC_ + strprint("\nAttempt to add new dimension (%s) with different size to the existing\n", dimname.c_str());
				throw(std::exception(msg.c_str()));
			}
		}

		NcVar var = getVar(dimname);
		if (var.isNull()) {
			var = addVar(dimname, getnctype(dimvals[0]), dim);
		}
		else {
			if (var.getDim(0).getSize() != dimsize) {
				std::string msg = _SRC_ + strprint("\nAttempt to add new dimension (%s) with different size to the existing\n", dimname.c_str());
				throw(std::exception(msg.c_str()));
			}
		}
		var.putVar(dimvals.data());
		return dim;
	}

	GVar getGeophysicsVar(const std::string& name) {
		if (getVarCount() == 0) return GVar(*this, NcVar());
		return GVar(*this, getVar(name));
	}

	GSampleVar getSampleVar(const std::string& name) {
		if (getVarCount() == 0) return GSampleVar(*this, NcVar());
		return GSampleVar(*this, getVar(name));
	}

	GLineVar getLineVar(const std::string& name) {
		if (getVarCount() == 0) return GLineVar(*this, NcVar());
		return GLineVar(*this, getVar(name));
	}

	bool addSampleVar(const std::string& name, const NcType& type, const std::vector<NcDim>& dims) {
		if (hasVarCaseInsensitive(name)) {
			return false;
		}

		std::vector<NcDim> vardims = { dim_sample() };
		for (size_t i = 0; i < dims.size(); i++) {
			vardims.push_back(dims[i]);
		}
		
		GSampleVar newvar(*this, addVar(name, type, vardims));
		if (newvar.isNull())return false;

		newvar.set_default_missingvalue();
		return true;
	}

	bool addSampleVar(const std::string& name, const NcType& type, const NcDim& banddim = NcDim()) {
		std::vector<NcDim> dims;
		if (banddim.isNull() == false) {
			dims.push_back(banddim);
		}
		return addSampleVar(name, type, dims);
	}

	GSampleVar addgetSampleVar(const std::string& name, const NcType& type, const std::vector<NcDim>& dims) {
		if (addSampleVar(name, type, dims)) {
			return GSampleVar(*this, getVar(name));
		}
		return GSampleVar(*this, NcVar());
	}

	GSampleVar addgetSampleVar(const std::string& name, const NcType& type, const NcDim& banddim = NcDim()) {
		if (addSampleVar(name, type, banddim)) {
			return GSampleVar(*this, getVar(name));
		}
		return GSampleVar(*this, NcVar());
	}

	bool addLineVar(const std::string& name, const NcType& type, const std::vector<NcDim>& dims) {
		if (hasVarCaseInsensitive(name)) {
			return false;
		}

		std::vector<NcDim> vardims = { dim_line() };
		for (size_t i = 0; i < dims.size(); i++) {
			vardims.push_back(dims[i]);
		}

		GLineVar newvar(*this, addVar(name, type, vardims));
		if (newvar.isNull()) return false;

		newvar.set_default_missingvalue();
		return true;
	}

	bool addLineVar(const std::string& name, const NcType& type, const NcDim& banddim = NcDim()) {
		std::vector<NcDim> dims;
		if (banddim.isNull() == false) {
			dims.push_back(banddim);
		}
		return addLineVar(name, type, dims);
	}

	GLineVar addgetLineVar(const std::string& name, const NcType& type, const std::vector<NcDim>& dims) {
		if (addLineVar(name, type, dims)) {
			return GLineVar(*this, getVar(name));
		}
		return GLineVar(*this, NcVar());
	}

	GLineVar addgetLineVar(const std::string& name, const NcType& type, const NcDim& banddim = NcDim()) {
		if (addLineVar(name, type, banddim)) {
			return GLineVar(*this, getVar(name));
		}
		return GLineVar(*this, NcVar());
	}


	bool findNonNullLineStartEndPoints(const std::string& xvar, const std::string& yvar,
		std::vector<double>& x1, std::vector<double>& x2, std::vector<double>& y1, std::vector<double>& y2) {
		x1.resize(nlines());
		x2.resize(nlines());
		y1.resize(nlines());
		y2.resize(nlines());

		GSampleVar vx = getSampleVar(xvar);
		GSampleVar vy = getSampleVar(yvar);
		double nvx = vx.missingvalue(nvx);
		double nvy = vy.missingvalue(nvy);
		for (size_t li = 0; li < nlines(); li++) {
			std::vector<double> x;
			std::vector<double> y;
			vx.getLine(li, x);
			vy.getLine(li, y);

			const size_t ns = x.size();

			x1[li] = nvx;
			y1[li] = nvy;
			for (size_t si = 0; si < ns; si++) {
				if (x[si] != nvx && y[si] != nvy) {
					x1[li] = x[si];
					y1[li] = y[si];
					break;
				}
			}

			x2[li] = nvx;
			y2[li] = nvy;
			for (size_t si = ns - 1; ; si--) {
				if (x[si] != nvx && y[si] != nvy) {
					x2[li] = x[si];
					y2[li] = y[si];
					break;
				}
				if (si == 0)break;//avoid endless loop
			}
		}
		return true;
	}

#ifdef ENABLE_GDAL
	bool addCRS(const int epsgcode) {

#ifdef _WIN32 
		NcVar v = NcFile::addVar("crs", ncByte);
#else 			
		std::vector<NcDim> d;
		NcVar v = NcFile::addVar("crs", ncByte, d);
#endif

		v.putAtt(AN_LONG_NAME, "coordinate_reference_system");
		OGRSpatialReference srs = getsrs(epsgcode);
		char* wkt = NULL;
		srs.exportToWkt(&wkt);
		srs.importFromEPSG(epsgcode);
		v.putAtt("spatial_ref", wkt);
		v.putAtt("inverse_flattening", ncDouble, srs.GetInvFlattening());
		v.putAtt("semi_major_axis", ncDouble, srs.GetSemiMajor());
		v.putAtt("longitude_of_prime_meridian", ncDouble, srs.GetPrimeMeridian());
		return true;
	}
#endif

	bool addLineVariableAndData(
		const NcType& type,
		const std::string name,
		const std::string longname,
		const std::string description,
		const std::string units,
		const std::vector<double> data
	) {
		bool status = addLineVar(name, type);
		if (status) {
			GLineVar v = getLineVar(name);
			v.add_missing_value(defaultmissingvalue(ncDouble));
			v.add_long_name(longname);
			v.add_description(description);
			if(units.size() > 0) v.add_units(units);
			v.putAll(data);
		}
	}

	bool addLineStartEndPointsLL() {
		if (hasVar("longitude_first")) {
			std::string msg = _SRC_ + strprint("\nWarning: Variable longitude_first already exists\n");
			//glog.logmsg(msg);
			return false;
		}

		if (hasVar("latitude_first")) {
			std::string msg = _SRC_ + strprint("\nWarning: Variable latitude_first already exists\n");
			//glog.logmsg(msg);
			return false;
		}

		std::vector<double> x1;
		std::vector<double> x2;
		std::vector<double> y1;
		std::vector<double> y2;

		std::string xvarname = getVarNameByLongName("longitude");
		std::string yvarname = getVarNameByLongName("latitude");
		findNonNullLineStartEndPoints(xvarname, yvarname, x1, x2, y1, y2);

		bool status1 = addLineVariableAndData(ncDouble,"longitude_first","longitude_first","first non-null longitude coordinate in the line","degree_east",x1);
		bool status2 = addLineVariableAndData(ncDouble, "longitude_last", "longitude_last", "last non-null longitude coordinate in the line", "degree_east", x2);
		bool status3 = addLineVariableAndData(ncDouble, "latitude_first", "latitude_first", "first non-null latitude coordinate in the line", "degree_north", y1);
		bool status4 = addLineVariableAndData(ncDouble, "latitude_last", "latitude_last", "last non-null latitude coordinate in the line", "degree_north", y2);
		return true;
	}

	bool addLineStartEndPointsEN() {
		if (hasVar("easting_first")) {
			std::string msg = _SRC_ + strprint("\nWarning: Variable easting_first already exists\n");
			//glog.logmsg(msg);
			return false;
		}

		if (hasVar("northing_first")) {
			std::string msg = _SRC_ + strprint("Warning: Variable northing_first already exists\n");
			//glog.logmsg(msg);
			return false;
		}

		std::vector<double> x1;
		std::vector<double> x2;
		std::vector<double> y1;
		std::vector<double> y2;
		if (hasVar("easting") && hasVar("northing")) {
			findNonNullLineStartEndPoints("easting", "northing", x1, x2, y1, y2);
		}
		else if (hasVar("Easting") && hasVar("Northing")) {
			findNonNullLineStartEndPoints("Easting", "Northing", x1, x2, y1, y2);
		}
		else {
			std::string msg = _SRC_ + strprint("\nWarning: Could not find easting/Easting/northin/Northing\n");
			//glog.logmsg(msg);
			return false;
		}

		bool status1 = addLineVariableAndData(ncDouble, "easting_first", "easting_first", "first non-null easting coordinate in the line", "m", x1);
		bool status2 = addLineVariableAndData(ncDouble, "easting_last", "easting_last", "last non-null easting coordinate in the line", "m", x2);
		bool status3 = addLineVariableAndData(ncDouble, "northing_first", "northing_first", "first non-null northing coordinate in the line", "m", y1);
		bool status4 = addLineVariableAndData(ncDouble, "northing_last", "northing_last", "last non-null northing coordinate in the line", "m", y2);
		return true;
	}

#ifdef ENABLE_CGAL
	bool addAlphaShapePolygon(const std::string xvarname, const std::string yvarname) {
		std::vector<double> x;
		std::vector<double> y;
		std::vector<double> px;
		std::vector<double> py;
		cSampleVar vx = getSampleVar(xvarname);
		cSampleVar vy = getSampleVar(yvarname);
		vx.getAll(x);
		vy.getAll(y);
		double nullx = vx.missingvalue(nullx);
		double nully = vy.missingvalue(nully);

		line_data_alpha_shape_polygon_ch(
			line_index_start, line_index_count,
			x, y, nullx, nully, 64, px, py);

		size_t nv = px.size();
		std::vector<double> poly(nv * 2);
		for (size_t i = 0; i < nv; i++) {
			poly[i * 2] = px[i];
			poly[i * 2 + 1] = py[i];
		}

		std::vector<NcDim> dims;
		dims.push_back(addDim("polygonvertex", nv));
		dims.push_back(addDim("polygonordinate", 2));

		cGeophysicsVar v(*this, addVar("bounding_polygon", ncDouble, dims));
		v.add_long_name("bounding_polygon");
		v.add_description("bounding polygon of survey");
		v.add_units("degree");
		v.putVar(poly.data());
		return true;
	}
#endif

	bool minmax(const std::string& varname, double& minval, double& maxval) {
		GSampleVar var = getSampleVar(varname);
		var.minmax(minval, maxval);
		return true;
	}

	bool addGeospatialMetadataItem(const std::string& varname, const std::string& label, const std::string& units) {

		if (hasVar(varname) == false)return false;

		double vmin, vmax;
		GSampleVar var = getSampleVar(varname);
		var.minmax(vmin, vmax);
		std::string s = "geospatial_" + label;
		putAtt(s + "_min", ncDouble, vmin);
		putAtt(s + "_max", ncDouble, vmax);
		putAtt(s + "_units", units);
		putAtt(s + "_resolution", "point");
		if (label == "vertical") {
			putAtt(s + "_positive", "up");
		}
		return true;
	}

#ifdef ENABLE_GDAL
	bool addGeospatialMetadataXY() {
		bool status;
		if (hasVar("longitude") && hasVar("latitude")) {
			status = addGeospatialMetadataItem("longitude", "lon", "degrees_east");
			if (status == false)return status;

			status = addGeospatialMetadataItem("latitude", "lat", "degrees_north");
			if (status == false)return status;
			addCRS(erm2epsgcode("GDA94"));
		}
		else if (hasVar("easting") && hasVar("northing")) {
			status = addGeospatialMetadataItem("easting", "east", "m");
			if (status == false)return status;

			status = addGeospatialMetadataItem("northing", "north", "m");
			if (status == false)return status;
			addCRS(erm2epsgcode("GDA94"));
		}
		else return false;

		return true;
	}
#endif

	bool addGeospatialMetadataVertical() {

		if (hasVar("height")) {
			bool status = addGeospatialMetadataItem("height", "vertical", "m");
			if (status == false)return status;
		}
		else return false;

		return true;

	}

	bool isScalarVar(const NcVar& var) const {
		if (var.getDimCount() == 0) return true;
		return false;
	}

	bool isLineVar(const NcVar& var) const {
		if (var.getDimCount() == 0) return false;
		if (var.getDim(0).getName() == DN_LINE) {
			return true;
		}
		return false;
	}

	bool isSampleVar(const NcVar& var) const {
		if (var.getDimCount() == 0) return false;
		if (var.getDim(0).getName() == DN_POINT) {
			return true;
		}
		return false;
	}

	std::vector<NcDim> getAllDims() {
		std::vector<NcDim> dims;
		std::multimap<std::string, NcDim> m = getDims();
		for (auto it = m.begin(); it != m.end(); it++) {
			dims.push_back(it->second);
		}
		return dims;
	};

	std::vector<NcVar> getAllVars() {
		std::vector<NcVar> vars;
		std::multimap<std::string, NcVar> vm = getVars();
		for (auto it = vm.begin(); it != vm.end(); it++) {
			vars.push_back(it->second);
		}
		return vars;
	};

	std::vector<GLineVar> getLineVars() {
		std::vector<GLineVar> vars;
		std::multimap<std::string, NcVar> vm = getVars();
		for (auto it = vm.begin(); it != vm.end(); it++) {
			if (isLineVar(it->second)) {
				GLineVar v(*this, it->second);
				vars.push_back(v);
			}
		}
		return vars;
	};

	std::vector<GSampleVar> getSampleVars() {
		std::vector<GSampleVar> vars;
		std::multimap<std::string, NcVar> vm = getVars();
		for (auto it = vm.begin(); it != vm.end(); it++) {
			if (isSampleVar(it->second)) {
				GSampleVar v(*this, it->second);
				vars.push_back(v);
			}
		}
		return vars;
	};

	bool export_ASEGGDF2(const std::string& datfilepath, const std::string& dfnfilepath) {
		std::ofstream of(datfilepath);
		of << std::fixed;

		std::vector<GLineVar>   lvars = getLineVars();
		std::vector<GSampleVar> svars = getSampleVars();
		std::vector<GVar> vars;
		for (size_t i = 0; i < lvars.size(); i++) {
			if (lvars[i].donotexport() == false) {
				vars.push_back(lvars[i]);
			}
		}
		for (size_t i = 0; i < svars.size(); i++) {
			std::string vname = svars[i].getName();
			if (svars[i].donotexport() == false) {
				vars.push_back(svars[i]);
			}
		}

		const size_t nvars = vars.size(); // number of vars to be exported
		std::vector<double>  mval(nvars); // missing value
		std::vector<bool>    islv(nvars); // is it a line var
		std::vector<cExportFormat> efmt(nvars);//format

		cOutputFileInfo I;
		for (size_t vi = 0; vi < nvars; vi++) {
			GVar& v = vars[vi];

			if (v.isLineVar()) islv[vi] = true;
			else islv[vi] = false;
			mval[vi] = v.missingvalue(double(0));
			efmt[vi] = v.defaultexportformat();

			int bands = (int)v.nbands();
			I.addfield(v.getName(), efmt[vi].form, efmt[vi].width, efmt[vi].decimals, bands);

			std::string units = v.getUnits();
			if (units != "1") I.setunits(units);

			std::string desc = v.getDescription();
			desc = v.getStringAtt(AN_LONG_NAME);
			I.setdescription(desc);
		}
		I.write_aseggdf_header(dfnfilepath);

		const size_t nl = nlines();
		for (size_t li = 0; li < nl; li++) {
			std::cout << "Exporting line " << line_number[li] << std::endl;
			const size_t ns = line_index_count[li];

			std::vector<andres::Marray<double>> A(nvars);
			for (size_t vi = 0; vi < nvars; vi++) {
				const GVar& v = vars[vi];
				v.getLine(li, A[vi]);
			}

			for (size_t si = 0; si < ns; si += 100) {
				for (size_t vi = 0; vi < nvars; vi++) {
					const GVar& v = vars[vi];
					of << std::setw(efmt[vi].width);
					of << std::setprecision(efmt[vi].decimals);

					const andres::Marray<double>& a = A[vi];
					double val;
					const size_t nb = v.nbands();

					double* b;
					if (islv[vi]) b = &(a(0));
					else          b = &(a(si));
					for (size_t bi = 0; bi < nb; bi++) {
						val = b[bi];
						if (val == mval[vi]) val = efmt[vi].nullvalue;
						of << val;
					}
				}
				of << std::endl;
			}
		}
		return true;
	};

};

// Defined here only because it needs to be after cGeophysicsNcFile definition
inline size_t GVar::line_index_start(const size_t& index) const {
	size_t start = Parent.get_line_index_start(index);
	return start;
}

// Defined here only because it needs to be after cGeophysicsNcFile definition
inline size_t GVar::line_index_count(const size_t& index) const {
	size_t count = Parent.get_line_index_count(index);
	return count;
}

};//endname space

