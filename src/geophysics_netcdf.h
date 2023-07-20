/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2016.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#ifndef _geophysics_netcdf_H
#define _geophysics_netcdf_H

#include "stacktrace.h"
#include "netcdf.h"

#include <cassert>
#include <stdexcept>
#include <map>
#include <algorithm>
#include <iomanip> 
#include <memory> 
#include <cfloat>
#include <netcdf>

using namespace netCDF;
using namespace netCDF::exceptions;


#include "general_utils.h"
#include "vector_utils.h"
#include "file_formats.h"
#include "stopwatch.h"
#include "logger.h"
#ifdef HAVE_GDAL
	#include "crs.h"
#endif
#ifdef HAVE_CGAL
	#include "cgal_utils.h"
#endif



#pragma warning (push)
#pragma warning (disable:858) //warning #858: type qualifier on return type is meaningless
#include "marray.hxx"
#pragma warning (pop)
using namespace andres;

extern cLogger glog;//The global instance of the log file manager

constexpr auto DN_POINT = "point";
constexpr auto DN_LINE = "line";

constexpr auto VN_LI_START   = "line_index_start";
constexpr auto VN_LI_COUNT   = "line_index_count";
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

inline NcType nctype(const uint8_t) { return ncUbyte; }
inline NcType nctype(const int8_t) { return ncByte; }
inline NcType nctype(const short) { return ncShort; }
inline NcType nctype(const int) { return ncInt; }
inline NcType nctype(const unsigned int) { return ncUint; }
inline NcType nctype(const float) { return ncFloat; }
inline NcType nctype(const double) { return ncDouble; }
inline NcType nctype(const std::string) { return ncString; }

inline uint8_t defaultmissingvalue(const NcUbyte&) { return (uint8_t) NC_FILL_UBYTE; }
inline int8_t defaultmissingvalue(const NcByte&) { return (int8_t)NC_FILL_BYTE; }
inline short defaultmissingvalue(const NcShort&) { return (short) NC_FILL_SHORT; }
inline int defaultmissingvalue(const NcInt&) { return (int) NC_FILL_INT; }
inline unsigned int defaultmissingvalue(const NcUint&) { return (unsigned int)NC_FILL_UINT; }
inline float defaultmissingvalue(const NcFloat&) { return (float) NC_FILL_FLOAT; }
inline double defaultmissingvalue(const NcDouble&) { return (double) NC_FILL_DOUBLE; }
inline std::string defaultmissingvalue(const NcString&) { return std::string(NC_FILL_STRING); }

class cExportFormat{

public:
	char   form='\0';
	int    width=0;
	int    decimals=0;
	double nullvalue=0;

	cExportFormat()
	{
		_GSTITEM_
	};

	cExportFormat(char _form, int _width, int _decimals, double _nullvalue){
		_GSTITEM_
		form      = _form;
		width     = _width;
		decimals  = _decimals;
		nullvalue = _nullvalue;
	}
	
};

class cGeophysicsNcFile;

class cGeophysicsVar : public NcVar{

private:	
	cGeophysicsNcFile* FilePtr = NULL;
	
public:	

	cGeophysicsVar() : NcVar() {};

	cGeophysicsVar(cGeophysicsNcFile* parentptr, const NcVar& var) : NcVar(var) {
		_GSTITEM_
		FilePtr = parentptr;		
	};

	cGeophysicsVar(const cGeophysicsVar& var) : NcVar(var) {
		_GSTITEM_
		FilePtr = var.FilePtr;
	};

	cGeophysicsVar operator=(const cGeophysicsVar& rhs) {
		_GSTITEM_
		NcVar::operator=(rhs);
		FilePtr = rhs.FilePtr;
		return *this;
	};

	cGeophysicsVar& operator=(cGeophysicsVar& rhs){
		_GSTITEM_
		NcVar::operator=(rhs);
		FilePtr = rhs.FilePtr;
		return *this;
	};

	size_t line_index_start(const size_t& index) const;
	
	size_t line_index_count(const size_t& index) const;

	size_t length(){
		_GSTITEM_
		std::vector<NcDim> dims = getDims();
		size_t len = dims[0].getSize();
		for (size_t di = 1; di < dims.size(); di++){
			len *= dims[di].getSize();
		}
		return len;
	}

	size_t sizeBytes(){
		_GSTITEM_
		return length() * getType().getSize();
	}

	size_t elementspersample() const {		
		_GSTITEM_
		std::vector<NcDim> dims = getDims();
		size_t len = 1;
		if (dims.size() == 1){
			return len;
		}
		else{
			len = dims[1].getSize();
			for (size_t di = 2; di < dims.size(); di++){
				len *= dims[di].getSize();
			}
		}
		return len;
	}

	size_t lineelements(const size_t& lineindex) const {
		_GSTITEM_
		return elementspersample() * line_index_count(lineindex);
	}

	size_t nbands() const {
		_GSTITEM_
		std::vector<NcDim> dims = getDims();
		if (dims.size() == 1) return 1;
		else if (dims.size() == 2) return dims[1].getSize();
		return 0;
		//this need fixing for case of extra dims
	}

	bool isLineVar() const {
		_GSTITEM_
		if (getDimCount() == 0) return false;
		if (getDim(0).getName() == DN_LINE) return true;		
		return false;
	}

	bool isSampleVar() const {
		_GSTITEM_
		if (getDimCount() == 0) return false;
		if (getDim(0).getName() == DN_POINT) return true;		
		return false;
	}

	NcVarAtt add_attribute(const std::string& att, std::string value){
		_GSTITEM_
		return putAtt(att, value);
	}

	NcVarAtt add_standard_name(const std::string& value) {
		_GSTITEM_
			return putAtt(AN_STANDARD_NAME, value);
	}

	NcVarAtt add_long_name(const std::string& value){
		_GSTITEM_
		return putAtt(AN_LONG_NAME, value);
	}

	NcVarAtt add_original_dataset_fieldname(const std::string& value){
		_GSTITEM_
		return putAtt(AN_ORIGINAL_DATASET_FIELDNAME, value);
	}

	NcVarAtt add_units(const std::string& value){
		_GSTITEM_
		return putAtt(AN_UNITS, value);
	}

	NcVarAtt add_description(const std::string& value){
		_GSTITEM_
		return putAtt(AN_DESCRIPTION, value);
	}

	static bool hasAtt(const NcVar& var, const std::string& name) {		
		_GSTITEM_
		nc_type vr_type;
		size_t  vr_len;
		int status = nc_inq_att(var.getParentGroup().getId(), var.getId(), name.c_str(), &vr_type, &vr_len);
		if (status == NC_NOERR) return true;
		else return false;				
	};

	bool hasAtt(const std::string& name) const {		
		_GSTITEM_
		return cGeophysicsVar::hasAtt(*this, name);		
	}

	std::string getStringAtt(const std::string& attname) const {
		_GSTITEM_
		std::string attvalue;
		if (hasAtt(attname)){
			NcVarAtt a = getAtt(attname);
			a.getValues(attvalue);
		}
		return attvalue;
	}

	std::string getUnits() const {
		_GSTITEM_
		return getStringAtt(AN_UNITS);
	}

	std::string getDescription() const {
		_GSTITEM_
		return getStringAtt(AN_DESCRIPTION);
	}

	//static float preferred_float_missing_value() {
	//	_GSTITEM_
	//	return NC_FILL_FLOAT;
	//}

	//static double preferred_double_missing_value() {
	//	_GSTITEM_
	//	return NC_FILL_DOUBLE;
	//}

	double lowest_possible_value() const {
		_GSTITEM_
		const NcType type = getType();
		if (type == ncShort) return SHRT_MIN;
		else if (type == ncUint) return 0;
		else if (type == ncInt) {
			double v = (double)INT_MIN;
			return v;
		}
		else if (type == ncFloat) return -FLT_MAX;
		else return -DBL_MAX;
	}

	double highest_possible_value(){
		_GSTITEM_
		const NcType type = getType();
		if (type == ncShort) return SHRT_MAX;
		else if (type == ncUint) return UINT_MAX;
		else if (type == ncInt) return INT_MAX;
		else if (type == ncFloat) return FLT_MAX;
		else return DBL_MAX;
	}

	void set_default_missingvalue(){		
		_GSTITEM_
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
					throw(std::runtime_error(msg));
				}
		}
		return;
	}

	template<typename T>
	bool add_missing_value(const T& value){
		_GSTITEM_
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
		_GSTITEM_
		if (hasAtt(AN_MISSINGVALUE)){
			T v;
			getAtt(AN_MISSINGVALUE).getValues(&v);
			return v;
		}		
		else{
			const NcType type = getType();
			if (type == ncShort) return (T) NC_FILL_SHORT;
			else if (type == ncInt)  return (T) NC_FILL_INT;
			else if (type == ncFloat) return (T) NC_FILL_FLOAT;
			else if (type == ncDouble) return (T) NC_FILL_DOUBLE;
			else return (T) NC_FILL_SHORT;
		}
	}

	template<typename T>
	bool getAll(std::vector<T>& vals){
		_GSTITEM_
		vals.resize(length());
		getVar(vals.data());
		return true;
	}	

	template<typename T>
	bool minmax(T& minval, T& maxval){
		_GSTITEM_
		std::vector<T> vals;
		getAll(vals);
		T nullv;
		nullv = missingvalue(nullv);
		minval = highest_possible_value();
		maxval = lowest_possible_value();
		for (size_t i = 0; i < vals.size(); i++){
			if (vals[i] == nullv) continue;
			if (vals[i] < minval) minval = vals[i];
			if (vals[i] > maxval) maxval = vals[i];
		}
		return true;
	}

	template<typename T>
	void getLine(const size_t& lineindex, Marray<T>& A) const {
		_GSTITEM_
		if (isNull()){
			std::string msg = _SRC_ + strprint("\nAttempt to read from a Null variable\n");
			throw(std::runtime_error(msg));
		}
		std::vector<NcDim>  dims = getDims();
		std::vector<size_t> start((size_t)getDimCount());
		std::vector<size_t> count((size_t)getDimCount());
		
		
		if (isLineVar()){
			start[0] = lineindex;
			count[0] = 1;
		}
		else{
			start[0] = line_index_start(lineindex);
			count[0] = line_index_count(lineindex);			
		}

		for (size_t i = 1; i < dims.size(); i++){
			start[i] = 0;
			count[i] = dims[i].getSize();
		}

		A.resize(count.data(), count.data() + count.size());
		getVar(start, count, &(A(0)));
	}

	template<typename T>
	bool getRecord(const size_t& record, std::vector<T>& v) const
	{	
		_GSTITEM_
		if (isNull()) { 
			std::string msg = _SRC_ + strprint("\nAttempt to read from a Null variable\n");			
			throw(std::runtime_error(msg));
			return false; 
		}
		std::vector<size_t> startp = { record, 0 };
		std::vector<size_t> countp = { 1,      nbands()};
		v.resize(nbands());
		getVar(startp, countp, v.data());
		return true;
	};

	template<typename T>
	bool putRecord(const size_t& record, const T& v) const
	{
		_GSTITEM_
		if (isNull()) {
			std::string msg = _SRC_ + strprint("\nAttempt to write to a Null variable\n");
			throw(std::runtime_error(msg));
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
		_GSTITEM_
		if (isNull()) {
			std::string msg = _SRC_ + strprint("\nAttempt to write to a Null variable\n");
			throw(std::runtime_error(msg));
			return false;
		}
		std::vector<size_t> startp = { record, 0 };
		std::vector<size_t> countp = { 1,      nbands() };
		assert(countp[1] == v.size());
		putVar(startp, countp, v.data());
		return true;
	};
		

	bool donotexport(){
		_GSTITEM_
		std::vector<std::string> s = {
			DN_POINT, VN_LI_START, VN_LI_COUNT,
			"longitude_first", "longitude_last",
			"latitude_first", "latitude_last",
			"easting_first", "easting_last",
			"northing_first", "northing_last"};

		auto it = std::find(s.begin(), s.end(), getName());
		if (it == s.end()) return false;
		else return true;
	};

	cExportFormat defaultexportformat() const {
		_GSTITEM_
		cExportFormat e;
		nc_type t = getType().getId();
		switch (t){
		case NC_SHORT: e = cExportFormat('I',8,0, -999); break;
		case NC_INT: e = cExportFormat('I',12,0, -999); break;
		case NC_UINT: e = cExportFormat('I', 12, 0, 999); break;
		case NC_FLOAT: e = cExportFormat('F', 10, 4, -999); break;
		case NC_DOUBLE: e = cExportFormat('F', 16, 6, -999); break;
		default: e = cExportFormat('F', 16, 6, -999);  break;
		}
		return e;
	}
};

class cLineVar : public cGeophysicsVar{

private:

public:

	cLineVar(cGeophysicsNcFile* parentptr, const NcVar& var)
		: cGeophysicsVar(parentptr, var)
	{
		_GSTITEM_
	}

	template<typename T>
	bool putAll(std::vector<T> vals){
		_GSTITEM_
		if (isNull()){
			std::string msg = _SRC_ + strprint("\nAttempt to write to a Null variable\n");			
			throw(std::runtime_error(msg));
		}

		if (vals.size() != length()){
			std::string msg = _SRC_ + strprint("\nAttempt to write variable (%s) with non-matching size\n", getName().c_str());			
			throw(std::runtime_error(msg));
		}

		putVar(vals.data());
		return true;
	}

	template<typename T>
	bool getLine(const size_t& lineindex, const size_t& bandindex, T& val){
		_GSTITEM_
		if (isNull()){ return false; }
		std::vector<size_t> startp = { lineindex, bandindex };
		std::vector<size_t> countp = { 1, 1 };
		getVar(startp, countp, &val);
		return true;
	}

	template<typename T>
	bool getLine(const size_t& lineindex, T& val) const {
		_GSTITEM_
		return getLine(lineindex, 0, val);
	}

	template<typename T>
	bool getLine(const size_t& lineindex, const size_t& bandindex, std::vector<T>& vals) const {
		_GSTITEM_
		if (isNull()){ return false; }
		std::vector<size_t> startp = { lineindex, bandindex };
		std::vector<size_t> countp = { 1, 1 };
		vals.resize(countp[0] * countp[1]);
		getVar(startp, countp, vals.data());
		return true;
	}

	template<typename T>
	bool getLine(const size_t& lineindex, std::vector<T>& vals) const {
		_GSTITEM_
		return getLine(lineindex, 0, vals);
	}

	template<typename T>
	T getSample(const size_t& lineindex, const size_t& sampleindex, const size_t& bandindex, T& val) const {
		_GSTITEM_
		if (isNull()){ return false; }
		std::vector<size_t> startp = { lineindex, bandindex };
		std::vector<size_t> countp = { 1, 1 };		
		getVar(startp, countp, &val);
		return val;
	};

};

class cSampleVar : public cGeophysicsVar{

private:

public:	

	cSampleVar() : cGeophysicsVar()
	{
		_GSTITEM_
	}

	cSampleVar(cGeophysicsNcFile* parentptr, const NcVar& var) 
		: cGeophysicsVar(parentptr, var)
	{
		_GSTITEM_
	}
	

	template<typename T>
	bool putAll(const std::vector<T>& vals){
		_GSTITEM_
		if (isNull()){
			std::string msg = _SRC_ + strprint("\nAttempt to write to a Null variable\n");			
			throw(std::runtime_error(msg));
		}

		if (vals.size() != length()){
			std::string msg = _SRC_ + strprint("\nAttempt to write variable (%s) with non-matching size\n", getName().c_str());			
			throw(std::runtime_error(msg));
		}

		putVar(vals.data());
		return true;
	}
	
	template<typename T>
	bool putLineBand(const size_t& lineindex, const size_t& bandindex, const std::vector<T>& vals){
		_GSTITEM_
		if (isNull()){
			std::string msg = _SRC_ + strprint("\nAttempt to write to a Null variable\n");			
			throw(std::runtime_error(msg));
		}

		std::vector<NcDim>  dims = getDims();
		if (dims.size() > 2){
			std::string msg = _SRC_ + strprint("\nAttempt to use putLineBand() to write to a variable with more than 2 dimensions\n");			
			throw(std::runtime_error(msg));
		}

		if (vals.size() != line_index_count(lineindex)){
			std::string msg = _SRC_ + strprint("\nAttempt to write line/band of variable (%s) with non-matching size\n", getName().c_str());
			throw(std::runtime_error(msg));
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
	bool putLine(const size_t& lineindex, const std::vector<T>& vals){		
		_GSTITEM_
		if (isNull()){
			std::string msg = _SRC_ + strprint("\nAttempt to write to a Null variable\n");			
			throw(std::runtime_error(msg));
		}

		size_t sz = lineelements(lineindex);
		if (vals.size() != sz){
			std::string msg = _SRC_ + strprint("\nAttempt to write line/band of variable (%s) with non-matching size\n", getName().c_str());			
			throw(std::runtime_error(msg));
		}
		
		std::vector<NcDim>  dims = getDims();
		std::vector<size_t> start(getDimCount());
		std::vector<size_t> count(getDimCount());
		start[0] = line_index_start(lineindex);
		count[0] = line_index_count(lineindex);
		for (size_t i = 1; i < dims.size(); i++){
			start[i] = 0;
			count[i] = dims[i].getSize();
		}		
		putVar(start, count, vals.data());
		return true;		
	}

	template<typename T>
	bool getLine(const size_t& lineindex, std::vector<T>& vals){				
		_GSTITEM_
		if (isNull()){
			std::string msg = _SRC_ + strprint("\nAttempt to read from a Null variable\n");			
			throw(std::runtime_error(msg));
		}
		
		std::vector<NcDim>  dims = getDims();
		std::vector<size_t> start((size_t)getDimCount());
		std::vector<size_t> count((size_t)getDimCount());
		start[0] = line_index_start(lineindex);
		count[0] = line_index_count(lineindex);		
		for (size_t i = 1; i < dims.size(); i++){
			start[i] = 0;
			count[i] = dims[i].getSize();
		}		
		size_t sz = lineelements(lineindex);		
		vals.resize(sz);		
		getVar(start, count, vals.data());		
		return true;		
	}

	template<typename T>	
	void getLine_temp(const size_t& lineindex, Marray<T>& A) const {
		_GSTITEM_
		if (isNull()){
			std::string msg = _SRC_ + strprint("\nAttempt to read from a Null variable\n");			
			throw(std::runtime_error(msg));
		}		
		std::vector<NcDim>  dims = getDims();				
		std::vector<size_t> start((size_t)getDimCount());
		std::vector<size_t> count((size_t)getDimCount());		
		start[0] = line_index_start(lineindex);
		count[0] = line_index_count(lineindex);
		for (size_t i = 1; i < dims.size(); i++){
			start[i] = 0;
			count[i] = dims[i].getSize();
		}		
		A.resize(count.data(), count.data() + count.size());
		size_t sz = lineelements(lineindex);				
		getVar(start, count, &(A(0)));		
	}

	template<typename T>
	bool getSample(const size_t& lineindex, const size_t& sampleindex, const size_t& bandindex, T& val) const {		
		_GSTITEM_		
		if (isNull()){ return false; }
		std::vector<size_t> startp = { line_index_start(lineindex)+sampleindex, bandindex };
		std::vector<size_t> countp = { 1, 1 };	
		getVar(startp, countp, &val);
		return true;
	};

	template<typename T>
	bool getPoint(const size_t& pointindex, std::vector<T>& v) const
	{		
		_GSTITEM_
		if (isNull()) return false;
		getRecord(pointindex, v);
		return true;
	};

	template<typename T>
	bool putPoint(const size_t& pointindex, const std::vector<T>& v) const
	{
		_GSTITEM_
		if (isNull()) return false;
		putRecord(pointindex, v);
		return true;
	};

};

class cGeophysicsNcFile : public NcFile {

private:

	std::vector<unsigned int> line_index_start;
	std::vector<unsigned int> line_index_count;
	std::vector<unsigned int> line_number;

	NcDim dim_sample() { return getDim(DN_POINT); }

	NcDim dim_line() { return getDim(DN_LINE); }

	bool InitialiseExisting() {
		_GSTITEM_
		if (readLineIndex() == false) return false;
		if (getLineNumbers(line_number) == false) return false;
		return true;
	}

	bool readLineIndex() {
		_GSTITEM_
		if (hasVar(VN_LI_COUNT)) {
			cLineVar vc = getLineVar(VN_LI_COUNT);
			if (vc.getAll(line_index_count) == false)return false;

			cLineVar vs = getLineVar(VN_LI_START);
			if (vs.getAll(line_index_start) == false)return false;
		}
		else if (hasVar(VN_LINE_INDEX)) {
			cLineVar vl = getLineVar(VN_LINE_INDEX);
			std::vector<unsigned int> line_index;
			if (vl.getAll(line_index) == false)return false;
			set_start_count(line_index);
		}

		cLineVar v = getLineVar(DN_LINE);
		v.getAll(line_number);
		return true;
	}

	static std::vector<unsigned int> compute_line_index(std::vector<unsigned int>& count)
	{
		_GSTITEM_
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
		_GSTITEM_
		line_index_start = compute_start_from_line_index(line_index);
		line_index_count = compute_count_from_start(line_index_start, line_index.size());
		return;
	}

	static std::vector<unsigned int> compute_start_from_line_index(const std::vector<unsigned int>& line_index)
	{
		_GSTITEM_
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
		_GSTITEM_
		const size_t nl = count.size();
		std::vector<unsigned int> start(nl);
		start[0] = 0;
		for (size_t i = 1; i < nl; i++) {
			start[i] = start[i - 1] + count[i-1];
		}
		return start;
	}

	static std::vector<unsigned int> compute_count_from_start(const std::vector<unsigned int>& start, size_t nsamplestotal)
	{
		_GSTITEM_
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
		_GSTITEM_
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
		_GSTITEM_
		size_t n = nelements(v);
		n *= v.getType().getSize();
		return n;
	}

	static size_t nbytes(const NcAtt& a)
	{
		_GSTITEM_
		size_t n = a.getAttLength() * a.getType().getSize();
		return n;
	}



public:

	//Do not allow implicit definition of copy constructor or assignment operators	
	cGeophysicsNcFile& operator=(const NcGroup & rhs) = delete;
	cGeophysicsNcFile& operator=(const NcFile & rhs) = delete;
	cGeophysicsNcFile& operator=(const cGeophysicsNcFile & rhs) = delete;

	cGeophysicsNcFile(const NcGroup& rhs) = delete;
	cGeophysicsNcFile(const NcFile& rhs) = delete;
	cGeophysicsNcFile(const cGeophysicsNcFile& rhs) = delete;
	
	// Constructor generates a null object.
	cGeophysicsNcFile() : NcFile() {} // invoke base class constructor	

	//Open existing file constructor
	cGeophysicsNcFile(const std::string& ncpath, const netCDF::NcFile::FileMode& filemode = netCDF::NcFile::FileMode::read)
		: netCDF::NcFile(ncpath, filemode)
	{
		_GSTITEM_
		open(ncpath,filemode);		
	};

	//Destructor
	~cGeophysicsNcFile(){};	

	size_t get_line_index_start(const size_t& li) const { return line_index_start[li]; }
	size_t get_line_index_count(const size_t& li) const { return line_index_count[li]; }

	std::string pathname() const {
		size_t len;
		nc_inq_path(getId(), &len, nullptr);
		char* p = new char[len + 1];
		nc_inq_path(getId(), nullptr, p);
		std::string path(p);
		delete p;
		return path;
	}

	bool isopen() const {
		_GSTITEM_
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
		_GSTITEM_
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

	bool InitialiseNew(const std::vector<size_t>& linenumbers, const std::vector<size_t>& linesamplecount){
		_GSTITEM_
		size_t n = linenumbers.size();
		std::vector<unsigned int> uint_linenumbers(n);
		std::vector<unsigned int> uint_linesamplecount(n);
		for (size_t i = 0; i < n; i++){
			uint_linenumbers[i] = (unsigned int)linenumbers[i];
			uint_linesamplecount[i] = (unsigned int)linesamplecount[i];
		}
		return InitialiseNew(uint_linenumbers, uint_linesamplecount);
	}

	bool InitialiseNew(const std::vector<unsigned int>& linenumbers, const std::vector<unsigned int>& linesamplecount){
		_GSTITEM_
		const size_t nl = linenumbers.size();
		line_number = linenumbers;
		line_index_count = linesamplecount;
		line_index_start = compute_start_from_count(line_index_count);
		
		size_t nsamples = sum(line_index_count);		
		NcDim ds = addDim(DN_POINT, nsamples);
		NcDim dl = addDim(DN_LINE, nl);

		std::vector<unsigned int> line_index = compute_line_index(line_index_count);
		add_line_index(line_index);
		//add_line_index_start(line_index);
		//add_line_index_count(line_index);
		//add_point_variable();		
		add_line_number_variable(linenumbers);
		return true;
	}

	void add_line_index(const std::vector<unsigned int> line_index)
	{
		_GSTITEM_
		cSampleVar v = addSampleVar(VN_LINE_INDEX, ncUint);
		v.putVar(line_index.data());
		v.add_long_name(LN_LINE_INDEX);		
		v.add_description("zero-based index of line associated with point");
		//v.add_units("1");		
	}

	//Deprecated
	void add_line_index_start(const std::vector<unsigned int> line_index_start)
	{
		_GSTITEM_
		cLineVar vstart = addLineVar(VN_LI_START, ncUint);		
		vstart.putVar(line_index_start.data());
		vstart.add_long_name(VN_LI_START);		
		vstart.add_description("zero-based index of the first sample in the line");
		//vstart.add_units("1");
	}

	//Deprecated
	void add_line_index_count(const std::vector<unsigned int> line_index_count)
	{
		_GSTITEM_
		cLineVar vcount = addLineVar(VN_LI_COUNT, ncUint);
		vcount.putVar(line_index_count.data());
		vcount.add_long_name(VN_LI_COUNT);
		vcount.add_description("number of samples in the line");
		//vcount.add_units("1");
	}

	//Deprecated
	void add_point_variable()
	{
		_GSTITEM_
		std::vector<unsigned int> sample = increment((unsigned int)ntotalsamples(), (unsigned int)0, (unsigned int)1);
		cSampleVar v = addSampleVar(DN_POINT, ncUint);		
		v.putVar(sample.data());
		v.add_long_name(LN_SAMPLE_NUMBER);
		v.add_description("sequential point number");
		//v.add_units("1");
	}

	void add_line_number_variable(const std::vector<unsigned int> linenumbers)
	{
		_GSTITEM_
		cLineVar vline = addLineVar(DN_LINE, ncUint);
		vline.putVar(linenumbers.data());
		vline.add_long_name(LN_LINE_NUMBER);
		vline.add_description("flight line number");
		//vline.add_units("1");		
	}

	template<typename T>
	bool isinlist(const std::vector<T>& vec, const T& val){		
		if (std::find(vec.begin(), vec.end(), val) == vec.end()) return false;
		return true;
	}
	
	bool subsample(const cGeophysicsNcFile& srcfile, const size_t subsamplerate, std::vector<std::string> include_varnames, std::vector<std::string> exclude_varnames)
	{
		_GSTITEM_
		unsigned int nsnew = (unsigned int) std::ceil((double)srcfile.ntotalsamples() / (double)subsamplerate);
		
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
		for (size_t li = 0; li < start.size(); li++){						
			linenumber[li] = srcfile.line_number[ss_srclineindex[start[li]]];
		}
						
		InitialiseNew(linenumber,count);
						
		auto dm = srcfile.getDims();
		for (auto dit = dm.begin(); dit != dm.end(); dit++) {
			NcDim& srcdim = dit->second;
			const std::string dimname = srcdim.getName();
			if (hasDim(dimname)) continue;
			else {
				addDim(dimname,srcdim.getSize());
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
			if (include_varnames.size()>0 && incstatus==false) {
				status = false;
			}

			if (exclude_varnames.size() > 0 && excstatus==true) {
				status = false;
			}
						
			if (status) {
				copy_var(subsamplerate, srcvar);
			}
		}
		return true;
	}

#ifdef HAVE_GDAL
	//Convert legacy file
	bool convert_legacy(const cGeophysicsNcFile& srcfile)
	{
		_GSTITEM_
		//bookmark
		InitialiseNew(srcfile.line_number, srcfile.line_index_count);
		copy_global_atts(srcfile);
		copy_dims(srcfile);
		auto vm = srcfile.getVars();
		for (auto vit = vm.begin(); vit != vm.end(); vit++){
			NcVar& srcvar = vit->second;			
			if (srcvar.getName() == VN_LI_START) continue;
			if (srcvar.getName() == VN_LI_COUNT) continue;
			if (srcvar.getName() == DN_POINT) continue;

			if (srcvar.getName() == "crs"){
				NcVarAtt a = srcvar.getAtt("epsg_code");
				std::string epsg_string;
				a.getValues(epsg_string);
				int epsgcode = cCRS::epsgcode(epsg_string);
				addCRS(epsgcode);
			}
			else{
				copy_var(1,srcvar);
				NcVar v = getVar(srcvar.getName());
				nc_rename_att(getId(), v.getId(), "standard_name", "long_name");
				
				if (v.getName() == "latitude"){
					v.putAtt("standard_name", "latitude");
				}

				if (v.getName() == "longitude"){
					v.putAtt("standard_name", "longitude");
				}
							
				//Remove units = "1"				
				if (cGeophysicsVar::hasAtt(v,AN_UNITS)){					
					NcVarAtt a = v.getAtt(AN_UNITS);
					if (!a.isNull()){
						std::string units;
						a.getValues(units);
						if (units == "1"){
							nc_del_att(getId(), v.getId(), AN_UNITS);
						}
					}
				}
			}					
		}		
		return true;
	}
#endif

	bool copy_global_atts(const cGeophysicsNcFile& srcfile)
	{
		_GSTITEM_
		auto m = srcfile.getAtts();
		for (auto it = m.begin(); it != m.end(); it++){
			NcGroupAtt& srcatt = it->second;			
			size_t len = nbytes(srcatt);
			if (len > 0){
				const std::vector<uint8_t> buf(len);
				srcatt.getValues((void*)buf.data());
				putAtt(srcatt.getName(), srcatt.getType(), srcatt.getAttLength(), (void*)buf.data());				
			}
			
		}
		return true;
	}

	bool copy_dims(const cGeophysicsNcFile& srcfile)
	{
		_GSTITEM_
		auto dm = srcfile.getDims();
		for (auto dit = dm.begin(); dit != dm.end(); dit++){
			NcDim& srcdim = dit->second;
			if (hasDim(srcdim.getName())) continue;
			addDim(srcdim.getName(), srcdim.getSize());
		}
		return true;
	}


	bool copy_var(const size_t& subsample, const NcVar& srcvar, std::string newname = std::string())
	{
		_GSTITEM_
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
				count[i]  = (size_t) std::ceil(srcdims[i].getSize() / (double)subsample);
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
		if (bufsize > 0){
			std::vector<uint8_t> buf(bufsize);
			srcvar.getVar(start, count, stride,(void*)buf.data());			
			v.putVar(start, count, stride_out,(void*)buf.data());
			//v.putVar((void*)buf.data());
		}
		return true;
	}

	bool copy_varatts(const NcVar& srcvar, const NcVar& dstvar)
	{
		_GSTITEM_
		auto am = srcvar.getAtts();
		for (auto ait = am.begin(); ait != am.end(); ait++){
			NcVarAtt& srcatt = ait->second;
			size_t len = nbytes(srcatt);
			if (len > 0){
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
		_GSTITEM_
		const size_t n = line_index_start.size() - 1;
		for (size_t k = 0; k < n; k++) {
			if ((int)line_index_start[k+1] > pointindex) return k;
		}
		if ((int)(line_index_start[n] + line_index_count[n]) > pointindex) return n;
		return undefinedvalue<size_t>();
	};

	size_t getLineIndex(const int& linenumber){
		_GSTITEM_
		auto it = std::find(line_number.begin(), line_number.end(), linenumber);
		return (size_t)(it - line_number.begin());
	}

	bool hasVar(const std::string& varname){
		_GSTITEM_
		NcVar v = getVar(varname);
		if (v.isNull()) return false;
		return true;
	}

	bool hasVarCaseInsensitive(std::string& varname) {
		_GSTITEM_
		NcVar v = getVarByNameCaseInsensitive(varname);
		if (v.isNull()) return false;
		return true;
	}

	NcVar getVarByNameCaseInsensitive(std::string& varname) {
		_GSTITEM_
		std::multimap<std::string, NcVar> vars = getVars();
		for (auto vit = vars.begin(); vit != vars.end(); vit++) {
			if (strcasecmp(varname, vit->first) == 0) {
				varname = vit->first;
				return vit->second;
			}
		}
		return NcVar();
	}

	bool hasDim(const std::string& dimname){
		_GSTITEM_
		NcDim d = getDim(dimname);
		if (d.isNull()) return false;
		return true;
	}

	NcDim addDim(const std::string& dimname, const size_t& dimsize) {
		_GSTITEM_
			if (hasDim(dimname)) {
				NcDim d = getDim(dimname);
				if (d.getSize() == dimsize) return d;
				else {
					std::string msg = _SRC_ + strprint("\nAttempt to add dimension (%s) with size (%zu) unequal to existing dimension (%zu) with same name\n", dimname.c_str(), dimsize, d.getSize());
					throw(std::runtime_error(msg));
				}
			}
			else {
				NcDim d = NcGroup::addDim(dimname, dimsize);
				return d;
			}
	}
	
	NcVar getVarByAtt(const std::string& att_name, const std::string& att_value){
		_GSTITEM_
		std::multimap<std::string, NcVar> vars = getVars();
		for (auto vit = vars.begin(); vit != vars.end(); vit++){
			std::map<std::string, NcVarAtt> atts = vit->second.getAtts();
			for (auto ait = atts.begin(); ait != atts.end(); ait++){
				std::string name = ait->second.getName();
				if (name == att_name){
					std::string value;
					ait->second.getValues(value);
					if (value == att_value){
						return vit->second;
					}
				}
			}
		}
		return NcVar();
	}

	NcVar getVarByLongName(const std::string& att_value){
		_GSTITEM_
		return getVarByAtt(AN_LONG_NAME, att_value);
	}

	std::string getVarNameByAtt(const std::string& att_name, const std::string& att_value){
		_GSTITEM_
		NcVar var = getVarByAtt(att_name, att_value);
		if (var.isNull()){
			return std::string();
		}
		else return var.getName();
	}

	std::string getVarNameByLongName(const std::string& att_value){
		_GSTITEM_
		return getVarNameByAtt(AN_LONG_NAME, att_value);
	}

	std::string getVarNameByStandardName(const std::string& att_value) {
		_GSTITEM_
		return getVarNameByAtt(AN_STANDARD_NAME, att_value);
	}

	template<typename T>
	bool getLineNumbers(std::vector<T>& vals){
		_GSTITEM_
		NcVar v = getVarByLongName(LN_LINE_NUMBER);
		if (v.isNull() == false){
			cLineVar var(this, v);
			return var.getAll(vals);
		}
		return false;
	}
	
	std::vector<int> getLineNumbers(){
		_GSTITEM_
		std::vector<int> num;
		getLineNumbers(num);
		return num;
	}

	template<typename T>
	bool getFlightNumbers(std::vector<T>& vals){
		_GSTITEM_
		NcVar v = getVarByLongName(LN_FLIGHT_NUMBER);
		if (v.isNull() == false){
			cLineVar var(this, v);
			return var.getAll(vals);
		}
		return false;
	}

	std::vector<int> getFlightNumbers(){
		_GSTITEM_
		std::vector<int> num;
		getFlightNumbers(num);
		return num;
	}

	template<typename T>
	bool getDataByLineIndex(const cSampleVar& var, const size_t& lineindex, std::vector<T>& vals){		
		_GSTITEM_
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
	bool getDataByLineNumber(const std::string& varname, const size_t& linenumber, std::vector<T>& vals){
		size_t index = getLineIndex(linenumber);
		if (index >= nlines())return false;
		return getDataByLineIndex(varname, index, vals);
	}

	template<typename T>
	bool getDataByLineIndex(const std::string& varname, const size_t& lineindex, std::vector<T>& vals){
		cSampleVar var = getSampleVar(varname);
		return getDataByLineIndex(var, lineindex, vals);		
	}

	template<typename T>
	bool getDataByLineIndex(const std::string& varname, const size_t& lineindex, std::vector<std::vector<T>>& vals){
		NcVar var = NcFile::getVar(varname);
		std::vector<NcDim> dims = var.getDims();
		size_t nd = dims.size();
		if (nd != 2) return false;

		size_t nsamples = line_index_count[lineindex];
		size_t nbands   = dims[1].getSize();

		std::vector<size_t> start(2);
		std::vector<size_t> count(2);
		start[0] = line_index_start[lineindex];
		count[0] = nsamples;

		vals.resize(nbands);
		for (size_t bi = 0; bi < nbands; bi++){
			start[1] = bi;
			count[1] = 1;
			vals[bi].resize(nsamples);
			var.getVar(start, count, vals[bi].data());
		}
		return true;
	}

	template<typename T>
	bool getDataByPointIndex(const std::string& varname, const size_t& pointindex, std::vector<T>& vals) {
		cGeophysicsVar var(this,getVar(varname));
		if (var.isNull()) {			
			std::string msg = _SRC_ + strprint("\nAttempt to read variable (%s)\n", varname.c_str());			
			throw(std::runtime_error(msg));
		}

		size_t record;
		if (isScalarVar(var)){			
			vals.resize(1);			
			NcVar v = getVar(varname);
			v.getVar(vals.data());
			return true;
		}
		else if(isSampleVar(var)) record = pointindex;
		else if(isLineVar(var))   record = getLineIndexByPointIndex((int)pointindex);
		else {
			std::string msg = _SRC_ + strprint("\nVariable %s is neither a \"point\" or a \"line\" variable\n", varname.c_str());
			throw(std::runtime_error(msg));
		}
		return var.getRecord(record,vals);
	}

	template<typename T>
	NcDim addDimVar(const std::string& dimname, const std::vector<T>& dimvals){

		size_t dimsize = dimvals.size();
		NcDim dim = getDim(dimname);
		if (dim.isNull()){
			dim = addDim(dimname, dimsize);
		}
		else{
			if (dim.getSize() != dimsize){
				std::string msg = _SRC_ + strprint("\nAttempt to add new dimension (%s) with different size to the existing\n", dimname.c_str());				
				throw(std::runtime_error(msg));
			}
		}

		NcVar var = getVar(dimname);
		if (var.isNull()){
			var = addVar(dimname, nctype(dimvals[0]), dim);
		}
		else{
			if (var.getDim(0).getSize() != dimsize){
				std::string msg = _SRC_ + strprint("\nAttempt to add new dimension (%s) with different size to the existing\n", dimname.c_str());				
				throw(std::runtime_error(msg));
			}
		}
		var.putVar(dimvals.data());
		return dim;
	}
		
	cSampleVar getSampleVar(const std::string& name){					
		if (getVarCount() == 0) return cSampleVar(this, NcVar());
		return cSampleVar(this, getVar(name));
	}

	cLineVar getLineVar(const std::string& name){
		if (getVarCount() == 0) return cLineVar(this, NcVar());
		return cLineVar(this, getVar(name));
	}
	
	cSampleVar addSampleVar(const std::string& name, const NcType& type, const std::vector<NcDim>& dims){

		cSampleVar var = getSampleVar(name);
		if (var.isNull() == false){
			return var;
		}

		std::vector<NcDim> vardims = { dim_sample() };
		for (size_t i = 0; i < dims.size(); i++){
			vardims.push_back(dims[i]);
		}
		var = cSampleVar(this, addVar(name, type, vardims));
		var.set_default_missingvalue();
		return var;
	}

	cSampleVar addSampleVar(const std::string& name, const NcType& type, const NcDim& banddim = NcDim()){
		std::vector<NcDim> dims;
		if (banddim.isNull() == false){
			dims.push_back(banddim);
		}
		return addSampleVar(name, type, dims);
	}

	cLineVar addLineVar(const std::string& name, const NcType& type, const std::vector<NcDim>& dims){

		cLineVar var = getLineVar(name);
		if (var.isNull() == false) return var;

		std::vector<NcDim> vardims = { dim_line() };
		for (size_t i = 0; i < dims.size(); i++){
			vardims.push_back(dims[i]);
		}
		var = cLineVar(this, addVar(name, type, vardims));
		var.set_default_missingvalue();			
		return var;		
	}

	cLineVar addLineVar(const std::string& name, const NcType& type, const NcDim& banddim = NcDim()){
		std::vector<NcDim> dims;
		if (banddim.isNull() == false){
			dims.push_back(banddim);
		}
		return addLineVar(name, type, dims);
	}
	
	bool findNonNullLineStartEndPoints(const std::string& xvar, const std::string& yvar,
		std::vector<double>& x1, std::vector<double>& x2, std::vector<double>& y1, std::vector<double>& y2){
		x1.resize(nlines());
		x2.resize(nlines());
		y1.resize(nlines());
		y2.resize(nlines());

		cSampleVar vx = getSampleVar(xvar);
		cSampleVar vy = getSampleVar(yvar);
		double nvx = vx.missingvalue(nvx);
		double nvy = vy.missingvalue(nvy);
		for (size_t li = 0; li < nlines(); li++){			
			std::vector<double> x;
			std::vector<double> y;			
			vx.getLine(li, x);			
			vy.getLine(li, y);
			
			const size_t ns = x.size();

			x1[li] = nvx;
			y1[li] = nvy;			
			for (size_t si = 0; si < ns; si++){
				if (x[si] != nvx && y[si] != nvy){
					x1[li] = x[si];
					y1[li] = y[si];
					break;
				}
			}			

			x2[li] = nvx;
			y2[li] = nvy;			
			for (size_t si = ns - 1; ; si--){
				if (x[si] != nvx && y[si] != nvy){
					x2[li] = x[si];
					y2[li] = y[si];					
					break;
				}
				if (si == 0)break;//avoid endless loop
			}			
		}		
		return true;
	}

#ifdef HAVE_GDAL
	bool addCRS(const int epsgcode){

		#ifdef _WIN32 
			NcVar v = NcFile::addVar("crs", ncByte);
		#else 			
			std::vector<NcDim> d;
			NcVar v = NcFile::addVar("crs", ncByte, d);
		#endif

		v.putAtt(AN_LONG_NAME,"coordinate_reference_system");
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

	bool addLineStartEndPointsLL(){

		if (hasVar("longitude_first")){
			std::string msg = _SRC_ + strprint("\nWarning: Variable longitude_first already exists\n");
			glog.logmsg(msg);
			return false;
		}

		if (hasVar("latitude_first")){
			std::string msg = _SRC_ + strprint("\nWarning: Variable latitude_first already exists\n");
			glog.logmsg(msg);
			return false;
		}

		std::vector<double> x1;
		std::vector<double> x2;
		std::vector<double> y1;
		std::vector<double> y2;

		std::string xvarname = getVarNameByLongName("longitude");
		std::string yvarname = getVarNameByLongName("latitude");		
		findNonNullLineStartEndPoints(xvarname,yvarname, x1, x2, y1, y2);

		cLineVar vx1 = addLineVar("longitude_first", ncDouble);
		vx1.add_missing_value(defaultmissingvalue(ncDouble));
		vx1.add_long_name("longitude_first");
		vx1.add_description("first non-null longitude coordinate in the line");
		vx1.add_units("degree_east");
		vx1.putAll(x1);

		cLineVar vx2 = addLineVar("longitude_last", ncDouble);
		vx2.add_missing_value(defaultmissingvalue(ncDouble));		
		vx2.add_long_name("longitude_last");
		vx2.add_description("last non-null longitude coordinate in the line");
		vx2.add_units("degree_east");
		vx2.putAll(x2);

		cLineVar vy1 = addLineVar("latitude_first", ncDouble);
		vy1.add_missing_value(defaultmissingvalue(ncDouble));
		vy1.add_long_name("latitude_first");
		vy1.add_description("first non-null latitude coordinate in the line");
		vy1.add_units("degree_north");
		vy1.putAll(y1);

		cLineVar vy2 = addLineVar("latitude_last", ncDouble);
		vy2.add_missing_value(defaultmissingvalue(ncDouble));
		vy2.add_long_name("latitude_last");
		vy2.add_description("last non-null latitude coordinate in the line");
		vy2.add_units("degree_north");
		vy2.putAll(y2);
		return true;
	}

	bool addLineStartEndPointsEN(){

		if (hasVar("easting_first")){
			std::string msg = _SRC_ + strprint("\nWarning: Variable easting_first already exists\n");
			glog.logmsg(msg);
			return false;
		}

		if (hasVar("northing_first")){
			std::string msg = _SRC_ + strprint("Warning: Variable northing_first already exists\n");
			glog.logmsg(msg);
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
			glog.logmsg(msg);
			return false;
		}		

		cLineVar vx1 = addLineVar("easting_first", ncDouble);
		vx1.add_missing_value(defaultmissingvalue(ncDouble));
		vx1.add_long_name("easting_first");
		vx1.add_description("first non-null easting coordinate in the line");
		vx1.add_units("m");
		vx1.putAll(x1);

		cLineVar vx2 = addLineVar("easting_last", ncDouble);
		vx2.add_missing_value(defaultmissingvalue(ncDouble));
		vx2.add_long_name("easting_last");
		vx2.add_description("last non-null easting coordinate in the line");
		vx2.add_units("m");
		vx2.putAll(x2);

		cLineVar vy1 = addLineVar("northing_first", ncDouble);
		vy1.add_missing_value(defaultmissingvalue(ncDouble));
		vy1.add_long_name("northing_first");
		vy1.add_description("first non-null northing coordinate in the line");
		vy1.add_units("m");
		vy1.putAll(y1);

		cLineVar vy2 = addLineVar("northing_last", ncDouble);
		vy2.add_missing_value(defaultmissingvalue(ncDouble));
		vy2.add_long_name("northing_last");
		vy2.add_description("last non-null northing coordinate in the line");
		vy2.add_units("m");
		vy2.putAll(y2);
		return true;
	}
	
#ifdef HAVE_CGAL
	bool addAlphaShapePolygon(const std::string xvarname, const std::string yvarname){
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
		for (size_t i = 0; i < nv; i++){
			poly[i * 2] = px[i];
			poly[i * 2 + 1] = py[i];
		}

		std::vector<NcDim> dims;
		dims.push_back(addDim("polygonvertex", nv));
		dims.push_back(addDim("polygonordinate", 2));

		cGeophysicsVar v(this, addVar("bounding_polygon", ncDouble, dims));
		v.add_long_name("bounding_polygon");
		v.add_description("bounding polygon of survey");
		v.add_units("degree");
		v.putVar(poly.data());
		return true;
	}
#endif
	
	bool minmax(const std::string& varname, double& minval, double& maxval){		
		cSampleVar var = getSampleVar(varname);
		var.minmax(minval, maxval);		
		return true;
	}
	
	bool addGeospatialMetadataItem(const std::string& varname, const std::string& label, const std::string& units){
				
		if (hasVar(varname) == false)return false;

		double vmin, vmax;
		cSampleVar var = getSampleVar(varname);
		var.minmax(vmin, vmax);
		std::string s = "geospatial_" + label;
		putAtt(s + "_min", ncDouble, vmin);
		putAtt(s + "_max", ncDouble, vmax);
		putAtt(s + "_units", units);
		putAtt(s + "_resolution", "point");
		if (label == "vertical"){
			putAtt(s + "_positive", "up");
		}		
		return true;		
	}

#ifdef HAVE_GDAL
	bool addGeospatialMetadataXY(){		
		bool status;
		if (hasVar("longitude") && hasVar("latitude")){
			status = addGeospatialMetadataItem("longitude", "lon", "degrees_east");
			if (status == false)return status;
			
			status = addGeospatialMetadataItem("latitude", "lat", "degrees_north");
			if (status == false)return status;			
			addCRS(erm2epsgcode("GDA94"));
		}
		else if (hasVar("easting") && hasVar("northing")){
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

	bool addGeospatialMetadataVertical(){

		if (hasVar("height")){
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
		if (var.getDim(0).getName() == DN_LINE){
			return true;
		}
		return false;
	}

	bool isSampleVar(const NcVar& var) const {
		if (var.getDimCount() == 0) return false;
		if (var.getDim(0).getName() == DN_POINT){
			return true;
		}
		return false;
	}

	std::vector<NcDim> getAllDims() {
		_GSTITEM_
		std::vector<NcDim> dims;
		std::multimap<std::string, NcDim> m = getDims();
		for (auto it = m.begin(); it != m.end(); it++){
			dims.push_back(it->second);
		}
		return dims;
	};

	std::vector<NcVar> getAllVars() {
		_GSTITEM_
		std::vector<NcVar> vars;
		std::multimap<std::string, NcVar> vm = getVars();
		for (auto it = vm.begin(); it != vm.end(); it++){						
			vars.push_back(it->second);
		}		
		return vars;
	};

	std::vector<cLineVar> getLineVars() {
		_GSTITEM_
		std::vector<cLineVar> vars;
		std::multimap<std::string, NcVar> vm = getVars();
		for (auto it = vm.begin(); it != vm.end(); it++){			
			if (isLineVar(it->second)){
				cLineVar v(this,it->second);
				vars.push_back(v);
			}			
		}
		return vars;
	};

	std::vector<cSampleVar> getSampleVars(){
		_GSTITEM_
		std::vector<cSampleVar> vars;
		std::multimap<std::string, NcVar> vm = getVars();
		for (auto it = vm.begin(); it != vm.end(); it++){			
			if (isSampleVar(it->second)){
				cSampleVar v(this, it->second);
				vars.push_back(v);
			}
		}
		return vars;
	};

	bool export_ASEGGDF2(const std::string& datfilepath, const std::string& dfnfilepath){		
		_GSTITEM_
		std::ofstream of(datfilepath);
		of << std::fixed;

		std::vector<cLineVar>   lvars = getLineVars();
		std::vector<cSampleVar> svars = getSampleVars();
		std::vector<cGeophysicsVar> vars;
		for (size_t i = 0; i < lvars.size(); i++){
			if (lvars[i].donotexport() == false){
				vars.push_back(lvars[i]);
			}
		}
		for (size_t i = 0; i < svars.size(); i++){
			std::string vname = svars[i].getName();
			if (svars[i].donotexport() == false){
				vars.push_back(svars[i]);
			}
		}

		const size_t nvars = vars.size(); // number of vars to be exported
		std::vector<double>  mval(nvars); // missing value
		std::vector<bool>    islv(nvars); // is it a line var
		std::vector<cExportFormat> efmt(nvars);//format
		
		cOutputFileInfo I;
		for (size_t vi = 0; vi < nvars; vi++){
			cGeophysicsVar& v = vars[vi];
			
			if (v.isLineVar()) islv[vi] = true;
			else islv[vi] = false;
			mval[vi] = v.missingvalue(double(0));
			efmt[vi] = v.defaultexportformat();

			int bands = (int) v.nbands();
			I.addfield(v.getName(), efmt[vi].form, efmt[vi].width, efmt[vi].decimals, bands);

			std::string units = v.getUnits();
			if (units != "1") I.setunits(units);

			std::string desc = v.getDescription();
			desc = v.getStringAtt(AN_LONG_NAME);
			I.setdescription(desc);
		}
		I.write_aseggdf_header(dfnfilepath);
		
		cStopWatch sw;
		const size_t nl = nlines();
		for (size_t li = 0; li < nl; li++){
			std::cout << "Exporting line " << line_number[li] << std::endl;
			const size_t ns = line_index_count[li];

			std::vector<Marray<double>> A(nvars);			
			for (size_t vi = 0; vi < nvars; vi++){
				const cGeophysicsVar& v = vars[vi];								
				v.getLine(li, A[vi]);				
			}
						
			for (size_t si = 0; si < ns; si += 100){
				for (size_t vi = 0; vi < nvars; vi++){
					const cGeophysicsVar& v = vars[vi];	
					of << std::setw(efmt[vi].width);
					of << std::setprecision(efmt[vi].decimals);

					const Marray<double>& a = A[vi];
					double val;
					const size_t nb = v.nbands();
					
					double *b;					
					if (islv[vi]) b = &(a(0));
					else          b = &(a(si));
					for (size_t bi = 0; bi < nb; bi++){
						val = b[bi];
						if (val == mval[vi]) val = efmt[vi].nullvalue;						
						of << val;
					}
				}
				of << std::endl;
			}
		}
		sw.reportnow();
		prompttocontinue();
		return true;
	};
	
};

inline size_t cGeophysicsVar::line_index_start(const size_t& index) const {
	size_t start = FilePtr->get_line_index_start(index);	
	return start;
}

inline size_t cGeophysicsVar::line_index_count(const size_t& index) const {
	size_t count = FilePtr->get_line_index_count(index);
	return count;
}



#endif

