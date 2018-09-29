/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2016.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#ifndef _geophysics_netcdf_H
#define _geophysics_netcdf_H

#ifdef USEGLOBALSTACKTRACE
#include "stacktrace.h"
#endif

#include <cassert>
#include <stdexcept>
#include <map>
#include <iomanip> 
#include <memory> 
#include "float.h"
#include "general_utils.h"
#include "vector_utils.h"
#include "crs.h"
#include "cgal_utils.h"
#include "file_formats.h"
#include "stopwatch.h"

#include <netcdf>
using namespace netCDF;
using namespace netCDF::exceptions;

#include "marray.hxx"
using namespace andres;

#define DN_POINT "point"
#define DN_LINE  "line"

#define VN_LI_START "index_line"
#define VN_LI_COUNT "index_count"
#define VN_LINE_INDEX "line_index"

#define SN_LINE_NUMBER   "line_number"
#define SN_SAMPLE_NUMBER "point_number"
#define SN_FLIGHT_NUMBER "flight_number"

#define AN_STANDARD_NAME "standard_name"
#define AN_UNITS "units"
#define AN_DESCRIPTION "description"
#define AN_MISSINGVALUE "_FillValue"
#define AN_ORIGINAL_NAME "original_database_name"

#define NcShortNull -32767
#define NcIntNull -2147483647
#define NcFloatNull  9.969209968386869e+36F
#define NcDoubleNull 9.969209968386869e+36

NcType nctype(const short dummy);
NcType nctype(const int dummy);
NcType nctype(const unsigned int dummy);
NcType nctype(const float dummy);
NcType nctype(const double dummy);
NcType nctype(const std::string dummy);

class cExportFormat{

public:
	char   form;
	size_t width;
	size_t decimals;
	double nullvalue;

	cExportFormat(){};

	cExportFormat(char _form, int _width, int _decimals, double _nullvalue){
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

	cGeophysicsVar(cGeophysicsNcFile* parentptr, const NcVar& var) : NcVar(var) {
		FilePtr = parentptr;
	};

	cGeophysicsVar(const cGeophysicsVar& var) : NcVar(var) {
		FilePtr = var.FilePtr;
	};

	cGeophysicsVar& operator=(const cGeophysicsVar& rhs){
		NcVar::operator=(rhs);
		FilePtr = rhs.FilePtr;
		return *this;
	};

	size_t line_index_start(const size_t& index) const;
	size_t line_index_count(const size_t& index) const;

	size_t length(){
		std::vector<NcDim> dims = getDims();
		size_t len = dims[0].getSize();
		for (size_t di = 1; di < dims.size(); di++){
			len *= dims[di].getSize();
		}
		return len;
	}

	size_t sizeBytes(){
		return length() * getType().getSize();
	}

	size_t elementspersample() const {
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

	NcVarAtt add_attribute(const std::string& att, std::string value){
		return putAtt(att, value);
	}

	NcVarAtt add_standard_name(const std::string& value){
		return putAtt(AN_STANDARD_NAME, value);
	}

	NcVarAtt add_original_name(const std::string& value){
		return putAtt(AN_ORIGINAL_NAME, value);
	}

	NcVarAtt add_units(const std::string& value){
		return putAtt(AN_UNITS, value);
	}

	NcVarAtt add_description(const std::string& value){
		return putAtt(AN_DESCRIPTION, value);
	}

	bool hasAtt(const std::string& name) const {
		std::map<std::string, NcVarAtt> m = NcVar::getAtts();
		auto it = m.find(name);
		if (it == m.end()){
			return false;
		}
		return true;
	}

	std::string getStringAtt(const std::string& attname) const {
		std::string attvalue;
		if (hasAtt(attname)){
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

	static float preferred_float_missing_value() {
		return NcFloatNull;
	}

	static double preferred_double_missing_value() {
		return NcDoubleNull;
	}

	double lowest_possible_value() const {
		const NcType type = getType();
		if (type == ncShort) return SHRT_MIN;
		else if (type == ncUint) return 0;
		else if (type == ncInt) return INT_MIN;
		else if (type == ncFloat) return -FLT_MAX;
		else return -DBL_MAX;
	}

	double highest_possible_value(){
		const NcType type = getType();
		if (type == ncShort) return SHRT_MAX;
		else if (type == ncUint) return UINT_MAX;
		else if (type == ncInt) return INT_MAX;
		else if (type == ncFloat) return FLT_MAX;
		else return DBL_MAX;
	}

	bool set_default_missingvalue(){
		const NcType type = getType();
		if (type == ncShort){
			putAtt(AN_MISSINGVALUE, ncShort, (short)NcShortNull);
		}
		else if (type == ncInt){
			putAtt(AN_MISSINGVALUE, ncInt, (int)NcIntNull);
		}
		else if (type == ncFloat){
			putAtt(AN_MISSINGVALUE, ncFloat, (float)NcFloatNull);
		}
		else if (type == ncDouble){
			putAtt(AN_MISSINGVALUE, ncDouble, (double)NcDoubleNull);
		}
		else return false;
		return true;
	}

	template<typename T>
	bool add_missing_value(const T& value){
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
		
		if (hasAtt(AN_MISSINGVALUE)){
			T v;
			getAtt(AN_MISSINGVALUE).getValues(&v);
			return v;
		}		
		else{
			const NcType type = getType();
			if (type == ncShort) return (T) NcShortNull;
			else if (type == ncInt)  return (T) NcIntNull;
			else if (type == ncFloat) return (T) NcFloatNull;
			else if (type == ncDouble) return (T) NcDoubleNull;
			else return (T) NcShortNull;
		}
	}

	template<typename T>
	bool getAll(std::vector<T>& vals){
		vals.resize(length());
		getVar(vals.data());
		return true;
	}	

	template<typename T>
	bool minmax(T& minval, T& maxval){

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
	bool getLine(const size_t& lineindex, const size_t& bandindex, T& val) const{
		return false;
	}

	template<typename T>
	bool getLine(const size_t& lineindex, T& val) const = 0;

	template<typename T>
	bool getLine(const size_t& lineindex, const size_t& bandindex, std::vector<T>& vals){
		return false;
	}
	
	template<typename T>
	bool getLine(const size_t& lineindex, std::vector<T>& vals) const = 0;
		
	template<typename T>
	void getLine(const size_t& lineindex, Marray<T>& A) const {
		if (isNull()){
			std::string msg = _SRC_ + strprint("Attempt to read from a Null variable\n");
			logmsg(msg);
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
		size_t sz = lineelements(lineindex);
		getVar(start, count, &(A(0)));
	}

	bool donotexport(){
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
		
		cExportFormat e;
		nc_type t = getType().getId();
		switch (t){
		case NC_SHORT: e = cExportFormat('I',8,0, -999); break;
		case NC_INT: e = cExportFormat('I',12,0, -999); break;
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
		: cGeophysicsVar(parentptr, var) {
	}

	template<typename T>
	bool putAll(std::vector<T> vals){

		if (isNull()){
			std::string msg = _SRC_ + strprint("Attempt to write to a Null variable\n");
			logmsg(msg);
			throw(std::runtime_error(msg));
		}

		if (vals.size() != length()){
			std::string msg = _SRC_ + strprint("Attempt to write variable (%s) with non-matching size\n", getName().c_str());
			logmsg(msg);
			throw(std::runtime_error(msg));
		}

		putVar(vals.data());
		return true;
	}

	template<typename T>
	bool getLine(const size_t& lineindex, const size_t& bandindex, T& val){
		if (isNull()){ return false; }
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
		if (isNull()){ return false; }
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

	cSampleVar(cGeophysicsNcFile* parentptr, const NcVar& var) 
		: cGeophysicsVar(parentptr, var) {
	}

	template<typename T>
	bool putAll(const std::vector<T>& vals){

		if (isNull()){
			std::string msg = _SRC_ + strprint("Attempt to write to a Null variable\n");
			logmsg(msg);
			throw(std::runtime_error(msg));
		}

		if (vals.size() != length()){
			std::string msg = _SRC_ + strprint("Attempt to write variable (%s) with non-matching size\n", getName().c_str());
			logmsg(msg);
			throw(std::runtime_error(msg));
		}

		putVar(vals.data());
		return true;
	}
	
	template<typename T>
	bool putLineBand(const size_t& lineindex, const size_t& bandindex, const std::vector<T>& vals){
		
		if (isNull()){
			std::string msg = _SRC_ + strprint("Attempt to write to a Null variable\n");
			logmsg(msg);
			throw(std::runtime_error(msg));
		}

		std::vector<NcDim>  dims = getDims();
		if (dims.size() > 2){
			std::string msg = _SRC_ + strprint("Attempt to use putLineBand() to write to a variable with more than 2 dimensions\n");
			logmsg(msg);
			throw(std::runtime_error(msg));
		}

		if (vals.size() != line_index_count(lineindex)){
			std::string msg = _SRC_ + strprint("Attempt to write line/band of variable (%s) with non-matching size\n", getName().c_str());
			logmsg(msg);
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
		if (isNull()){
			std::string msg = _SRC_ + strprint("Attempt to write to a Null variable\n");
			logmsg(msg);
			throw(std::runtime_error(msg));
		}

		size_t sz = lineelements(lineindex);
		if (vals.size() != sz){
			std::string msg = _SRC_ + strprint("Attempt to write line/band of variable (%s) with non-matching size\n", getName().c_str());
			logmsg(msg);
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
		if (isNull()){
			std::string msg = _SRC_ + strprint("Attempt to read from a Null variable\n");
			logmsg(msg);
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
		if (isNull()){
			std::string msg = _SRC_ + strprint("Attempt to read from a Null variable\n");
			logmsg(msg);
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
	T getSample(const size_t& lineindex, const size_t& sampleindex, const size_t& bandindex, T& val) const {
		if (isNull()){ return false; }
		std::vector<size_t> startp = { lineindex+sampleindex, bandindex };
		std::vector<size_t> countp = { 1, 1 };	
		getVar(startp, countp, &val);
		return val;
	};

};

class cGeophysicsNcFile : public NcFile {

private:
		
	std::vector<unsigned int> line_index_start;
	std::vector<unsigned int> line_index_count;
	std::vector<unsigned int> line_number;

	NcDim dim_sample() { return getDim(DN_POINT); }

	NcDim dim_line() { return getDim(DN_LINE); }

	bool InitialiseExisting(){		
		if (readLineIndex() == false) return false;
		if (getLineNumbers(line_number) == false) return false;
		return true;		
	}

	bool readLineIndex(){

		if (hasVar(VN_LI_COUNT)){
			cLineVar vc = getLineVar(VN_LI_COUNT);
			if (vc.getAll(line_index_count) == false)return false;

			cLineVar vs = getLineVar(VN_LI_START);			
			if (vs.getAll(line_index_start) == false)return false;
		}
		else if (hasVar(VN_LINE_INDEX)){
			cLineVar vl = getLineVar(VN_LINE_INDEX);
			std::vector<unsigned int> line_index;
			if (vl.getAll(line_index) == false)return false;
			set_start_count(line_index);
		}				
		return true;
	}

	static std::vector<unsigned int> compute_line_index(std::vector<unsigned int>& count)
	{
		std::vector<unsigned int> index(sum(count));
		size_t k = 0;		
		for (size_t i = 0; i < count.size(); i++){						
			for (size_t j = 0; j < count[i]; j++){
				index[k] = (unsigned int)i;
				k++;
			}
		}
		return index;
	}

	void set_start_count(const std::vector<unsigned int>& line_index)
	{
		line_index_start = compute_start_from_line_index(line_index);		
		line_index_count = compute_count_from_start(line_index_start,line_index.size());
		return;
	}

	static std::vector<unsigned int> compute_start_from_line_index(const std::vector<unsigned int>& line_index)
	{
		size_t ns = line_index.size();
		std::vector<unsigned int> start;		
		start.push_back(0);
		for (size_t i = 1; i<ns; i++){
			if (line_index[i] != line_index[i - 1]){
				start.push_back((unsigned int)i);
			}
		}		
		return start;
	}

	static std::vector<unsigned int> compute_start_from_count(const std::vector<unsigned int>& count)
	{
		size_t nl = count.size();
		std::vector<unsigned int> start(nl);
		start[0] = 0;
		for (size_t i = 1; i<nl; i++){
			start[i] = start[i-1] + count[i];
		}		
		return start;
	}

	static std::vector<unsigned int> compute_count_from_start(const std::vector<unsigned int>& start, size_t nsamplestotal)
	{
		size_t nl = start.size();
		std::vector<unsigned int> count(nl);
		for (size_t i = 0; i<nl - 1; i++){
			count[i] = start[i + 1] - start[i];
		}
		count[nl - 1] = nsamplestotal - start[nl - 1];		
		return count;
	}

	static size_t nelements(const NcVar& v)
	{
		std::vector<NcDim> dims = v.getDims();
		size_t n = 0;
		for (size_t di = 0; di < dims.size(); di++){
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
	cGeophysicsNcFile& operator=(const cGeophysicsNcFile & rhs);
	cGeophysicsNcFile& operator=(const NcGroup & rhs);
	cGeophysicsNcFile& operator=(const NcFile & rhs);
	cGeophysicsNcFile(const cGeophysicsNcFile& rhs);
	cGeophysicsNcFile(const NcGroup& rhs);
	cGeophysicsNcFile(const NcFile& rhs);

	size_t get_line_index_start(const size_t& i) const { return line_index_start[i]; }
	size_t get_line_index_count(const size_t& i) const { return line_index_count[i]; }
	
	//Open existing file constructor
	cGeophysicsNcFile(const std::string& ncpath, const FileMode& filemode = NcFile::FileMode::read)
		: NcFile(ncpath, filemode)
	{
		if (filemode == NcFile::read){
			InitialiseExisting();
		}
		else if (filemode == NcFile::write){
			InitialiseExisting();
		}
		else if (filemode == NcFile::replace){

		}
		else if (filemode == NcFile::newFile){

		}
		else{

		}
	};

	//Convert legacy file constructor
	cGeophysicsNcFile(const std::string& filename, const cGeophysicsNcFile& srcfile)
		: NcFile(filename, NcFile::FileMode::newFile)
	{
		bool status  = convert_legacy(srcfile);
	};
	
	//Destructor
	~cGeophysicsNcFile(){};	

	bool InitialiseNew(const std::vector<size_t>& linenumbers, const std::vector<size_t>& linesamplecount){
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
		size_t nl = linenumbers.size();
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
		cSampleVar vindex = addSampleVar(VN_LINE_INDEX, ncUint);
		vindex.putVar(line_index.data());
		vindex.add_standard_name(VN_LINE_INDEX);
		vindex.add_description("zero based index of line associated with point");
		vindex.add_units("1");
	}

	//Deprecated
	void add_line_index_start(const std::vector<unsigned int> line_index_start)
	{
		cLineVar vstart = addLineVar(VN_LI_START, ncUint);		
		vstart.putVar(line_index_start.data());
		vstart.add_standard_name(VN_LI_START);		
		vstart.add_description("zero based index of the first sample in the line");
		vstart.add_units("1");
	}

	//Deprecated
	void add_line_index_count(const std::vector<unsigned int> line_index_count)
	{
		cLineVar vcount = addLineVar(VN_LI_COUNT, ncUint);
		vcount.putVar(line_index_count.data());
		vcount.add_standard_name(VN_LI_COUNT);
		vcount.add_description("number of samples in the line");
		vcount.add_units("1");
	}

	//Deprecated
	void add_point_variable()
	{
		std::vector<unsigned int> sample = increment((unsigned int)ntotalsamples(), (unsigned int)0, (unsigned int)1);
		cSampleVar v = addSampleVar(DN_POINT, ncUint);		
		v.putVar(sample.data());
		v.add_standard_name(SN_SAMPLE_NUMBER);
		v.add_description("sequential point number");
		v.add_units("1");
	}

	void add_line_number_variable(const std::vector<unsigned int> linenumbers)
	{
		cLineVar vline = addLineVar(DN_LINE, ncUint);
		vline.putVar(linenumbers.data());
		vline.add_standard_name(SN_LINE_NUMBER);
		vline.add_description("flight line number");
		vline.add_units("1");		
	}

	
	//Convert legacy file
	static bool convert_legacy(const std::string& srcpathname, const std::string& dstpathname)
	{
		cGeophysicsNcFile srcfile(srcpathname);
		cGeophysicsNcFile dstfile(dstpathname, NcFile::FileMode::newFile);
		bool status = dstfile.convert_legacy(srcfile);
		return status;
	}

	//Convert legacy file
	bool convert_legacy(const cGeophysicsNcFile& srcfile)
	{
		//bookmark
		InitialiseNew(srcfile.line_number, srcfile.line_index_count);
		copy_global_atts(srcfile);
		copy_dims(srcfile);
		auto vm = srcfile.getVars();
		for (auto vit = vm.begin(); vit != vm.end(); vit++){
			NcVar& srcvar = vit->second;			
			if (srcvar.getName() == VN_LI_START) continue;
			if (srcvar.getName() == VN_LI_COUNT) continue;			
			copy_var(srcvar);
		}		
		return true;
	}

	bool copy_global_atts(const cGeophysicsNcFile& srcfile)
	{
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
		auto dm = srcfile.getDims();
		for (auto dit = dm.begin(); dit != dm.end(); dit++){
			NcDim& srcdim = dit->second;
			if (hasDim(srcdim.getName())) continue;
			addDim(srcdim.getName(), srcdim.getSize());
		}
		return true;
	}

	bool copy_var(const NcVar& srcvar, std::string newname=std::string())
	{				
		std::string name = newname;
		if (name.size() == 0){
			name = srcvar.getName();
		}

		if (hasVar(name)) return false;
		NcVar v = addVar(name, srcvar.getType(), srcvar.getDims());
		copy_varatts(srcvar, v);

		size_t len = nbytes(srcvar);
		if (len > 0){
			std::vector<uint8_t> buf(len);
			srcvar.getVar((void*)buf.data());
			v.putVar((void*)buf.data());
		}
		return true;
	}
		
	bool copy_varatts(const NcVar& srcvar, const NcVar& dstvar)
	{
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

	size_t nlines(){ return line_index_start.size(); }
	size_t ntotalsamples(){ return sum(line_index_count); }
	size_t nlinesamples(const size_t lineindex){ return line_index_count[lineindex]; }

	size_t getLineIndex(const int& linenumber){
		auto it = std::find(line_number.begin(), line_number.end(), linenumber);
		return (size_t)(it - line_number.begin());
	}

	bool hasVar(const std::string& varname){
		NcVar v = getVar(varname);
		if (v.isNull()) return false;
		return true;
	}

	bool hasDim(const std::string& dimname){
		NcDim d = getDim(dimname);
		if (d.isNull()) return false;
		return true;
	}

	NcVar getVarByAtt(const std::string& att_name, const std::string& att_value){
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

	NcVar getVarByStandardName(const std::string& att_value){
		return getVarByAtt(AN_STANDARD_NAME, att_value);
	}

	std::string getVarNameByAtt(const std::string& att_name, const std::string& att_value){
		NcVar var = getVarByAtt(att_name, att_value);
		if (var.isNull()){
			return std::string();
		}
		else return var.getName();
	}

	std::string getVarNameByStandardName(const std::string& att_value){
		return getVarNameByAtt(AN_STANDARD_NAME, att_value);
	}

	template<typename T>
	bool getLineNumbers(std::vector<T>& vals){
		NcVar v = getVarByStandardName(SN_LINE_NUMBER);
		if (v.isNull() == false){
			cLineVar var(this, v);
			return var.getAll(vals);
		}
		return false;
	}
	
	std::vector<int> getLineNumbers(){
		std::vector<int> num;
		getLineNumbers(num);
		return num;
	}

	template<typename T>
	bool getFlightNumbers(std::vector<T>& vals){
		NcVar v = getVarByStandardName(SN_FLIGHT_NUMBER);
		if (v.isNull() == false){
			cLineVar var(this, v);
			return var.getAll(vals);
		}
		return false;
	}

	std::vector<int> getFlightNumbers(){
		std::vector<int> num;
		getFlightNumbers(num);
		return num;
	}

	template<typename T>
	bool getDataByLineIndex(const cSampleVar& var, const size_t& lineindex, std::vector<T>& vals){		
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
	NcDim addDimVar(const std::string& dimname, const std::vector<T>& dimvals){

		size_t dimsize = dimvals.size();
		NcDim dim = getDim(dimname);
		if (dim.isNull()){
			dim = addDim(dimname, dimsize);
		}
		else{
			if (dim.getSize() != dimsize){
				std::string msg = _SRC_ + strprint("Attempt to add new dimension (%s) with different size to the existing\n", dimname.c_str());
				logmsg(msg);
				throw(std::runtime_error(msg));
			}
		}

		NcVar var = getVar(dimname);
		if (var.isNull()){
			var = addVar(dimname, nctype(dimvals[0]), dim);
		}
		else{
			if (var.getDim(0).getSize() != dimsize){
				std::string msg = _SRC_ + strprint("Attempt to add new dimension (%s) with different size to the existing\n", dimname.c_str());
				logmsg(msg);
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
		std::vector<double>& x1, std::vector<double>& x2,
		std::vector<double>& y1, std::vector<double>& y2){
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

	bool addCRS(const cCRS& crs){

		#ifdef _WIN32 
			NcVar v = NcFile::addVar("crs", ncInt);
		#else 			
			std::vector<NcDim> d;
			NcVar v = NcFile::addVar("crs", ncInt, d);
		#endif
		v.putAtt("grid_mapping_name", "latitude_longitude");
		v.putAtt("epsg_code", crs.epsg_code.c_str());
		v.putAtt("semi_major_axis", ncDouble, crs.semi_major_axis);
		v.putAtt("inverse_flattening", ncDouble, crs.inverse_flattening);
		return true;
	}

	bool addLineStartEndPointsLL(){

		if (hasVar("longitude_first")){
			std::string msg = _SRC_ + strprint("Warning: Variable longitude_first already exists\n");
			logmsg(msg);			
			return false;
		}

		if (hasVar("latitude_first")){
			std::string msg = _SRC_ + strprint("Warning: Variable latitude_first already exists\n");
			logmsg(msg);
			return false;
		}

		std::vector<double> x1;
		std::vector<double> x2;
		std::vector<double> y1;
		std::vector<double> y2;
		findNonNullLineStartEndPoints("longitude", "latitude", x1, x2, y1, y2);

		cLineVar vx1 = addLineVar("longitude_first", ncDouble);
		vx1.add_missing_value(vx1.preferred_double_missing_value());
		vx1.add_standard_name("longitude_first");
		vx1.add_description("first non-null longitude coordinate in the line");
		vx1.add_units("degree_east");
		vx1.putAll(x1);

		cLineVar vx2 = addLineVar("longitude_last", ncDouble);
		vx2.add_missing_value(vx2.preferred_double_missing_value());
		vx2.add_standard_name("longitude_last");
		vx2.add_description("last non-null longitude coordinate in the line");
		vx2.add_units("degree_east");
		vx2.putAll(x2);

		cLineVar vy1 = addLineVar("latitude_first", ncDouble);
		vy1.add_missing_value(vy1.preferred_double_missing_value());
		vy1.add_standard_name("latitude_first");
		vy1.add_description("first non-null latitude coordinate in the line");
		vy1.add_units("degree_north");
		vy1.putAll(y1);

		cLineVar vy2 = addLineVar("latitude_last", ncDouble);
		vy2.add_missing_value(vy2.preferred_double_missing_value());
		vy2.add_standard_name("latitude_last");
		vy2.add_description("last non-null latitude coordinate in the line");
		vy2.add_units("degree_north");
		vy2.putAll(y2);
		return true;
	}

	bool addLineStartEndPointsEN(){

		if (hasVar("easting_first")){
			std::string msg = _SRC_ + strprint("Warning: Variable easting_first already exists\n");
			logmsg(msg);
			return false;
		}

		if (hasVar("northing_first")){
			std::string msg = _SRC_ + strprint("Warning: Variable northing_first already exists\n");
			logmsg(msg);
			return false;
		}

		std::vector<double> x1;
		std::vector<double> x2;
		std::vector<double> y1;
		std::vector<double> y2;
		findNonNullLineStartEndPoints("easting", "northing", x1, x2, y1, y2);
				
		

		cLineVar vx1 = addLineVar("easting_first", ncDouble);
		vx1.add_missing_value(vx1.preferred_double_missing_value());
		vx1.add_standard_name("easting_first");
		vx1.add_description("first non-null easting coordinate in the line");
		vx1.add_units("m");
		vx1.putAll(x1);

		cLineVar vx2 = addLineVar("easting_last", ncDouble);
		vx2.add_missing_value(vx2.preferred_double_missing_value());
		vx2.add_standard_name("easting_last");
		vx2.add_description("last non-null easting coordinate in the line");
		vx2.add_units("m");
		vx2.putAll(x2);

		cLineVar vy1 = addLineVar("northing_first", ncDouble);
		vy1.add_missing_value(vy1.preferred_double_missing_value());
		vy1.add_standard_name("northing_first");
		vy1.add_description("first non-null northing coordinate in the line");
		vy1.add_units("m");
		vy1.putAll(y1);

		cLineVar vy2 = addLineVar("northing_last", ncDouble);
		vy2.add_missing_value(vy2.preferred_double_missing_value());
		vy2.add_standard_name("northing_last");
		vy2.add_description("last non-null northing coordinate in the line");
		vy2.add_units("m");
		vy2.putAll(y2);
		return true;
	}
	
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
		v.add_standard_name("bounding_polygon");
		v.add_description("bounding polygon of survey");
		v.add_units("degree");
		v.putVar(poly.data());
		return true;
	}
	
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

	bool addGeospatialMetadataXY(){
		
		bool status;
		if (hasVar("longitude") && hasVar("latitude")){
			status = addGeospatialMetadataItem("longitude", "lon", "degrees_east");
			if (status == false)return status;
			
			status = addGeospatialMetadataItem("latitude", "lat", "degrees_north");
			if (status == false)return status;
			cCRS crs("GDA94");
			addCRS(crs);
		}
		else if (hasVar("easting") && hasVar("northing")){
			status = addGeospatialMetadataItem("easting", "east", "m");
			if (status == false)return status;
			
			status = addGeospatialMetadataItem("northing", "north", "m");
			if (status == false)return status;

			cCRS crs("GDA94");
			addCRS(crs);			
		}
		else return false;

		return true;
	}

	bool addGeospatialMetadataVertical(){

		if (hasVar("height")){
			bool status = addGeospatialMetadataItem("height", "vertical", "m");
			if (status == false)return status;
		}
		else return false;

		return true;

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
		std::vector<NcDim> dims;
		std::multimap<std::string, NcDim> m = getDims();
		for (auto it = m.begin(); it != m.end(); it++){
			dims.push_back(it->second);
		}
		return dims;
	};

	std::vector<NcVar> getAllVars() {
		std::vector<NcVar> vars;
		std::multimap<std::string, NcVar> vm = getVars();
		for (auto it = vm.begin(); it != vm.end(); it++){						
			vars.push_back(it->second);
		}		
		return vars;
	};

	std::vector<cLineVar> getLineVars() {
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

			size_t bands = v.nbands();
			I.addfield(v.getName(), efmt[vi].form, efmt[vi].width, efmt[vi].decimals, bands);

			std::string units = v.getUnits();
			if (units != "1") I.setunits(units);

			std::string desc = v.getDescription();
			desc = v.getStringAtt(AN_STANDARD_NAME);
			I.setcomment(desc);
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

#endif

