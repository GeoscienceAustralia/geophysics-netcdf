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

#include <stdexcept>
#include "float.h"
#include "general_utils.h"
#include "vector_utils.h"
#include "crs.h"
#include "cgal_utils.h"

#include <netcdf>
using namespace netCDF;
using namespace netCDF::exceptions;

#define DN_POINT "point"
#define DN_LINE  "line"
#define VN_LI_START "index_line"
#define VN_LI_COUNT "index_count"

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

class cGeophysicsNcFile;

class cGeophysicsVar : public NcVar{

private:	

	cGeophysicsNcFile* parent=NULL;

public:
	
	cGeophysicsVar(cGeophysicsNcFile* _parent, const NcVar& var) : NcVar(var) {
		parent = _parent;
	}

	cGeophysicsNcFile* get_parent(){ return parent; }

	size_t line_index_start(const size_t& index);
	size_t line_index_count(const size_t& index);
	
	size_t length(){		
		std::vector<NcDim> dims = getDims();
		size_t len = dims[0].getSize();
		for (size_t di = 1; di < dims.size(); di++){
			len *= dims[di].getSize();
		}
		return len;
	}

	size_t nbands(){
		std::vector<NcDim> dims = getDims();
		if (dims.size() == 1) return 1;
		else if (dims.size() == 2) return dims[1].getSize();
		return 0;
		//this need fixing for case of extra dims
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
	
	static float preferred_float_missing_value(){		
		return NcFloatNull;
	}

	static double preferred_double_missing_value(){		
		return NcDoubleNull;
	}

	double lowest_possible_value(){
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
			putAtt(AN_MISSINGVALUE, ncShort, (short) NcShortNull);
		}
		else if (type == ncInt){
			putAtt(AN_MISSINGVALUE, ncInt, (int) NcIntNull);
		}
		else if (type == ncFloat){
			putAtt(AN_MISSINGVALUE, ncFloat, (float) NcFloatNull);
		}
		else if (type == ncDouble){
			putAtt(AN_MISSINGVALUE, ncDouble, (double) NcDoubleNull);
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
	T missingvalue(const T& value){
		T v;
		getAtt(AN_MISSINGVALUE).getValues(&v);
		return v;
	}

	template<typename T>
	bool getAll(std::vector<T>& vals){				
		vals.resize(length());
		getVar(vals.data());
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
	bool getLine(const size_t& lineindex, T& val){
		return getLine(lineindex, 0, val);
	}

	template<typename T>
	bool getLine(const size_t& lineindex, const size_t& bandindex, std::vector<T>& vals){
		if (isNull()){ return false; }
		std::vector<size_t> startp = { lineindex, bandindex };
		std::vector<size_t> countp = { 1, 1 };
		vals.resize(countp[0]*countp[1]);
		getVar(startp, countp, vals.data());
		return true;
	}

	template<typename T>
	bool getLine(const size_t& lineindex, std::vector<T>& vals){
		return getLine(lineindex, 0, vals);
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
};

class cSampleVar : public cGeophysicsVar{

private:

public:

	cSampleVar(cGeophysicsNcFile* _parent, const NcVar& var) 
		: cGeophysicsVar(_parent, var) {		
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
	bool putLine(const size_t& lineindex, const size_t& bandindex, const std::vector<T>& vals){
		if (isNull()){ return false; }		
		std::vector<size_t> start = { line_index_start(lineindex), bandindex };
		std::vector<size_t> count = { line_index_count(lineindex), 1 };
		putVar(start, count, vals.data());
		return true;
	}

	template<typename T>
	bool putLine(const size_t& lineindex, const std::vector<T>& vals){		
		return putLine(lineindex, 0, vals);	
	}

	template<typename T>
	bool getLine(const size_t& lineindex, const size_t& bandindex, std::vector<T>& vals){
		if (isNull()){ return false; }
		std::vector<size_t> start = { line_index_start(lineindex), bandindex };
		std::vector<size_t> count = { line_index_count(lineindex), 1 };
		vals.resize(count[0]*count[1]);
		getVar(start, count, vals.data());
		return true;
	}

	template<typename T>
	bool getLine(const size_t& lineindex, std::vector<T>& vals){		
		return getLine(lineindex, 0, vals);
	}

};

class cLineVar : public cGeophysicsVar{

private:

public:

	cLineVar(cGeophysicsNcFile* _parent, const NcVar& var)
		: cGeophysicsVar(_parent, var) {
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
};

class cGeophysicsNcFile : public NcFile {

private:
	const std::string classname = "cGeophysicsNcFile";	
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

		cLineVar vs = getLineVar(VN_LI_START);
		if (vs.getAll(line_index_start) == false)return false;
		
		cLineVar vc = getLineVar(VN_LI_COUNT);
		if(vc.getAll(line_index_count) == false)return false;
		
		return true;
	}

public:	

	size_t get_line_index_start(const size_t& i){ return line_index_start[i]; }
	size_t get_line_index_count(const size_t& i){ return line_index_count[i]; }

	//Open existing file
	cGeophysicsNcFile(const std::string& ncpath, const FileMode& filemode)
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

	~cGeophysicsNcFile()
	{
		
	};

	bool InitialiseNew(const std::vector<size_t>& linenumbers, const std::vector<size_t>& linesamplecount){
		size_t n = linenumbers.size();
		std::vector<unsigned int> uint_linenumbers(n);
		std::vector<unsigned int> uint_linesamplecount(n);
		for (size_t i = 0; i < n; i++){
			uint_linenumbers[i] = linenumbers[i];
			uint_linesamplecount[i] = linesamplecount[i];
		}
		return InitialiseNew(uint_linenumbers, uint_linesamplecount);
	}

	bool InitialiseNew(const std::vector<unsigned int>& linenumbers, const std::vector<unsigned int>& linesamplecount){

		size_t nl = linenumbers.size();
		line_number = linenumbers;
		line_index_count = linesamplecount;
		line_index_start.resize(nl);
		size_t nsamples = 0;
		for (size_t i = 0; i < line_number.size(); i++){
			line_index_start[i] = nsamples;
			nsamples += line_index_count[i];
		}

		NcDim ds = addDim(DN_POINT, nsamples);
		NcDim dl = addDim(DN_LINE, nl);

		cLineVar vstart = addLineVar(VN_LI_START, ncUint);
		vstart.putVar(line_index_start.data());
		vstart.add_standard_name(VN_LI_START);		
		vstart.add_description("zero based index of the first sample in the line");
		vstart.add_units("1");

		cLineVar vcount = addLineVar(VN_LI_COUNT, ncUint);
		vcount.putVar(line_index_count.data());
		vcount.add_standard_name(VN_LI_COUNT);
		vcount.add_description("number of samples in the line");
		vcount.add_units("1");

		cLineVar vline = addLineVar(DN_LINE, ncUint);
		vline.putVar(line_number.data());
		vline.add_standard_name(SN_LINE_NUMBER);
		vline.add_description("flight line number");
		vline.add_units("1");

		std::vector<unsigned int> sample = increment((unsigned int)ntotalsamples(), (unsigned int)0, (unsigned int)1);
		cSampleVar vsample = addSampleVar(DN_POINT, ncUint);
		vsample.putVar(sample.data());
		vsample.add_standard_name(SN_SAMPLE_NUMBER);
		vsample.add_description("sequential point number");
		vsample.add_units("1");
		return true;
	}

	size_t nlines(){ return line_index_start.size(); }
	size_t ntotalsamples(){ return sum(line_index_count); }
	size_t nlinesamples(const size_t lineindex){ return line_index_count[lineindex]; }

	size_t getLineIndex(const int& linenumber){
		auto it = std::find(line_number.begin(), line_number.end(), linenumber);
		return it - line_number.begin();
	}

	bool hasVar(const std::string& varname){
		NcVar v = getVar(varname);
		if (v.isNull()) return false;
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
		if (var.isNull()){
			std::vector<NcDim> vardims = { dim_sample() };
			for (size_t i = 0; i < dims.size(); i++){
				vardims.push_back(dims[i]);
			}
			var = cSampleVar(this, addVar(name, type, vardims));
			//var.set_parent(this);
			var.set_default_missingvalue();
		}
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
		if (var.isNull()){
			std::vector<NcDim> vardims = { dim_line() };
			for (size_t i = 0; i < dims.size(); i++){
				vardims.push_back(dims[i]);
			}
			var = cLineVar(this, addVar(name, type, vardims));
			var.set_default_missingvalue();
		}
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
			for (size_t si = ns - 1; si >= 0; si--){
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
};

#endif

