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

#define AN_STANDARD_NAME "standard_name"
#define SN_LINE_NUMBER   "line_number"
#define SN_SAMPLE_NUMBER "point_number"
#define SN_FLIGHT_NUMBER "flight_number"

#define AN_UNITS "units"
#define AN_DESCRIPTION "description"
#define AN_MISSINGVALUE "_FillValue"
#define AN_ORIGINAL_NAME "original_database_name"

#define NcShortNull -32767
#define NcIntNull -2147483647

//#define NcFloatNull  9.969209968386869e+36F
//#define NcDoubleNull 9.969209968386869e+36

#define NcFloatNull  -9999.0F
#define NcDoubleNull 9.969209968386869e+36


//#define NcFloatNull  -3.4E+38F
//#define NcDoubleNull -5.0E+75
//#define NcFloatNull  (float)-std::pow(2,128)
//#define NcDoubleNull -std::pow(2,1024)

/*
netCDF4.default_fillvals
{ 'S1': '\x00',
'f4' : 9.969209968386869e+36,
'f8' : 9.969209968386869e+36,
'i1' : -127,
'i2' : -32767,
'i4' : -2147483647,
'i8' : -9223372036854775806,
'u1' : 255,
'u2' : 65535,
'u4' : 4294967295L,
'u8' : 18446744073709551614L }
*/

std::string errormsg(const char* file, const int& linenumber, const char* function);

NcType nctype(const short dummy);
NcType nctype(const int dummy);
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
		
};

class cSampleVar : public cGeophysicsVar{

private:

public:

	cSampleVar(cGeophysicsNcFile* _parent, const NcVar& var) 
		: cGeophysicsVar(_parent, var) {		
	}

	template<typename T>
	bool putAll(std::vector<T> vals){

		if (isNull()){
			std::string msg = strprint("Attempt to write to a Null variable\n");
			msg += errormsg(__FILE__, __LINE__, __FUNCTION__);
			throw(std::runtime_error(msg.c_str()));
		}

		if (vals.size() != length()){
			std::string msg = strprint("Attempt to write variable (%s) with non-matching size\n", getName().c_str());
			msg += errormsg(__FILE__, __LINE__, __FUNCTION__);
			throw(std::runtime_error(msg));
		}

		putVar(vals.data());
		return true;
	}

	template<typename T>
	bool putLine(std::vector<T> vals, size_t lineindex, size_t bandindex = 0){
		if (isNull()){ return false; }		
		std::vector<size_t> startp = { line_index_start(lineindex), bandindex };
		std::vector<size_t> countp = { line_index_count(lineindex), 1 };
		putVar(startp, countp, vals.data());
		return true;
	}

	template<typename T>
	bool getLine(std::vector<T>& vals, const size_t lineindex, const size_t bandindex = 0){
		if (isNull()){ return false; }
		std::vector<size_t> startp = { line_index_start(lineindex), bandindex };
		std::vector<size_t> countp = { line_index_count(lineindex), 1 };
		vals.resize(countp[0]);
		getVar(startp, countp, vals.data());
		return true;
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
			std::string msg = strprint("Attempt to write to a Null variable\n");
			msg += errormsg(__FILE__, __LINE__, __FUNCTION__);
			throw(std::runtime_error(msg));
		}

		if (vals.size() != length()){
			std::string msg = strprint("Attempt to write variable (%s) with non-matching size\n", getName().c_str());
			msg += errormsg(__FILE__, __LINE__, __FUNCTION__);
			throw(std::runtime_error(msg.c_str()));
		}

		putVar(vals.data());
		return true;
	}
};

class cGeophysicsNcFile : public NcFile {

private:
	const std::string classname = "cGeophysicsNcFile";		
	std::vector<size_t> line_index_start;
	std::vector<size_t> line_index_count;
	std::vector<size_t> line_number;
		
	NcDim dim_sample() { return getDim(DN_POINT); }
	
	NcDim dim_line() { return getDim(DN_LINE); }

	bool InitialiseExisting(){
		bool status;
		status = readLineIndex();
		status = getLineNumbers(line_number);		
		return true;
	}
	
	bool readLineIndex(){

		bool status;
		size_t ns = getDim(DN_POINT).getSize();
		cLineVar v = getLineVar(VN_LI_START);
		status = v.getAll(line_index_start);
		//status = getVarAll(VN_LI_START, line_index_start);
		
		size_t nl = nlines();
		line_index_count.resize(nl);
		for (size_t i = 0; i < nl - 1; i++){
			line_index_count[i] = line_index_start[i + 1] - line_index_start[i];
		}
		line_index_count[nl - 1] = ns - line_index_start[nl - 1];
		return true;
	}

public:

	size_t get_line_index_start(const size_t& i){ return line_index_start[i];}
	size_t get_line_index_count(const size_t& i){ return line_index_count[i];}	

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

	//Create new file
	//cGeophysicsNcFile(const std::string& ncpath, const std::vector<size_t>& linenumbers, const std::vector<size_t>& nsamples)
	//	: NcFile(ncpath, FileMode::newFile, NcFile::FileFormat::nc4)
	//{
	//	InitialiseNew(linenumbers,nsamples);
	//};

	bool InitialiseNew(const std::vector<size_t>& linenumbers, const std::vector<size_t>& linesamplecount){
		
		bool status;
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

		std::vector<size_t> sample = increment(ntotalsamples(), (size_t)0, (size_t)1);
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

	NcVar findVarByAttribute(const std::string& attribute_name, const std::string& attribute_value){

		std::multimap<std::string, NcVar> vars = getVars();
		for (auto vit = vars.begin(); vit != vars.end(); vit++){
			std::map<std::string, NcVarAtt> atts = vit->second.getAtts();
			for (auto ait = atts.begin(); ait != atts.end(); ait++){
				std::string attname = ait->second.getName();
				if (attname == attribute_name){
					std::string attvalue;
					ait->second.getValues(attvalue);
					if (attvalue == attribute_value){
						return vit->second;						
					}
				}
			}
		}
		return NcVar();
	}

	template<typename T>
	bool getLineNumbers(std::vector<T>& vals){		
		NcVar v = findVarByAttribute(AN_STANDARD_NAME, SN_LINE_NUMBER);
		if (!v.isNull()){
			cLineVar var(this, v);
			return var.getAll(vals);
		}
		return false;
	}

	template<typename T>
	bool getFlightNumbers(std::vector<T>& vals){
		NcVar v = findVarByAttribute(AN_STANDARD_NAME, SN_FLIGHT_NUMBER);
		if (!v.isNull()){
			cLineVar var(this, v);
			return var.getAll(vals);
		}
		return false;		
	}
	
	template<typename T>
	bool getVarByLineNumber(const std::string& varname, const size_t& linenumber, std::vector<T>& vals){
		size_t index = getLineIndex(linenumber);
		if (index >= nlines())return false;
		return getVarByLineIndex(varname, index, vals);
	}

	template<typename T>
	bool getVarByLineIndex(const std::string& varname, const size_t& lineindex, std::vector<T>& vals){
		NcVar var = NcFile::getVar(varname);
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
	bool getVarByLineIndex(const std::string& varname, const size_t& lineindex, std::vector<std::vector<T>>& vals){
		NcVar var = NcFile::getVar(varname);
		std::vector<NcDim> dims = var.getDims();
		size_t nd = dims.size();
		if (nd != 2) return false;

		
		size_t nsamples = line_index_count[lineindex];
		size_t nbands   = dims[2].getSize();

		std::vector<size_t> start(2);
		std::vector<size_t> count(2);
		start[0] = line_index_start[lineindex];
		count[0] = nsamples;

		vals.resize(nbands);
		for(size_t i = 0; i < nbands; i++){
			start[1] = 0;
			count[1] = nbands;
			vals[i].resize(nsamples);
			var.getVar(start, count, vals[i].data());
		}		
		return true;
	}
	
	template<typename T>
	NcDim addDimVar(const std::string& dimname, const std::vector<T> dimvals){

		size_t dimsize = dimvals.size();
		NcDim dim = getDim(dimname);
		if (dim.isNull()){
			dim = addDim(dimname, dimsize);			
		}
		else{
			if(dim.getSize() != dimsize){				
				std::string msg = strprint("Attempt to add new dimension (%s) with different size to the existing\n",dimname.c_str());
				msg += errormsg(__FILE__, __LINE__, __FUNCTION__);				
				throw(std::runtime_error(msg));
			}
		}
		
		NcVar var = getVar(dimname);		
		if (var.isNull()){
			var = addVar(dimname, nctype(dimvals[0]), dim);
		}
		else{
			if (var.getDim(0).getSize() != dimsize){
				std::string msg = strprint("Attempt to add new dimension (%s) with different size to the existing\n", dimname.c_str());
				msg += errormsg(__FILE__, __LINE__, __FUNCTION__);
				throw(std::runtime_error(msg));
			}
		}
		var.putVar(dimvals.data());
		return dim;
	}

	cSampleVar getSampleVar(const std::string& name){
		if (getVarCount() == 0) return cSampleVar(this, NcVar());
		return cSampleVar(this,getVar(name));
	}

	cLineVar getLineVar(const std::string& name){
		if (getVarCount() == 0) return cLineVar(this,NcVar());
		return cLineVar(this,getVar(name));		
	}
	
	cSampleVar addSampleVar(const std::string& name, const NcType& type, const std::vector<NcDim>& dims){

		cSampleVar var = getSampleVar(name);
		if (var.isNull()){
			std::vector<NcDim> vardims = { dim_sample() };
			for (size_t i = 0; i<dims.size(); i++){
				vardims.push_back(dims[i]);
			}
			var = cSampleVar(this,addVar(name, type, vardims));
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
			for (size_t i = 0; i<dims.size(); i++){
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

	bool findNonNullLineStartEndPoints(std::vector<double>& x1, std::vector<double>& x2, std::vector<double>& y1, std::vector<double>& y2){
		x1.resize(nlines());
		x2.resize(nlines());
		y1.resize(nlines());
		y2.resize(nlines());

		cSampleVar vx = getSampleVar("longitude");
		cSampleVar vy = getSampleVar("latitude");				
		double nvx = vx.missingvalue(nvx);
		double nvy = vy.missingvalue(nvy);
		for (size_t li = 0; li < nlines(); li++){
			std::vector<double> x;
			std::vector<double> y;
			vx.getLine(x, li);
			vy.getLine(y, li);
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
			for (size_t si = ns-1; si >= 0; si--){
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

	bool addLineStartEndPoints(){

		std::vector<double> x1;
		std::vector<double> x2;
		std::vector<double> y1;
		std::vector<double> y2;
		findNonNullLineStartEndPoints(x1, x2, y1, y2);

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

	
	bool addAlphaShapePolygon(){
		std::vector<double> x;
		std::vector<double> y;
		double maxdistance;
		std::vector<double> px;
		std::vector<double> py;		
		cSampleVar vx = getSampleVar("longitude");
		cSampleVar vy = getSampleVar("latitude");
		vx.getAll(x);
		vy.getAll(y);
		double nullx = vx.missingvalue(nullx);
		double nully = vy.missingvalue(nully);

		bool status = line_data_alpha_shape_polygon_ch(
			line_index_start, line_index_count,
			x, y, nullx, nully, 64, px, py);
		
		size_t nv = px.size();
		std::vector<double> poly(nv*2);
		for (size_t i = 0; i < nv; i++){
			poly[i*2] = px[i];
			poly[i*2 + 1] = py[i];
		}
			
		std::vector<NcDim> dims;
		dims.push_back(addDim("polygonvertex",nv));
		dims.push_back(addDim("polygonordinate",2));
	
		cGeophysicsVar v(this,addVar("bounding_polygon", ncDouble, dims));
		v.add_standard_name("bounding_polygon");
		v.add_description("bounding polygon of survey");
		v.add_units("degree");
		v.putVar(poly.data());		
		return true;
	}
	

};

#endif

