/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2016.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#ifndef _geophysics_netcdf_H
#define _geophysics_netcdf_H

#include "general_utils.h"
#include "vector_utils.h"
#include "crs.h"
#include <netcdf>
using namespace netCDF;
using namespace netCDF::exceptions;

#define DN_POINT "point"
#define DN_LINE  "line"
#define VN_LI_START "_index_line"
#define VN_LI_COUNT "_index_count"

#define AN_STANDARD_NAME "standard_name"
#define SN_LINE_NUMBER   "line_number"
#define SN_SAMPLE_NUMBER "point_number"
#define SN_FLIGHT_NUMBER "flight_number"

#define AN_UNITS "units"
#define AN_DESCRIPTION "description"
#define AN_FILLVALUE "_FillValue"
#define AN_ORIGINAL_NAME "original_database_name"

#define NcShortNull -32767
#define NcIntNull -2147483647
#define NcFloatNull  -3.4E+38F
#define NcDoubleNull -5.0E+75

std::string errormsg(const char* file, const int& linenumber, const char* function);

NcType nctype(const short dummy);
NcType nctype(const int dummy);
NcType nctype(const double dummy);
NcType nctype(const std::string dummy);

class cGeophysicsNcFile;

class cGeophysicsVar : public NcVar{

private:	
	cGeophysicsNcFile* parent=NULL;

public:
	
	cGeophysicsVar(const NcVar& var) : NcVar(var) { }

	void set_parent(cGeophysicsNcFile* _parent){ parent = _parent; }
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
		else if (dims.size() == 2) dims[1].getSize();
		else return 0;
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
	
	static float preferred_float_fillvalue(){		
		return NcFloatNull;
	}

	static double preferred_double_fillvalue(){		
		return NcDoubleNull;
	}

	bool set_default_fillvalue(){
		const NcType type = getType();
		if (type == ncShort) putAtt(AN_FILLVALUE, ncShort, NcShortNull);
		else if (type == ncInt) putAtt(AN_FILLVALUE, ncInt, NcIntNull);
		else if (type == ncFloat) putAtt(AN_FILLVALUE, ncFloat, NcFloatNull);
		else if (type == ncDouble) putAtt(AN_FILLVALUE, ncDouble, NcDoubleNull);
		else return false;
		return true;
	}

	template<typename T>
	bool add_fillvalue(const T& value){
		const NcType type = getType();
		if (type == ncShort) putAtt(AN_FILLVALUE, ncShort, value);
		else if (type == ncInt) putAtt(AN_FILLVALUE, ncInt, value);
		else if (type == ncFloat) putAtt(AN_FILLVALUE, ncFloat, value);
		else if (type == ncDouble) putAtt(AN_FILLVALUE, ncDouble, value);
		else return false;
		return true;
	}
	
};

class cSampleVar : public cGeophysicsVar{

private:

public:

	cSampleVar(const NcVar& var) : cGeophysicsVar(var) {	}

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

};

class cLineVar : public cGeophysicsVar{

private:

public:

	cLineVar(const NcVar& var) : cGeophysicsVar(var) {	}

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
			throw(std::runtime_error(msg));
		}

		putVar(vals.data());
		return true;
	}
};

class cGeophysicsNcFile : public NcFile {

private:
	const std::string classname = "cGeophysicsNcFile";	
	NcDim dim_sample;
	NcDim dim_line;	
	std::vector<size_t> line_index_start;
	std::vector<size_t> line_index_count;	
	std::vector<size_t> flight_number;
	std::vector<size_t> line_number;
	   
	bool InitialiseExisting(){
		bool status;
		status = readLineIndex();
		status = getLineNumbers(line_number);
		status = getFlightNumbers(flight_number);
		return true;
	}
	
	bool readLineIndex(){

		bool status;
		size_t ns = getDim(DN_POINT).getSize();
		status = getVarAll(VN_LI_START, line_index_start);
		//line_index_start.pop_back();

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

		}
		else if (filemode == NcFile::replace){

		}
		else if (filemode == NcFile::newFile){

		}
		else{

		}
	};

	//Create new file
	cGeophysicsNcFile(const std::string& ncpath, const std::vector<size_t>& linenumbers, const std::vector<size_t>& nsamples)
		: NcFile(ncpath, FileMode::newFile)
	{
		InitialiseNew(linenumbers,nsamples);
	};

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

		dim_sample = addDim(DN_POINT, nsamples);
		dim_line   = addDim(DN_LINE, nl);
				
		cLineVar vstart = addLineVar(VN_LI_START, ncInt);
		vstart.putVar(line_index_start.data());
		vstart.add_standard_name(VN_LI_START);
		vstart.add_description("zero based index of the first sample in the line");
		vstart.add_units("1");
		
		
		cLineVar vcount = addLineVar(VN_LI_COUNT, ncInt);
		vcount.putVar(line_index_count.data());		
		vcount.add_standard_name(VN_LI_COUNT);
		vcount.add_description("number of samples in the line");
		vcount.add_units("1");

		cLineVar vline = addLineVar(DN_LINE, ncInt);
		vline.putVar(line_number.data());
		vline.add_standard_name(SN_LINE_NUMBER);
		vline.add_description("flight line number");
		vline.add_units("1");

		std::vector<int> sample = increment((int)ntotalsamples(), 0, 1);
		cSampleVar vsample = addSampleVar(DN_POINT, ncInt);
		vsample.putVar(sample.data());
		vsample.add_standard_name(SN_SAMPLE_NUMBER);
		vsample.add_description("sequential point number");
		vsample.add_units("1");
		return true;
	}

	size_t nlines(){ return line_index_start.size(); }
	size_t ntotalsamples(){ return sum(line_index_count); }
	size_t nlinesamples(const size_t lineindex){ return line_index_count[lineindex]; }

	bool getNcVarByAttribute(const std::string& attribute_name, const std::string& attribute_value, NcVar& var){

		std::multimap<std::string, NcVar> vars = getVars();
		for (auto vit = vars.begin(); vit != vars.end(); vit++){
			std::map<std::string, NcVarAtt> atts = vit->second.getAtts();
			for (auto ait = atts.begin(); ait != atts.end(); ait++){
				std::string attname = ait->second.getName();
				if (attname == attribute_name){
					std::string attvalue;
					ait->second.getValues(attvalue);
					if (attvalue == attribute_value){
						var = vit->second;
						return true;
					}
				}
			}
		}
		return false;
	}
	
	size_t getLineIndex(const int& linenumber){
		auto it = std::find(line_number.begin(), line_number.end(), linenumber);				
		return it - line_number.begin();		
	}

	template<typename T>
	bool getLineNumbers(std::vector<T>& vals){
		NcVar var;
		bool status = getNcVarByAttribute(AN_STANDARD_NAME, SN_LINE_NUMBER, var);
		if (status) return getVarAll(var.getName(), vals);
		return false;
	}

	template<typename T>
	bool getFlightNumbers(std::vector<T>& vals){
		NcVar var;
		bool status = getNcVarByAttribute(AN_STANDARD_NAME, SN_FLIGHT_NUMBER, var);
		if (status) return getVarAll(var.getName(), vals);
		return false;
	}
	
	template<typename T>
	bool getVarAll(const std::string& varname, std::vector<T>& vals){
		NcVar var = NcFile::getVar(varname);
		std::vector<NcDim> dims = var.getDims();
		size_t ne = 1;
		for (size_t i = 0; i < dims.size(); i++) ne *= dims[i].getSize();
		vals.resize(ne);
		var.getVar(vals.data());
		return true;
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
		return cSampleVar(getVar(name));
	}

	cLineVar getLineVar(const std::string& name){
		if (getVarCount() == 0) return cLineVar(NcVar());
		return cLineVar(getVar(name));		
	}
	
	cSampleVar addSampleVar(const std::string& name, const NcType& type, const std::vector<NcDim>& dims){

		cSampleVar var = getSampleVar(name);
		if (var.isNull()){
			std::vector<NcDim> vardims = { dim_sample };
			for (size_t i = 0; i<dims.size(); i++){
				vardims.push_back(dims[i]);
			}
			var = addVar(name, type, vardims);
			var.set_parent(this);
			var.set_default_fillvalue();
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
			std::vector<NcDim> vardims = { dim_line };
			for (size_t i = 0; i<dims.size(); i++){
				vardims.push_back(dims[i]);
			}
			var = addVar(name, type, vardims);
			var.set_parent(this);
			var.set_default_fillvalue();
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

	bool addLineStartEndPoints(const std::vector<double>& x1, const std::vector<double>& x2, const std::vector<double>& y1, const std::vector<double>& y2){
		cLineVar vx1 = addLineVar("longitude_first", ncDouble);
		vx1.add_fillvalue(vx1.preferred_double_fillvalue());
		vx1.add_standard_name("longitude_first");
		vx1.add_description("first non-null longitude coordinate in the line");
		vx1.add_units("degree_east");
		vx1.putAll(x1);

		cLineVar vx2 = addLineVar("longitude_last", ncDouble);
		vx2.add_fillvalue(vx2.preferred_double_fillvalue());
		vx2.add_standard_name("longitude_last");
		vx2.add_description("last non-null longitude coordinate in the line");
		vx2.add_units("degree_east");
		vx2.putAll(x2);

		cLineVar vy1 = addLineVar("latitude_first", ncDouble);
		vy1.add_fillvalue(vy1.preferred_double_fillvalue());
		vy1.add_standard_name("latitude_first");
		vy1.add_description("first non-null latitude coordinate in the line");
		vy1.add_units("degree_north");
		vy1.putAll(y1);

		cLineVar vy2 = addLineVar("latitude_last", ncDouble);
		vy2.add_fillvalue(vy2.preferred_double_fillvalue());
		vy2.add_standard_name("latitude_last");
		vy2.add_description("last non-null latitude coordinate in the line");
		vy2.add_units("degree_north");
		vy2.putAll(y2);
		return true;
	}

	bool addCRS(const cCRS& crs){					
		NcVar varcrs = addVar("crs", ncInt);
		varcrs.putAtt("grid_mapping_name", "latitude_longitude");
		varcrs.putAtt("epsg_code", crs.epsg_code.c_str());
		varcrs.putAtt("semi_major_axis", ncDouble, crs.semi_major_axis);
		varcrs.putAtt("inverse_flattening", ncDouble, crs.inverse_flattening);		
		return true;
	}

};

#endif

