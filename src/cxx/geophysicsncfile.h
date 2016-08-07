/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2016.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

//#include <cstdio>
//#include <vector>
//#include <limits>
#include <netcdf>


#ifndef _geophysicsncfile_H
#define _geophysicsncfile_H

std::string errormsg(const char* file, const int& linenumber, const char* function){
	char s[500];	
	sprintf(s,"Error: %s (line %d) in %s()",file,linenumber,function);
	std::string msg(s);
	return msg;
}

class cGeophysicsNcFile : public NcFile{

private:
	const std::string classname = "cGeophysicsNcFile";
	const std::string DN_SAMPLE = "sample";
	const std::string DN_LINE   = "line";

	const std::string VN_LI_START = "_index_lines";
	const std::string VN_LI_COUNT = "_index_count";

	const std::string AN_STANDARD_NAME = "standard_name";
	const std::string AN_LINE_NUMBER   = "line_number";
	const std::string AN_FLIGHT_NUMBER = "flight_number";

	NcDim dim_sample;
	NcDim dim_line;	
	std::vector<size_t> line_index_start;
	std::vector<size_t> line_index_count;	
	std::vector<size_t> flight_number;
	std::vector<size_t> line_number;

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
		
		dim_sample = addDim(DN_SAMPLE, nsamples);		
		dim_line   = addDim(DN_LINE, nl);				
		
		NcVar vstart = addVar(VN_LI_START, ncInt, dim_line);
		vstart.putVar(line_index_start.data());
		
		NcVar vcount = addVar(VN_LI_COUNT, ncInt, dim_line);
		vcount.putVar(line_index_count.data());
		
		NcVar vline = addVar(DN_LINE, ncInt, dim_line);
		vline.putVar(line_number.data());

		std::vector<int> sample = increment((int)ntotalsamples(),0,1);		
		NcVar vsample = addVar(DN_SAMPLE, ncInt, dim_sample);
		vsample.putVar(sample.data());
		return true;
	}

	bool InitialiseExisting(){
		bool status;
		status = readLineIndex();
		status = getLineNumbers(line_number);
		status = getFlightNumbers(flight_number);
		return true;
	}
	
	bool readLineIndex(){

		bool status;
		size_t ns = getDim(DN_SAMPLE).getSize();
		status = getVarAll(VN_LI_START, line_index_start);
		line_index_start.pop_back();

		size_t nl = nlines();
		line_index_count.resize(nl);
		for (size_t i = 0; i < nl - 1; i++){
			line_index_count[i] = line_index_start[i + 1] - line_index_start[i];
		}
		line_index_count[nl - 1] = ns - line_index_start[nl - 1];
		return true;
	}

public:

	

	static NcType nctype(const short dummy){ return ncShort; }
	static NcType nctype(const int dummy){ return ncInt; }
	static NcType nctype(const double dummy){ return ncDouble; }
	static NcType nctype(const std::string dummy){ return ncString; }

	//Open existing file
	cGeophysicsNcFile(const std::string& ncpath, const FileMode& filemode)
		: NcFile(ncpath, filemode)
	{
		InitialiseExisting();		
	};

	//Create new file
	cGeophysicsNcFile(const std::string& ncpath, const std::vector<size_t>& linenumbers, const std::vector<size_t>& nsamples)
		: NcFile(ncpath, FileMode::newFile)
	{
		InitialiseNew(linenumbers,nsamples);
	};

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
		bool status = getNcVarByAttribute(AN_STANDARD_NAME, AN_LINE_NUMBER, var);
		if (status) return getVarAll(var.getName(), vals);
		return false;
	}

	template<typename T>
	bool getFlightNumbers(std::vector<T>& vals){
		NcVar var;
		bool status = getNcVarByAttribute(AN_STANDARD_NAME, AN_FLIGHT_NUMBER, var);
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
		vals.resize(dims[0].getSize());
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
	
	size_t getVarLength(const NcVar& var){
		std::vector<NcDim> dims = var.getDims();
		size_t len = dims[0].getSize();
		for (size_t di = 1; di < dims.size(); di++){
			len *= dims[di].getSize();
		}
		return len;
	}

	size_t getVarBands(const NcVar& var){
		std::vector<NcDim> dims = var.getDims();
		if (dims.size() == 1) return 1;
		else if (dims.size() == 2) dims[1].getSize();
		else return -1;
	}

	template<typename T>
	NcDim addDimAndVar(const std::string& dimname, const std::vector<T> dimvals){

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

	NcVar addSampleVar(const std::string& name, const NcType& type, const NcDim& banddim = NcDim()){

		NcVar var = getVar(name);
		if (var.isNull()){
			std::vector<NcDim> vardims = { dim_sample };
			if (!banddim.isNull()){				
				vardims.push_back(banddim);
			}
			var = addVar(name, type, vardims);
		}		
		return var;
	}

	NcVar addSampleVar(const std::string& name, const NcType& type, const std::vector<NcDim>& dims){

		NcVar var = getVar(name);
		if (var.isNull()){
			std::vector<NcDim> vardims = { dim_sample };
			for(size_t i=0; i<dims.size(); i++){
				vardims.push_back(dims[i]);
			}
			var = addVar(name, type, vardims);
		}
		return var;
	}

	template<typename T>
	bool putSampleVarAll(const NcVar& var, std::vector<T> vals){
		
		if (var.isNull()){
			std::string msg = strprint("Attempt to write to a Null variable\n");
			msg += errormsg(__FILE__, __LINE__, __FUNCTION__);
			throw(std::runtime_error(msg));
		}
		
		size_t len = getVarLength(var);				
		if (vals.size() != len){
			std::string msg = strprint("Attempt to write variable (%s) with non-matching size\n", var.getName().c_str());
			msg += errormsg(__FILE__, __LINE__, __FUNCTION__);
			throw(std::runtime_error(msg));
		}	

		var.putVar(vals.data());
		return true;		
	}

	template<typename T>
	bool putSampleVarLine(const NcVar& var, std::vector<T> vals, size_t lineindex, size_t bandindex=0){		
		if (var.isNull()){ return false; }
		std::vector<size_t> startp = { line_index_start[lineindex], bandindex};
		std::vector<size_t> countp = { line_index_count[lineindex], 1 };
		var.putVar(startp,countp,vals.data());
		return true;		
	}
};


#endif
