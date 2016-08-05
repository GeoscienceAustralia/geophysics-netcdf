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

class cGeophysicsNcFile : public NcFile{

private:

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
		
		dim_sample = addDim(dim_name_sample, nsamples);
		dim_line   = addDim(dim_name_line, nl);
		NcDim lineindexdim = addDim("line_index", nl+1);

		size_t nwindows  = 45;
		NcDim dim_window = addDim("window", nwindows);

		std::vector<NcDim> ldims = { dim_line };
		NcVar lvar = addVar("line", ncInt, ldims);
		lvar.putVar(line_number.data());

		NcVar svar = addVar("line_index_start", ncInt, ldims);
		svar.putVar(line_index_start.data());
		
		NcVar cvar = addVar("line_index_count", ncInt, ldims);
		cvar.putVar(line_index_count.data());

		std::vector<NcDim> sdims = { dim_sample };
		std::vector<double> fid(nsamples);		
		for (size_t i = 0; i < fid.size(); i++) fid[i] = i;
		addSampleVarAll("fiducial", fid);
		fid *= 2;
		addSampleVarAll("fiducial2", fid);
		addSampleVar("fiducial3",ncDouble, 2);
		

		std::vector<NcDim>  emdims = { dim_sample, dim_window };
		std::vector<double> em(nsamples*nwindows);
		for (size_t i = 0; i < em.size(); i++) em[i] = i;

		NcVar v2 = addVar("em_z", ncDouble, emdims);
		v2.putVar(em.data());

		return true;
	}

	NcType nctype(short dummy){ return ncShort; }
	NcType nctype(int dummy){ return ncInt; }
	NcType nctype(double dummy){ return ncDouble; }

	bool addSampleVar(const std::string& name, const NcType& type, const size_t nbands = 1){

		NcVar var = getVar(name);
		if (var.isNull()){						
			NcDim dim = addDim("fred", nbands);			
			std::vector<NcDim> dims = { dim_sample, dim };
			var = addVar(name, type, dims);
			return true;
		}
		return false;
	}

	template<typename T>
	bool addSampleVarAll(const std::string& name, std::vector<T> vals){
		
		NcVar var = getVar(name);
		if (var.isNull()){
			NcType type(vals[0]);
			std::vector<NcDim> dims = { dim_sample };
			var = addVar(name, nctype(vals[0]), dims);
		}
		var.putVar(vals.data());
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
		size_t ns = getDim(dim_name_sample).getSize();
		status = getVarAll(var_name_line_index_start, line_index_start);
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

	const std::string dim_name_sample = "sample";
	const std::string dim_name_line = "line";
	const std::string var_name_line_index_start = "index_lines";
	const std::string att_name_standard_name = "my_standard_name";
	const std::string att_name_line_number = "line_number";
	const std::string att_name_flight_number = "flight_number";

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
		bool status = getNcVarByAttribute(att_name_standard_name, att_name_line_number, var);
		if (status) return getVarAll(var.getName(), vals);
		return false;
	}

	template<typename T>
	bool getFlightNumbers(std::vector<T>& vals){
		NcVar var;
		bool status = getNcVarByAttribute(att_name_standard_name, att_name_flight_number, var);
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
	
};


#endif
