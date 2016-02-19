/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
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

	std::vector<size_t> line_index_start;
	std::vector<size_t> line_index_count;	
	std::vector<size_t> flight_number;
	std::vector<size_t> line_number;

	bool readLineIndex(){

		bool status;
		size_t ns = getDim(dim_name_sample).getSize();
		status = getVar(var_name_line_index_start, line_index_start);
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

	cGeophysicsNcFile(const std::string& ncpath, const FileMode& filemode)
		: NcFile(ncpath, filemode)
	{
		readLineIndex();
	};

	size_t nlines(){ return line_index_start.size(); }

	bool getVarByAttribute(const std::string& attribute_name, const std::string& attribute_value, NcVar& var){

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

	template<typename T>
	bool getLineIndex(const int& linenumber){

	}

	template<typename T>
	bool getFlightNumbers(std::vector<T>& vals){
		NcVar var;
		bool status = getVarByAttribute(att_name_standard_name, att_name_flight_number, var);
		if (status) return getVar(var.getName(), vals);
		return false;
	}

	template<typename T>
	bool getVar(const std::string& varname, std::vector<T>& vals){
		NcVar var = NcFile::getVar(varname);
		std::vector<NcDim> dims = var.getDims();
		size_t ne = 1;
		for (size_t i = 0; i < dims.size(); i++) ne *= dims[i].getSize();
		vals.resize(ne);
		var.getVar(vals.data());
		return true;
	}

	template<typename T>
	bool getVarByLineIndex(const std::string& varname, std::vector<T>& vals, const size_t& lineindex){
		NcVar var = NcFile::getVar(varname);
		std::vector<NcDim> dims = var.getDims();

		std::vector<size_t> start(1);
		std::vector<size_t> count(1);
		start[0] = line_index_start[lineindex];
		count[0] = line_index_count[lineindex];

		size_t ne = count[0];
		for (size_t i = 1; i < dims.size(); i++){
			ne *= dims[i].getSize();
		}
		vals.resize(ne);
		var.getVar(start, count, vals.data());
		return true;
	}

	template<typename T>
	bool getVarByLineNumber(const std::string& varname, std::vector<T>& vals, const size_t& linenumber){
		
	}
};


#endif
